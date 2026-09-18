!>@brief Cross-solver (FUSS <-> MOSE) MPI exchange for type(103) coupling BCs.
!>
!> For each BC of type(103 the "primary" block (Bm) belongs to one solver
!> while the "secondary" block (Bs) belongs to the OTHER solver. When the
!> two blocks are owned by different MPI ranks (in their respective solver
!> partitions), the primary-side rank needs data about the secondary-side
!> cell (Is,Js,Ks) and its surrounding coupling stencil.
!>
!> Stencil layout (matches Tg(6) and Pg(1,6) in Advanced_Types_m):
!>   indices 1..gc_       : normal-direction ghost cells (gc_=2 by default)
!>   indices gc_+1..gc_+4 : 4 tangential corner cells
!>   Total NSTENC = gc_+4 = 6 entries.
!>
!> Payloads bundled per BC entry:
!>   - FUSS recv (data from MOSE-side Bs):
!>       Tg(1..NSTENC)        [NSTENC reals]
!>       P(1:nsc, Is,Js,Ks)   [nsc reals, needed by FUSS wall-T routine]
!>   - MOSE recv (data from FUSS-side Bs):
!>       Pg(1, 1..NSTENC)     [NSTENC reals]
!>       matID(Is,Js,Ks)      [1 real, needed by MOSE wall-T routine]
!>
!> The schedule is built once per MG level (both solvers share BC ordering
!> because both read the same hydra inifile on every rank). Exchanges are
!> aggregated per remote rank with non-blocking point-to-point MPI.
!>
!> In serial mode (USE_MPI not defined), all routines are cheap no-ops; the
!> existing Fill_Coupled_Ghost_Cell and wall-T routines continue to work
!> because every BC is local-local.
module hydraMF_Mod_CouplingExchange
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I4 => int32
  use FUSS_Advanced_Types_m, only: FUSS_domain_type
  use MOSE_Advanced_Types_m, only: MOSE_domain_type
  use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner, &
                          FUSS_is_local_block => is_local_block
  use MOSE_Mod_MPI, only: MOSE_block_owner => block_owner, &
                          MOSE_is_local_block => is_local_block, &
                          mpi_rank_, mpi_size_, check_mpi_error
  use hydra_Global_m, only: gc_

  implicit none
  private

  public :: build_coupling_schedule
  public :: cleanup_coupling_schedule
  public :: exchange_coupling_Tg_P
  public :: exchange_coupling_Pg_matID
  public :: exchange_coupling_reverse_flux_FUSS
  public :: exchange_coupling_reverse_flux_MOSE
  public :: partner_of_FUSS_bc, partner_of_MOSE_bc
  public :: is_FUSS_bc_local_local, is_MOSE_bc_local_local
  public :: get_F_remote_P, get_M_remote_matID
  public :: coupled_secondary_node_indices
  public :: exchange_coupling_metric_nodes, free_coupling_metric_nodes
  public :: own103_FUSS, own103_MOSE
  public :: n_F_meta_entries, n_M_meta_entries
  public :: get_F_meta_entry, get_M_meta_entry

  ! Coupling stencil size: gc_ ghost cells (normal) + 4 corner cells (tangential)
  ! Matches Tg(6) and Pg(1,6) in Advanced_Types_m for the default gc_=2.
  integer, parameter, public :: NSTENC = 6

  ! MPI message tags (distinct per direction + quantity)
  integer, parameter :: TAG_TG_P    = 7001
  integer, parameter :: TAG_PG_MID  = 7002
  integer, parameter :: TAG_RFLUX_F = 7003  ! reverse flux into MOSE%ext_flux(np)
  integer, parameter :: TAG_RFLUX_M = 7004  ! reverse flux into FUSS%ext_flux
  integer, parameter :: TAG_META_F   = 7005  ! secondary node coords, MOSE -> FUSS
  integer, parameter :: TAG_META_M   = 7006  ! secondary node coords, FUSS -> MOSE

  !> Node coordinates a primary rank needs from the secondary block to build the
  !> coupling metric tensor: 4 nodes per ghost level, 3 components each.
  integer, parameter :: NMETA_NODE = 4


  !> Entry on the primary-side ("recv") rank: local BC needing remote stencil data.
  type :: recv_entry_t
    integer :: bc_idx          !< local BC index (into *_domain%bc(:))
    integer :: remote_rank     !< rank owning the secondary block (other solver)
    integer :: remote_bs       !< secondary block id (global, in other solver)
  end type recv_entry_t

  !> Entry on the secondary-side ("send") rank: local Bs block serving a remote BC.
  !> Holds the metadata needed to compute the stencil locally (no blkm access).
  type :: send_entry_t
    integer :: global_bc_idx   !< index into the "primary" solver's BC array (for matching)
    integer :: remote_rank     !< rank owning the primary block (the destination)
    integer :: local_bs        !< secondary block id (local to this rank, in this solver)
    integer :: is, js, ks, fs
    integer :: d11, d12, d21, d22
  end type send_entry_t

  !> Per-rank aggregation (contiguous run of entries sharing the same remote rank).
  type :: rank_group_t
    integer :: rank   = -1
    integer :: offset = 0      !< 1-based start in entry list
    integer :: count  = 0
  end type rank_group_t

  !> One direction of a paired exchange (send-side OR recv-side view of the same flow).
  type :: dir_t
    integer :: n_entries = 0
    integer :: payload_per_entry = 0      !< reals per entry in this direction
    type(rank_group_t), allocatable :: groups(:)
    integer :: n_groups = 0
    real(R8), allocatable :: buf(:)
    integer,  allocatable :: reqs(:)
  end type dir_t

  type :: level_sched_t
    ! FUSS-primary BCs: FUSS owns Bm, needs data about MOSE Bs.
    type(dir_t) :: F_recv                 !< recv on FUSS-owner rank (Tg + P)
    type(recv_entry_t), allocatable :: F_recv_ent(:)
    type(dir_t) :: F_send                 !< send on MOSE-owner rank (pack from MOSE block)
    type(send_entry_t), allocatable :: F_send_ent(:)

    ! MOSE-primary BCs: MOSE owns Bm, needs data about FUSS Bs.
    type(dir_t) :: M_recv                 !< recv on MOSE-owner rank (Pg + matID)
    type(recv_entry_t), allocatable :: M_recv_ent(:)
    type(dir_t) :: M_send                 !< send on FUSS-owner rank (pack from FUSS block)
    type(send_entry_t), allocatable :: M_send_ent(:)

    ! Remote-received single-cell data, stored in BC-index-sized arrays so the
    ! wall-T routine can reuse them without touching the remote block directly.
    ! Allocated only on the FUSS-owner (for P_cell) and MOSE-owner (for matID_cell).
    real(R8), allocatable :: F_remote_P(:,:)     !< (nsc, FUSS%nbound), for type(103te only
    real(R8), allocatable :: M_remote_matID(:)   !< (MOSE%nbound),       for type(103remote only

    ! Partner indices: for local BCs whose Bs is ALSO locally owned (same rank owns both
    ! the FUSS Bm and the partner MOSE Bm'), the reverse-flux update into partner%ext_flux
    ! can be done locally. Store the partner BC index; 0 means partner is remote (MPI path).
    integer, allocatable :: partner_of_FUSS(:)   !< (FUSS%nbound) -> MOSE bc idx (or 0)
    integer, allocatable :: partner_of_MOSE(:)   !< (MOSE%nbound) -> FUSS bc idx (or 0)

    ! Compact lists of the type(103) BC indices whose PRIMARY block this rank
    ! owns -- i.e. exactly the coupling work this rank has to do.
    !
    ! The coupling loops used to sweep the whole bc array on every rank and
    ! `cycle` past everything they did not own. That array is global, and it
    ! GROWS with the block count: decomposing the 6.9M-cell benchmark from 1
    ! block to 80 takes it from 658,080 entries to 2,475,360, because every cut
    ! manufactures type(101 split faces. The type(103 work inside it is constant
    ! at 115,200 -- the interface is the same surface however the mesh is cut.
    ! So adding ranks lengthened the sweep 3.8x while shrinking each rank's share
    ! of the useful work, and the coupling got slower with rank count instead of
    ! faster. Threads did help, because the !$omp do was on the sweep itself.
    ! Iterating these lists makes the cost proportional to what the rank owns.
    integer, allocatable :: own103_F(:)          !< FUSS bc indices, primary block local
    integer, allocatable :: own103_M(:)          !< MOSE bc indices, primary block local
    integer :: n_own103_F = 0                    !< used length of own103_F
    integer :: n_own103_M = 0                    !< used length of own103_M

    ! Per-entry reverse-flux scratch, in bc-index space because
    ! exchange_coupling_reverse_flux_* addresses it by bc index. Persistent: it
    ! used to be allocated and zeroed at global size on every call (14 MB for
    ! MOSE at 80 blocks, four times per coupling step, >95% of it untouched).
    ! Every index the exchange reads is written by the same call, so it needs no
    ! clearing between steps.
    real(R8), allocatable :: rflux_F(:)          !< (FUSS%nbound)
    real(R8), allocatable :: rflux_M(:)          !< (MOSE%nbound)

    ! Remote reverse-flux exchange (same pairing as forward, but reversed direction and
    ! using a scalar accumulator per entry).
    type(dir_t) :: RF_FUSS_send   !< FUSS-owner sends -0.5*F_new_flux to MOSE-owner
    integer, allocatable :: RF_FUSS_send_peer(:)  !< peer MOSE bc idx (for matching); not always needed
    type(dir_t) :: RF_FUSS_recv   !< MOSE-owner receives (applied to MOSE partner ext_flux(np))
    integer, allocatable :: RF_FUSS_recv_mosebc(:) !< maps recv slot -> MOSE bc idx

    type(dir_t) :: RF_MOSE_send   !< MOSE-owner sends -0.5*M_new_flux(np) to FUSS-owner
    type(dir_t) :: RF_MOSE_recv   !< FUSS-owner receives (applied to FUSS partner ext_flux)
    integer, allocatable :: RF_MOSE_recv_fussbc(:) !< maps recv slot -> FUSS bc idx

    ! One-off exchange of the secondary-block node coordinates that
    ! BC_Coupled_Metrics_* needs to build Mg. Same rank pairing as the forward
    ! directions, so the groups are copied from them. Allocated during setup
    ! and released again as soon as the metrics have been built.
    type(dir_t) :: F_meta_recv, F_meta_send
    type(dir_t) :: M_meta_recv, M_meta_send

    logical :: built = .false.
  end type level_sched_t

  type(level_sched_t), allocatable, target :: levels(:)

contains


  !---------------------------------------------------------------------------
  !> Indices, inside the SECONDARY block, of the four nodes that feed ghost
  !> level g of a type(103) coupling BC.
  !>
  !> They depend only on secondary-side data (Is,Js,Ks,Fs and the four
  !> orientation flags d11s..d22s) and never on the primary block, which is what
  !> lets the rank owning the secondary block pack these coordinates for a
  !> primary block it knows nothing about.  BC_Coupled_Metrics_* and the
  !> metric-node packers all call this, so they cannot drift apart.
  !>
  !> Node order: 1 = the connected node itself, 2 = its (i1-1,i2) neighbour,
  !> 3 = (i1,i2-1), 4 = (i1-1,i2-1).
  pure subroutine coupled_secondary_node_indices(g, Is, Js, Ks, Fs, &
                                                 d11s, d12s, d21s, d22s, idx)
    implicit none
    integer, intent(in)  :: g, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    integer, intent(out) :: idx(3,NMETA_NODE)
    ! Local
    integer :: Is1, Js1, Ks1, II, JJ, KK
    integer :: guides(6,3), guides2(6,3)

    guides(1,1)= 0 ; guides(1,2)= 0 ; guides(1,3)= 0 ; guides(2,1)=-1 ; guides(2,2)= 0 ; guides(2,3)= 0
    guides(3,1)= 0 ; guides(3,2)= 0 ; guides(3,3)= 0 ; guides(4,1)= 0 ; guides(4,2)=-1 ; guides(4,3)= 0
    guides(5,1)= 0 ; guides(5,2)= 0 ; guides(5,3)= 0 ; guides(6,1)= 0 ; guides(6,2)= 0 ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0 ; guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0 ; guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1 ; guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2

    ! connected cell node index
    Is1 = Is + (gc_-g)*guides(Fs,1) + (g-1)*guides2(Fs,1)
    Js1 = Js + (gc_-g)*guides(Fs,2) + (g-1)*guides2(Fs,2)
    Ks1 = Ks + (gc_-g)*guides(Fs,3) + (g-1)*guides2(Fs,3)

    select case(Fs)
      case(1:2) ! i-faces
        II = Is1
        JJ = (2*Js1+d11s+d21s)/2
        KK = (2*Ks1+d12s+d22s)/2
      case(3:4) ! j-faces
        II = (2*Is1+d11s+d21s)/2
        JJ = Js1
        KK = (2*Ks1+d12s+d22s)/2
      case default ! k-faces
        II = (2*Is1+d11s+d21s)/2
        JJ = (2*Js1+d12s+d22s)/2
        KK = Ks1
    end select
    idx(1,1) = II ; idx(2,1) = JJ ; idx(3,1) = KK

    ! node i1-1,i2
    select case(Fs)
      case(1:2)
        idx(1,2) = II            ; idx(2,2) = JJ-d11s       ; idx(3,2) = KK-d12s
      case(3:4)
        idx(1,2) = II-d11s       ; idx(2,2) = JJ            ; idx(3,2) = KK-d12s
      case default
        idx(1,2) = II-d11s       ; idx(2,2) = JJ-d12s       ; idx(3,2) = KK
    end select

    ! node i1,i2-1
    select case(Fs)
      case(1:2)
        idx(1,3) = II            ; idx(2,3) = JJ-d21s       ; idx(3,3) = KK-d22s
      case(3:4)
        idx(1,3) = II-d21s       ; idx(2,3) = JJ            ; idx(3,3) = KK-d22s
      case default
        idx(1,3) = II-d21s       ; idx(2,3) = JJ-d22s       ; idx(3,3) = KK
    end select

    ! node i1-1,i2-1
    select case(Fs)
      case(1:2)
        idx(1,4) = II            ; idx(2,4) = JJ-d11s-d21s  ; idx(3,4) = KK-d12s-d22s
      case(3:4)
        idx(1,4) = II-d11s-d21s  ; idx(2,4) = JJ            ; idx(3,4) = KK-d12s-d22s
      case default
        idx(1,4) = II-d11s-d21s  ; idx(2,4) = JJ-d12s-d22s  ; idx(3,4) = KK
    end select

  end subroutine coupled_secondary_node_indices


  !---------------------------------------------------------------------------
  !> Build the cross-solver coupling schedule for one MG level.
  !> Idempotent on built=.true.  Depends on the two solvers' block_owner arrays
  !> having already been populated by partition_blocks().
  !> Type-103 interface entries whose solid block belongs to this rank.
  integer function count_owned_interface(FUSS_domain) result(n)
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    integer :: i

    n = 0
    do i = 1, FUSS_domain%nbound
      if (FUSS_domain%bc(i)%type /= 103) cycle
      if (FUSS_block_owner(FUSS_domain%bc(i)%b) /= mpi_rank_) cycle
      n = n + 1
    end do
  end function count_owned_interface


  !> Collect the type(103) BC indices whose primary block this rank owns, for
  !> each phase, plus the persistent reverse-flux scratch addressed by bc index.
  !>
  !> This is the only place that still walks the full bc array, and it runs once
  !> per multigrid level at setup rather than four times per coupling step.
  subroutine build_own103_lists(FUSS_domain, MOSE_domain, L)
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(level_sched_t), intent(inout) :: L
    integer :: i, n

    ! --- FUSS ---
    n = 0
    do i = 1, FUSS_domain%nbound
      if (FUSS_domain%bc(i)%type /= 103) cycle
      if (.not. owns_FUSS_primary(FUSS_domain, i)) cycle
      n = n + 1
    end do
    allocate(L%own103_F(max(n,1))); L%own103_F = 0
    L%n_own103_F = n
    n = 0
    do i = 1, FUSS_domain%nbound
      if (FUSS_domain%bc(i)%type /= 103) cycle
      if (.not. owns_FUSS_primary(FUSS_domain, i)) cycle
      n = n + 1
      L%own103_F(n) = i
    end do

    ! --- MOSE ---
    n = 0
    do i = 1, MOSE_domain%nbound
      if (MOSE_domain%bc(i)%type /= 103) cycle
      if (.not. owns_MOSE_primary(MOSE_domain, i)) cycle
      n = n + 1
    end do
    allocate(L%own103_M(max(n,1))); L%own103_M = 0
    L%n_own103_M = n
    n = 0
    do i = 1, MOSE_domain%nbound
      if (MOSE_domain%bc(i)%type /= 103) cycle
      if (.not. owns_MOSE_primary(MOSE_domain, i)) cycle
      n = n + 1
      L%own103_M(n) = i
    end do

    allocate(L%rflux_F(max(FUSS_domain%nbound,1))); L%rflux_F = 0.0d0
    allocate(L%rflux_M(max(MOSE_domain%nbound,1))); L%rflux_M = 0.0d0
  end subroutine build_own103_lists


  !> Ownership of a BC's primary block. With no partition (serial build, or
  !> before partition_blocks has run) every block is this rank's.
  logical function owns_FUSS_primary(FUSS_domain, i) result(own)
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    integer, intent(in) :: i
    own = .true.
    if (.not. allocated(FUSS_block_owner)) return
    own = ( FUSS_block_owner(FUSS_domain%bc(i)%b) == mpi_rank_ )
  end function owns_FUSS_primary

  logical function owns_MOSE_primary(MOSE_domain, i) result(own)
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in) :: i
    own = .true.
    if (.not. allocated(MOSE_block_owner)) return
    own = ( MOSE_block_owner(MOSE_domain%bc(i)%b) == mpi_rank_ )
  end function owns_MOSE_primary


  !> The FUSS coupling work this rank owns: count, bc indices, and the
  !> reverse-flux scratch. `n` is 0 and the pointers are unassociated when the
  !> schedule for `level` has not been built.
  subroutine own103_FUSS(level, n, idx, rflux)
    implicit none
    integer, intent(in) :: level
    integer, intent(out) :: n
    integer, pointer, intent(out) :: idx(:)
    real(R8), pointer, intent(out), optional :: rflux(:)
    type(level_sched_t), pointer :: L
    n = 0
    nullify(idx)
    if (present(rflux)) nullify(rflux)
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    L => levels(level)
    if (.not. allocated(L%own103_F)) return
    n = L%n_own103_F
    idx => L%own103_F
    if (present(rflux) .and. allocated(L%rflux_F)) rflux => L%rflux_F
  end subroutine own103_FUSS

  subroutine own103_MOSE(level, n, idx, rflux)
    implicit none
    integer, intent(in) :: level
    integer, intent(out) :: n
    integer, pointer, intent(out) :: idx(:)
    real(R8), pointer, intent(out), optional :: rflux(:)
    type(level_sched_t), pointer :: L
    n = 0
    nullify(idx)
    if (present(rflux)) nullify(rflux)
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    L => levels(level)
    if (.not. allocated(L%own103_M)) return
    n = L%n_own103_M
    idx => L%own103_M
    if (present(rflux) .and. allocated(L%rflux_M)) rflux => L%rflux_M
  end subroutine own103_MOSE


  subroutine build_coupling_schedule(FUSS_domain, MOSE_domain, level)
    use MOSE_Global_m, only: nsc
    use hydra_Mod_Timers, only: timer_set_interface_counts
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in) :: level
    ! Local
    integer :: i, Bm, Bs, fuss_o, mose_o
    integer :: n_F_recv, n_F_send, n_M_recv, n_M_send
    integer :: n_iface_owned
    type(level_sched_t), pointer :: L

    if (.not. allocated(levels)) allocate(levels(max(level,1)))
    if (level > size(levels)) then
      block
        type(level_sched_t), allocatable :: tmp(:)
        integer :: n_old, k
        n_old = size(levels)
        allocate(tmp(level))
        do k = 1, n_old
          tmp(k) = levels(k)
        end do
        call move_alloc(tmp, levels)
      end block
    end if
    L => levels(level)
    if (L%built) return

    call reset_level(L)

    ! Partner index always useful (even in serial); allocate sized to bc arrays.
    allocate(L%partner_of_FUSS(FUSS_domain%nbound)); L%partner_of_FUSS = 0
    allocate(L%partner_of_MOSE(MOSE_domain%nbound)); L%partner_of_MOSE = 0
    call build_partner_indices(FUSS_domain, MOSE_domain, L)

    ! Built before the serial/single-rank early returns below: the coupling loops
    ! iterate these lists in every build, MPI or not.
    call build_own103_lists(FUSS_domain, MOSE_domain, L)

#ifndef USE_MPI
    call timer_set_interface_counts(0, count_owned_interface(FUSS_domain))
    L%built = .true.
    return
#else
    if (mpi_size_ <= 1) then
      call timer_set_interface_counts(0, count_owned_interface(FUSS_domain))
      L%built = .true.
      return
    end if

    ! --- PASS 1: count entries per direction ---
    n_F_recv = 0; n_F_send = 0
    do i = 1, FUSS_domain%nbound
      if (FUSS_domain%bc(i)%type /= 103) cycle
      Bm = FUSS_domain%bc(i)%b
      Bs = FUSS_domain%bc(i)%bs
      fuss_o = FUSS_block_owner(Bm)
      mose_o = MOSE_block_owner(Bs)
      if (fuss_o == mose_o) cycle  ! local-local: no MPI
      if (mpi_rank_ == fuss_o) n_F_recv = n_F_recv + 1
      if (mpi_rank_ == mose_o) n_F_send = n_F_send + 1
    end do

    n_M_recv = 0; n_M_send = 0
    do i = 1, MOSE_domain%nbound
      if (MOSE_domain%bc(i)%type /= 103) cycle
      Bm = MOSE_domain%bc(i)%b
      Bs = MOSE_domain%bc(i)%bs
      mose_o = MOSE_block_owner(Bm)
      fuss_o = FUSS_block_owner(Bs)
      if (mose_o == fuss_o) cycle
      if (mpi_rank_ == mose_o) n_M_recv = n_M_recv + 1
      if (mpi_rank_ == fuss_o) n_M_send = n_M_send + 1
    end do

    ! How much of this rank's share of the interface has to go over MPI. MOSE
    ! and FUSS partition their own phase independently, so nothing keeps a
    ! solid block and the fluid block above it on the same rank; this ratio is
    ! what turns interface work into messages as the rank count grows.
    n_iface_owned = count_owned_interface(FUSS_domain)
    call timer_set_interface_counts(n_F_recv, n_iface_owned)

    L%F_recv%n_entries = n_F_recv
    L%F_send%n_entries = n_F_send
    L%M_recv%n_entries = n_M_recv
    L%M_send%n_entries = n_M_send

    L%F_recv%payload_per_entry = NSTENC + nsc
    L%F_send%payload_per_entry = NSTENC + nsc
    L%M_recv%payload_per_entry = NSTENC + 1
    L%M_send%payload_per_entry = NSTENC + 1

    allocate(L%F_recv_ent(max(n_F_recv,1)))
    allocate(L%F_send_ent(max(n_F_send,1)))
    allocate(L%M_recv_ent(max(n_M_recv,1)))
    allocate(L%M_send_ent(max(n_M_send,1)))

    ! --- PASS 2: fill entries ---
    call fill_FUSS_entries(FUSS_domain, L)
    call fill_MOSE_entries(MOSE_domain, L)

    ! Sort each list by remote_rank and build rank-groups
    call sort_recv_by_rank(L%F_recv_ent, n_F_recv)
    call sort_recv_by_rank(L%M_recv_ent, n_M_recv)
    call sort_send_by_rank(L%F_send_ent, n_F_send)
    call sort_send_by_rank(L%M_send_ent, n_M_send)
    call build_groups_recv(L%F_recv_ent, n_F_recv, L%F_recv)
    call build_groups_recv(L%M_recv_ent, n_M_recv, L%M_recv)
    call build_groups_send(L%F_send_ent, n_F_send, L%F_send)
    call build_groups_send(L%M_send_ent, n_M_send, L%M_send)

    ! Allocate per-direction buffers and request arrays
    call alloc_dir_buf(L%F_recv)
    call alloc_dir_buf(L%F_send)
    call alloc_dir_buf(L%M_recv)
    call alloc_dir_buf(L%M_send)

    ! Per-BC remote storage for the wall-T routines
    if (n_F_recv > 0) then
      allocate(L%F_remote_P(nsc, FUSS_domain%nbound)); L%F_remote_P = 0.0d0
    end if
    if (n_M_recv > 0) then
      allocate(L%M_remote_matID(MOSE_domain%nbound)); L%M_remote_matID = 0.0d0
    end if

    ! Reverse-flux schedules mirror the forward pairings, but the direction is
    ! swapped: forward FUSS_recv becomes backward FUSS->MOSE send, etc.
    call build_reverse_flux_schedule(FUSS_domain, MOSE_domain, L)

    L%built = .true.
#endif
  end subroutine build_coupling_schedule


  !---------------------------------------------------------------------------
  subroutine cleanup_coupling_schedule()
    implicit none
    integer :: k
    if (.not. allocated(levels)) return
    do k = 1, size(levels)
      call reset_level(levels(k))
    end do
    deallocate(levels)
  end subroutine cleanup_coupling_schedule


  !---------------------------------------------------------------------------
  !> Build partner_of_* arrays (purely local/serial bookkeeping).
  !> partner_of_FUSS(i) points to the MOSE bc index j whose (Bm_M,Im_M,Jm_M,Km_M,Fm_M)
  !> equals FUSS's (Bs,Is,Js,Ks,Fs).
  subroutine build_partner_indices(FUSS_domain, MOSE_domain, L)
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(level_sched_t), intent(inout) :: L
    integer :: i, j
    do i = 1, FUSS_domain%nbound
      if (FUSS_domain%bc(i)%type /= 103) cycle
      do j = 1, MOSE_domain%nbound
        if (MOSE_domain%bc(j)%type /= 103) cycle
        if (FUSS_domain%bc(i)%bs == MOSE_domain%bc(j)%b  .and. &
            FUSS_domain%bc(i)%is == MOSE_domain%bc(j)%i  .and. &
            FUSS_domain%bc(i)%js == MOSE_domain%bc(j)%j  .and. &
            FUSS_domain%bc(i)%ks == MOSE_domain%bc(j)%k  .and. &
            FUSS_domain%bc(i)%fs == MOSE_domain%bc(j)%f) then
          L%partner_of_FUSS(i) = j
          L%partner_of_MOSE(j) = i
          exit
        end if
      end do
    end do
  end subroutine build_partner_indices


#ifdef USE_MPI
  !---------------------------------------------------------------------------
  subroutine fill_FUSS_entries(FUSS_domain, L)
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(level_sched_t), intent(inout) :: L
    integer :: i, Bm, Bs, fuss_o, mose_o
    integer :: nr, ns
    nr = 0; ns = 0
    do i = 1, FUSS_domain%nbound
      if (FUSS_domain%bc(i)%type /= 103) cycle
      Bm = FUSS_domain%bc(i)%b
      Bs = FUSS_domain%bc(i)%bs
      fuss_o = FUSS_block_owner(Bm)
      mose_o = MOSE_block_owner(Bs)
      if (fuss_o == mose_o) cycle
      if (mpi_rank_ == fuss_o) then
        nr = nr + 1
        L%F_recv_ent(nr)%bc_idx      = i
        L%F_recv_ent(nr)%remote_rank = mose_o
        L%F_recv_ent(nr)%remote_bs   = Bs
      end if
      if (mpi_rank_ == mose_o) then
        ns = ns + 1
        L%F_send_ent(ns)%global_bc_idx = i
        L%F_send_ent(ns)%remote_rank   = fuss_o
        L%F_send_ent(ns)%local_bs      = Bs
        L%F_send_ent(ns)%is  = FUSS_domain%bc(i)%is
        L%F_send_ent(ns)%js  = FUSS_domain%bc(i)%js
        L%F_send_ent(ns)%ks  = FUSS_domain%bc(i)%ks
        L%F_send_ent(ns)%fs  = FUSS_domain%bc(i)%fs
        L%F_send_ent(ns)%d11 = FUSS_domain%bc(i)%d11
        L%F_send_ent(ns)%d12 = FUSS_domain%bc(i)%d12
        L%F_send_ent(ns)%d21 = FUSS_domain%bc(i)%d21
        L%F_send_ent(ns)%d22 = FUSS_domain%bc(i)%d22
      end if
    end do
  end subroutine fill_FUSS_entries


  !---------------------------------------------------------------------------
  subroutine fill_MOSE_entries(MOSE_domain, L)
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(level_sched_t), intent(inout) :: L
    integer :: i, Bm, Bs, fuss_o, mose_o
    integer :: nr, ns
    nr = 0; ns = 0
    do i = 1, MOSE_domain%nbound
      if (MOSE_domain%bc(i)%type /= 103) cycle
      Bm = MOSE_domain%bc(i)%b
      Bs = MOSE_domain%bc(i)%bs
      mose_o = MOSE_block_owner(Bm)
      fuss_o = FUSS_block_owner(Bs)
      if (mose_o == fuss_o) cycle
      if (mpi_rank_ == mose_o) then
        nr = nr + 1
        L%M_recv_ent(nr)%bc_idx      = i
        L%M_recv_ent(nr)%remote_rank = fuss_o
        L%M_recv_ent(nr)%remote_bs   = Bs
      end if
      if (mpi_rank_ == fuss_o) then
        ns = ns + 1
        L%M_send_ent(ns)%global_bc_idx = i
        L%M_send_ent(ns)%remote_rank   = mose_o
        L%M_send_ent(ns)%local_bs      = Bs
        L%M_send_ent(ns)%is  = MOSE_domain%bc(i)%is
        L%M_send_ent(ns)%js  = MOSE_domain%bc(i)%js
        L%M_send_ent(ns)%ks  = MOSE_domain%bc(i)%ks
        L%M_send_ent(ns)%fs  = MOSE_domain%bc(i)%fs
        L%M_send_ent(ns)%d11 = MOSE_domain%bc(i)%d11
        L%M_send_ent(ns)%d12 = MOSE_domain%bc(i)%d12
        L%M_send_ent(ns)%d21 = MOSE_domain%bc(i)%d21
        L%M_send_ent(ns)%d22 = MOSE_domain%bc(i)%d22
      end if
    end do
  end subroutine fill_MOSE_entries


  !---------------------------------------------------------------------------
  !> Reverse-flux schedule: pairs entries by global_bc_idx matching forward dirs.
  !> FUSS forward recv -> reverse send from FUSS to MOSE; uses same rank pairs.
  subroutine build_reverse_flux_schedule(FUSS_domain, MOSE_domain, L)
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(level_sched_t), intent(inout) :: L
    integer :: i, j, n, k
    integer :: Bm, Bs, fuss_o, mose_o
    integer :: nrf_s, nrf_r, nrm_s, nrm_r

    ! FUSS-primary: forward Tg flowed MOSE->FUSS. Reverse flux flows FUSS->MOSE
    ! and is applied into MOSE bc(partner)%ext_flux(np).
    nrf_s = L%F_recv%n_entries           ! FUSS owner sends reverse flux
    nrf_r = L%F_send%n_entries           ! MOSE owner receives reverse flux
    L%RF_FUSS_send%n_entries = nrf_s
    L%RF_FUSS_send%payload_per_entry = 1
    L%RF_FUSS_recv%n_entries = nrf_r
    L%RF_FUSS_recv%payload_per_entry = 1
    ! The groups mirror the forward direction's rank-group layout (entries already sorted).
    call copy_groups(L%F_recv, L%RF_FUSS_send)
    call copy_groups(L%F_send, L%RF_FUSS_recv)
    call alloc_dir_buf(L%RF_FUSS_send)
    call alloc_dir_buf(L%RF_FUSS_recv)

    ! Map RF_FUSS_recv slot -> local MOSE bc idx (partner of the global FUSS bc idx)
    allocate(L%RF_FUSS_recv_mosebc(max(nrf_r,1)))
    L%RF_FUSS_recv_mosebc = 0
    do k = 1, nrf_r
      i = L%F_send_ent(k)%global_bc_idx
      ! partner: the MOSE bc whose (b,i,j,k,f) matches (bs,is,js,ks,fs) of FUSS bc i
      do j = 1, MOSE_domain%nbound
        if (MOSE_domain%bc(j)%type /= 103) cycle
        if (FUSS_domain%bc(i)%bs == MOSE_domain%bc(j)%b  .and. &
            FUSS_domain%bc(i)%is == MOSE_domain%bc(j)%i  .and. &
            FUSS_domain%bc(i)%js == MOSE_domain%bc(j)%j  .and. &
            FUSS_domain%bc(i)%ks == MOSE_domain%bc(j)%k  .and. &
            FUSS_domain%bc(i)%fs == MOSE_domain%bc(j)%f) then
          L%RF_FUSS_recv_mosebc(k) = j
          exit
        end if
      end do
    end do

    ! MOSE-primary: forward Pg flowed FUSS->MOSE. Reverse flux flows MOSE->FUSS.
    nrm_s = L%M_recv%n_entries
    nrm_r = L%M_send%n_entries
    L%RF_MOSE_send%n_entries = nrm_s
    L%RF_MOSE_send%payload_per_entry = 1
    L%RF_MOSE_recv%n_entries = nrm_r
    L%RF_MOSE_recv%payload_per_entry = 1
    call copy_groups(L%M_recv, L%RF_MOSE_send)
    call copy_groups(L%M_send, L%RF_MOSE_recv)
    call alloc_dir_buf(L%RF_MOSE_send)
    call alloc_dir_buf(L%RF_MOSE_recv)

    allocate(L%RF_MOSE_recv_fussbc(max(nrm_r,1)))
    L%RF_MOSE_recv_fussbc = 0
    do k = 1, nrm_r
      i = L%M_send_ent(k)%global_bc_idx
      do j = 1, FUSS_domain%nbound
        if (FUSS_domain%bc(j)%type /= 103) cycle
        if (MOSE_domain%bc(i)%bs == FUSS_domain%bc(j)%b  .and. &
            MOSE_domain%bc(i)%is == FUSS_domain%bc(j)%i  .and. &
            MOSE_domain%bc(i)%js == FUSS_domain%bc(j)%j  .and. &
            MOSE_domain%bc(i)%ks == FUSS_domain%bc(j)%k  .and. &
            MOSE_domain%bc(i)%fs == FUSS_domain%bc(j)%f) then
          L%RF_MOSE_recv_fussbc(k) = j
          exit
        end if
      end do
    end do
    n = 0  ! silence unused-var warnings in some compilers
  end subroutine build_reverse_flux_schedule
#endif


  !---------------------------------------------------------------------------
  !> Exchange Tg(6) + P(1:nsc,Is,Js,Ks) for every FUSS-primary BC with remote Bs.
  !> MOSE-owner of Bs packs and sends; FUSS-owner of Bm receives and writes into
  !> bc(i)%Tg and L%F_remote_P(:, bc_idx).  Must be called BEFORE the local
  !> Fill_Coupled_Ghost_Cell loop, which is then guarded to skip remote entries.
  subroutine exchange_coupling_Tg_P(FUSS_domain, MOSE_domain, level)
    use MOSE_Global_m, only: nsc
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    integer, intent(in) :: level
#ifdef USE_MPI
    type(level_sched_t), pointer :: L
    integer :: g, k, off, i
    integer :: payload

    if (mpi_size_ <= 1) return
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    L => levels(level)
    if (.not. L%built) return
    payload = L%F_send%payload_per_entry
    if (L%F_recv%n_entries == 0 .and. L%F_send%n_entries == 0) return

    ! Post all receives first
    call post_recvs(L%F_recv, TAG_TG_P)

    ! Pack send buffer: loop over local F_send entries, compute Tg from MOSE block,
    ! append P(1:nsc, Is,Js,Ks).
    do k = 1, L%F_send%n_entries
      associate (e => L%F_send_ent(k))
        call pack_Tg_P_entry(MOSE_domain%blk(e%local_bs), &
                             e%is, e%js, e%ks, e%fs, &
                             e%d11, e%d12, e%d21, e%d22, &
                             L%F_send%buf((k-1)*payload + 1 : k*payload))
      end associate
    end do

    call post_sends(L%F_send, TAG_TG_P)
    call wait_all(L%F_recv)

    ! Unpack recv buffer into FUSS bc(i)%Tg and L%F_remote_P
    payload = L%F_recv%payload_per_entry
    do k = 1, L%F_recv%n_entries
      i = L%F_recv_ent(k)%bc_idx
      off = (k-1)*payload
      do g = 1, NSTENC
        FUSS_domain%bc(i)%Tg(g) = L%F_recv%buf(off + g)
      end do
      L%F_remote_P(1:nsc, i) = L%F_recv%buf(off + NSTENC + 1 : off + NSTENC + nsc)
    end do

    call wait_all(L%F_send)
#endif
  end subroutine exchange_coupling_Tg_P


#ifdef USE_MPI
  !---------------------------------------------------------------------------
  !> Pack one F_send entry: [Tg(1..6), P(1:nsc, Is,Js,Ks)].
  !> Uses MOSE primitive field to compute Tg via the same formulas as
  !> Ghost_Coupling_FUSS + Fill_BC_Ghost_Coupling_FUSS (local-only access).
  subroutine pack_Tg_P_entry(blks, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, buf_out)
    use FLINT_Lib_Thermodynamic
    use MOSE_Global_m, only: nsc, np
    use hydra_Parameters_m, only: guide_
    use hydra_Global_m, only: gc_
    use MOSE_Advanced_Types_m, only: MOSE_block_type
    implicit none
    type(MOSE_block_type), intent(in) :: blks
    integer, intent(in)  :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    real(R8), intent(out) :: buf_out(:)
    ! Local
    integer :: g, Is_, Js_, Ks_
    real(R8) :: rho, Rgas
    integer :: i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4

    ! Normal-direction ghost cells (indices 1..gc_ of Tg)
    do g = 1, gc_
      Is_ = Is + guide_(Fs,1)*(g-1)
      Js_ = Js + guide_(Fs,2)*(g-1)
      Ks_ = Ks + guide_(Fs,3)*(g-1)
      call co_rotot_Rtot ( blks%P(1:nsc,Is_,Js_,Ks_), rho, Rgas )
      buf_out(g) = EOS( p=blks%P(np,Is_,Js_,Ks_), rho=rho, R=Rgas )
    end do

    ! Tangential corner cells (indices gc_+1..gc_+4)
    select case (Fs)
      case (1:2)
        i1=Is;       j1=Js-d11s;  k1=Ks-d12s
        i2=Is;       j2=Js+d11s;  k2=Ks+d12s
        i3=Is;       j3=Js-d21s;  k3=Ks-d22s
        i4=Is;       j4=Js+d21s;  k4=Ks+d22s
      case (3:4)
        i1=Is-d11s;  j1=Js;       k1=Ks-d12s
        i2=Is+d11s;  j2=Js;       k2=Ks+d12s
        i3=Is-d21s;  j3=Js;       k3=Ks-d22s
        i4=Is+d21s;  j4=Js;       k4=Ks+d22s
      case (5:6)
        i1=Is-d11s;  j1=Js-d12s;  k1=Ks
        i2=Is+d11s;  j2=Js+d12s;  k2=Ks
        i3=Is-d21s;  j3=Js-d22s;  k3=Ks
        i4=Is+d21s;  j4=Js+d22s;  k4=Ks
    end select

    call co_rotot_Rtot ( blks%P(1:nsc,i1,j1,k1), rho, Rgas )
    buf_out(gc_+1) = EOS( p=blks%P(np,i1,j1,k1), rho=rho, R=Rgas )
    call co_rotot_Rtot ( blks%P(1:nsc,i2,j2,k2), rho, Rgas )
    buf_out(gc_+2) = EOS( p=blks%P(np,i2,j2,k2), rho=rho, R=Rgas )
    call co_rotot_Rtot ( blks%P(1:nsc,i3,j3,k3), rho, Rgas )
    buf_out(gc_+3) = EOS( p=blks%P(np,i3,j3,k3), rho=rho, R=Rgas )
    call co_rotot_Rtot ( blks%P(1:nsc,i4,j4,k4), rho, Rgas )
    buf_out(gc_+4) = EOS( p=blks%P(np,i4,j4,k4), rho=rho, R=Rgas )

    ! Trailing nsc reals: primitive composition at the BC interior cell (Is,Js,Ks)
    buf_out(NSTENC+1 : NSTENC+nsc) = blks%P(1:nsc, Is, Js, Ks)
  end subroutine pack_Tg_P_entry
#endif


  !---------------------------------------------------------------------------
  !> Exchange Pg(1,1..6) + matID(Is,Js,Ks) for every MOSE-primary BC with remote Bs.
  subroutine exchange_coupling_Pg_matID(FUSS_domain, MOSE_domain, level)
    implicit none
    type(FUSS_domain_type), intent(in)    :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level
#ifdef USE_MPI
    type(level_sched_t), pointer :: L
    integer :: g, k, i, off, payload

    if (mpi_size_ <= 1) return
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    L => levels(level)
    if (.not. L%built) return
    if (L%M_recv%n_entries == 0 .and. L%M_send%n_entries == 0) return

    call post_recvs(L%M_recv, TAG_PG_MID)

    payload = L%M_send%payload_per_entry
    do k = 1, L%M_send%n_entries
      associate (e => L%M_send_ent(k))
        call pack_Pg_matID_entry(FUSS_domain%blk(e%local_bs), &
                                 e%is, e%js, e%ks, e%fs, &
                                 e%d11, e%d12, e%d21, e%d22, &
                                 L%M_send%buf((k-1)*payload + 1 : k*payload))
      end associate
    end do

    call post_sends(L%M_send, TAG_PG_MID)
    call wait_all(L%M_recv)

    payload = L%M_recv%payload_per_entry
    do k = 1, L%M_recv%n_entries
      i = L%M_recv_ent(k)%bc_idx
      off = (k-1)*payload
      do g = 1, NSTENC
        MOSE_domain%bc(i)%Pg(1,g) = L%M_recv%buf(off + g)
      end do
      L%M_remote_matID(i) = L%M_recv%buf(off + NSTENC + 1)
    end do

    call wait_all(L%M_send)
#endif
  end subroutine exchange_coupling_Pg_matID


#ifdef USE_MPI
  !---------------------------------------------------------------------------
  subroutine pack_Pg_matID_entry(blks, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, buf_out)
    use hydra_Parameters_m, only: guide_
    use hydra_Global_m, only: gc_
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(in) :: blks
    integer, intent(in) :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    real(R8), intent(out) :: buf_out(:)
    ! Local
    integer :: g, Is_, Js_, Ks_
    integer :: i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4

    do g = 1, gc_
      Is_ = Is + guide_(Fs,1)*(g-1)
      Js_ = Js + guide_(Fs,2)*(g-1)
      Ks_ = Ks + guide_(Fs,3)*(g-1)
      buf_out(g) = blks%T(Is_,Js_,Ks_)
    end do

    select case (Fs)
      case (1:2)
        i1=Is;       j1=Js-d11s;  k1=Ks-d12s
        i2=Is;       j2=Js+d11s;  k2=Ks+d12s
        i3=Is;       j3=Js-d21s;  k3=Ks-d22s
        i4=Is;       j4=Js+d21s;  k4=Ks+d22s
      case (3:4)
        i1=Is-d11s;  j1=Js;       k1=Ks-d12s
        i2=Is+d11s;  j2=Js;       k2=Ks+d12s
        i3=Is-d21s;  j3=Js;       k3=Ks-d22s
        i4=Is+d21s;  j4=Js;       k4=Ks+d22s
      case (5:6)
        i1=Is-d11s;  j1=Js-d12s;  k1=Ks
        i2=Is+d11s;  j2=Js+d12s;  k2=Ks
        i3=Is-d21s;  j3=Js-d22s;  k3=Ks
        i4=Is+d21s;  j4=Js+d22s;  k4=Ks
    end select

    buf_out(gc_+1) = blks%T(i1,j1,k1)
    buf_out(gc_+2) = blks%T(i2,j2,k2)
    buf_out(gc_+3) = blks%T(i3,j3,k3)
    buf_out(gc_+4) = blks%T(i4,j4,k4)

    buf_out(NSTENC+1) = blks%matID(Is,Js,Ks)
  end subroutine pack_Pg_matID_entry
#endif


  !---------------------------------------------------------------------------
  !> One-off exchange of the secondary-block node coordinates needed to build
  !> the coupling metric tensor Mg of BCs whose two blocks live on different
  !> ranks.
  !>
  !> Setup_Coupled_Metrics can only build Mg for entries whose secondary block
  !> is local, because it reads that block's node coordinates directly. Without
  !> this exchange the remaining entries keep Mg = 0, and the wall-temperature
  !> routines then reduce to Tint = Num_loc/Den_loc: the secondary side drops
  !> out of the interface temperature altogether, silently and with no error.
  !>
  !> Call once per MG level from Setup_Coupled_Metrics and release the buffers
  !> with free_coupling_metric_nodes as soon as the metrics have been built.
  subroutine exchange_coupling_metric_nodes(FUSS_domain, MOSE_domain, level)
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in) :: level
#ifdef USE_MPI
    type(level_sched_t), pointer :: L
    integer :: k, payload

    if (mpi_size_ <= 1) return
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    L => levels(level)
    if (.not. L%built) return

    payload = 3 * NMETA_NODE * gc_

    L%F_meta_recv%n_entries = L%F_recv%n_entries
    L%F_meta_send%n_entries = L%F_send%n_entries
    L%M_meta_recv%n_entries = L%M_recv%n_entries
    L%M_meta_send%n_entries = L%M_send%n_entries
    L%F_meta_recv%payload_per_entry = payload
    L%F_meta_send%payload_per_entry = payload
    L%M_meta_recv%payload_per_entry = payload
    L%M_meta_send%payload_per_entry = payload
    call copy_groups(L%F_recv, L%F_meta_recv)
    call copy_groups(L%F_send, L%F_meta_send)
    call copy_groups(L%M_recv, L%M_meta_recv)
    call copy_groups(L%M_send, L%M_meta_send)
    call alloc_dir_buf(L%F_meta_recv)
    call alloc_dir_buf(L%F_meta_send)
    call alloc_dir_buf(L%M_meta_recv)
    call alloc_dir_buf(L%M_meta_send)

    ! FUSS-primary entries: the MOSE owner of Bs packs its node coordinates.
    if (L%F_meta_recv%n_entries > 0 .or. L%F_meta_send%n_entries > 0) then
      call post_recvs(L%F_meta_recv, TAG_META_F)
      do k = 1, L%F_meta_send%n_entries
        associate (e => L%F_send_ent(k))
          call pack_metric_nodes_MOSE(MOSE_domain%blk(e%local_bs), &
                                      e%is, e%js, e%ks, e%fs, &
                                      e%d11, e%d12, e%d21, e%d22, &
                                      L%F_meta_send%buf((k-1)*payload + 1 : k*payload))
        end associate
      end do
      call post_sends(L%F_meta_send, TAG_META_F)
      call wait_all(L%F_meta_recv)
      call wait_all(L%F_meta_send)
    end if

    ! MOSE-primary entries: the FUSS owner of Bs packs its node coordinates.
    if (L%M_meta_recv%n_entries > 0 .or. L%M_meta_send%n_entries > 0) then
      call post_recvs(L%M_meta_recv, TAG_META_M)
      do k = 1, L%M_meta_send%n_entries
        associate (e => L%M_send_ent(k))
          call pack_metric_nodes_FUSS(FUSS_domain%blk(e%local_bs), &
                                      e%is, e%js, e%ks, e%fs, &
                                      e%d11, e%d12, e%d21, e%d22, &
                                      L%M_meta_send%buf((k-1)*payload + 1 : k*payload))
        end associate
      end do
      call post_sends(L%M_meta_send, TAG_META_M)
      call wait_all(L%M_meta_recv)
      call wait_all(L%M_meta_send)
    end if
#endif
    ! Dummy usage silences unused-arg warnings in serial builds
    if (.false.) then
      if (FUSS_domain%nbound < 0) return
      if (MOSE_domain%nbound < 0) return
      if (level < 0) return
    end if
  end subroutine exchange_coupling_metric_nodes


  !---------------------------------------------------------------------------
  !> Release the metric-node buffers once Mg has been built.
  subroutine free_coupling_metric_nodes(level)
    implicit none
    integer, intent(in) :: level
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    call reset_dir(levels(level)%F_meta_recv)
    call reset_dir(levels(level)%F_meta_send)
    call reset_dir(levels(level)%M_meta_recv)
    call reset_dir(levels(level)%M_meta_send)
  end subroutine free_coupling_metric_nodes


  !---------------------------------------------------------------------------
  !> Number of FUSS-primary / MOSE-primary BCs on this rank whose secondary node
  !> coordinates arrived over MPI. Zero in serial, where every BC is local-local.
  integer function n_F_meta_entries(level) result(n)
    implicit none
    integer, intent(in) :: level
    n = 0
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    if (.not. levels(level)%built) return
    n = levels(level)%F_meta_recv%n_entries
  end function n_F_meta_entries

  integer function n_M_meta_entries(level) result(n)
    implicit none
    integer, intent(in) :: level
    n = 0
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    if (.not. levels(level)%built) return
    n = levels(level)%M_meta_recv%n_entries
  end function n_M_meta_entries


  !---------------------------------------------------------------------------
  !> Entry k of those lists: the local BC index it belongs to, and the received
  !> node coordinates as snode(3, NMETA_NODE, gc_), laid out exactly as
  !> BC_Coupled_Metrics_* expects them.
  subroutine get_F_meta_entry(level, k, bc_idx, snode)
    implicit none
    integer,  intent(in)  :: level, k
    integer,  intent(out) :: bc_idx
    real(R8), intent(out) :: snode(:,:,:)
    integer :: g, n, off, payload

    bc_idx = 0
    snode  = 0.0d0
    if (k < 1 .or. k > n_F_meta_entries(level)) return
#ifdef USE_MPI
    bc_idx  = levels(level)%F_recv_ent(k)%bc_idx
    payload = levels(level)%F_meta_recv%payload_per_entry
    do g = 1, gc_
      do n = 1, NMETA_NODE
        off = (k-1)*payload + ((g-1)*NMETA_NODE + (n-1))*3
        snode(1:3, n, g) = levels(level)%F_meta_recv%buf(off+1:off+3)
      end do
    end do
#endif
  end subroutine get_F_meta_entry

  subroutine get_M_meta_entry(level, k, bc_idx, snode)
    implicit none
    integer,  intent(in)  :: level, k
    integer,  intent(out) :: bc_idx
    real(R8), intent(out) :: snode(:,:,:)
    integer :: g, n, off, payload

    bc_idx = 0
    snode  = 0.0d0
    if (k < 1 .or. k > n_M_meta_entries(level)) return
#ifdef USE_MPI
    bc_idx  = levels(level)%M_recv_ent(k)%bc_idx
    payload = levels(level)%M_meta_recv%payload_per_entry
    do g = 1, gc_
      do n = 1, NMETA_NODE
        off = (k-1)*payload + ((g-1)*NMETA_NODE + (n-1))*3
        snode(1:3, n, g) = levels(level)%M_meta_recv%buf(off+1:off+3)
      end do
    end do
#endif
  end subroutine get_M_meta_entry


#ifdef USE_MPI
  !---------------------------------------------------------------------------
  !> Pack the gc_ * NMETA_NODE secondary node coordinates of one entry.
  subroutine pack_metric_nodes_MOSE(blks, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, buf_out)
    use MOSE_Advanced_Types_m, only: MOSE_block_type
    implicit none
    type(MOSE_block_type), intent(in) :: blks
    integer,  intent(in)  :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    real(R8), intent(out) :: buf_out(:)
    ! Local
    integer :: g, n, off, idx(3,NMETA_NODE)

    do g = 1, gc_
      call coupled_secondary_node_indices(g, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, idx)
      do n = 1, NMETA_NODE
        off = ((g-1)*NMETA_NODE + (n-1)) * 3
        buf_out(off+1:off+3) = blks%node(idx(1,n), idx(2,n), idx(3,n))%c
      end do
    end do
  end subroutine pack_metric_nodes_MOSE

  subroutine pack_metric_nodes_FUSS(blks, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, buf_out)
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(in) :: blks
    integer,  intent(in)  :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    real(R8), intent(out) :: buf_out(:)
    ! Local
    integer :: g, n, off, idx(3,NMETA_NODE)

    do g = 1, gc_
      call coupled_secondary_node_indices(g, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, idx)
      do n = 1, NMETA_NODE
        off = ((g-1)*NMETA_NODE + (n-1)) * 3
        buf_out(off+1:off+3) = blks%node(idx(1,n), idx(2,n), idx(3,n))%c
      end do
    end do
  end subroutine pack_metric_nodes_FUSS
#endif


  !---------------------------------------------------------------------------
  !> After Compute_FUSS_Flux has written FUSS%bc(i)%ext_flux = 0.5*FUSS_new_flux
  !> for each FUSS type(103) entry, this routine packs the FUSS_new_flux value
  !> for each ENTRY WITH REMOTE Bs and sends it to the MOSE-owner, which
  !> subtracts 0.5*FUSS_new_flux from partner MOSE%bc(j)%ext_flux(np).
  !> Expects the caller to provide the per-entry FUSS_new_flux in a temp array
  !> matching FUSS%bc index space.
  subroutine exchange_coupling_reverse_flux_FUSS(FUSS_new_flux_arr, FUSS_domain, MOSE_domain, level)
    use MOSE_Global_m, only: nprim, np
    implicit none
    real(R8), intent(in) :: FUSS_new_flux_arr(:)   !< (FUSS%nbound), 0 where unused
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level
#ifdef USE_MPI
    type(level_sched_t), pointer :: L
    integer :: k, i, j

    if (mpi_size_ <= 1) return
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    L => levels(level)
    if (.not. L%built) return
    if (L%RF_FUSS_send%n_entries == 0 .and. L%RF_FUSS_recv%n_entries == 0) return

    call post_recvs(L%RF_FUSS_recv, TAG_RFLUX_F)

    ! Pack: L%RF_FUSS_send is grouped identically to L%F_recv; entry k maps to
    ! L%F_recv_ent(k)%bc_idx which is the local FUSS bc index.
    do k = 1, L%RF_FUSS_send%n_entries
      i = L%F_recv_ent(k)%bc_idx
      L%RF_FUSS_send%buf(k) = FUSS_new_flux_arr(i)
    end do

    call post_sends(L%RF_FUSS_send, TAG_RFLUX_F)
    call wait_all(L%RF_FUSS_recv)

    ! Apply: for each recv slot, locate partner MOSE bc and subtract 0.5*flux from ext_flux(np).
    do k = 1, L%RF_FUSS_recv%n_entries
      j = L%RF_FUSS_recv_mosebc(k)
      if (j > 0) then
        MOSE_domain%bc(j)%ext_flux(np) = MOSE_domain%bc(j)%ext_flux(np) &
                                       - 0.5d0 * L%RF_FUSS_recv%buf(k)
      end if
    end do

    call wait_all(L%RF_FUSS_send)
#endif
    ! Dummy usage silences unused-arg warnings under some compilers
    if (.false.) then
      if (size(FUSS_new_flux_arr) < 0) return
      if (FUSS_domain%nbound < 0) return
      if (MOSE_domain%nbound < 0) return
    end if
  end subroutine exchange_coupling_reverse_flux_FUSS


  !---------------------------------------------------------------------------
  subroutine exchange_coupling_reverse_flux_MOSE(MOSE_new_flux_np_arr, FUSS_domain, MOSE_domain, level)
    implicit none
    real(R8), intent(in) :: MOSE_new_flux_np_arr(:)  !< (MOSE%nbound), 0 where unused
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in) :: level
#ifdef USE_MPI
    type(level_sched_t), pointer :: L
    integer :: k, i, j

    if (mpi_size_ <= 1) return
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    L => levels(level)
    if (.not. L%built) return
    if (L%RF_MOSE_send%n_entries == 0 .and. L%RF_MOSE_recv%n_entries == 0) return

    call post_recvs(L%RF_MOSE_recv, TAG_RFLUX_M)

    do k = 1, L%RF_MOSE_send%n_entries
      i = L%M_recv_ent(k)%bc_idx
      L%RF_MOSE_send%buf(k) = MOSE_new_flux_np_arr(i)
    end do

    call post_sends(L%RF_MOSE_send, TAG_RFLUX_M)
    call wait_all(L%RF_MOSE_recv)

    do k = 1, L%RF_MOSE_recv%n_entries
      j = L%RF_MOSE_recv_fussbc(k)
      if (j > 0) then
        FUSS_domain%bc(j)%ext_flux = FUSS_domain%bc(j)%ext_flux &
                                   - 0.5d0 * L%RF_MOSE_recv%buf(k)
      end if
    end do

    call wait_all(L%RF_MOSE_send)
#endif
    if (.false.) then
      if (size(MOSE_new_flux_np_arr) < 0) return
      if (FUSS_domain%nbound < 0) return
      if (MOSE_domain%nbound < 0) return
    end if
  end subroutine exchange_coupling_reverse_flux_MOSE


  !---------------------------------------------------------------------------
  !> Accessors for refactored consumers
  pure function partner_of_FUSS_bc(level, i) result(j)
    integer, intent(in) :: level, i
    integer :: j
    j = 0
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    if (.not. allocated(levels(level)%partner_of_FUSS)) return
    if (i < 1 .or. i > size(levels(level)%partner_of_FUSS)) return
    j = levels(level)%partner_of_FUSS(i)
  end function partner_of_FUSS_bc

  pure function partner_of_MOSE_bc(level, i) result(j)
    integer, intent(in) :: level, i
    integer :: j
    j = 0
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    if (.not. allocated(levels(level)%partner_of_MOSE)) return
    if (i < 1 .or. i > size(levels(level)%partner_of_MOSE)) return
    j = levels(level)%partner_of_MOSE(i)
  end function partner_of_MOSE_bc

  !> True iff both Bm and Bs for this FUSS BC are owned locally (partner is local).
  function is_FUSS_bc_local_local(FUSS_domain, i) result(lr)
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    integer, intent(in) :: i
    logical :: lr
    integer :: Bm, Bs
#ifdef USE_MPI
    if (.not. allocated(FUSS_block_owner) .or. .not. allocated(MOSE_block_owner)) then
      lr = .true.; return
    end if
    Bm = FUSS_domain%bc(i)%b
    Bs = FUSS_domain%bc(i)%bs
    lr = (FUSS_block_owner(Bm) == mpi_rank_) .and. (MOSE_block_owner(Bs) == mpi_rank_)
#else
    lr = .true.
    if (FUSS_domain%bc(i)%b < 0) Bm = 0  ! silence unused
    if (FUSS_domain%bc(i)%bs < 0) Bs = 0
#endif
  end function is_FUSS_bc_local_local

  function is_MOSE_bc_local_local(MOSE_domain, i) result(lr)
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in) :: i
    logical :: lr
    integer :: Bm, Bs
#ifdef USE_MPI
    if (.not. allocated(MOSE_block_owner) .or. .not. allocated(FUSS_block_owner)) then
      lr = .true.; return
    end if
    Bm = MOSE_domain%bc(i)%b
    Bs = MOSE_domain%bc(i)%bs
    lr = (MOSE_block_owner(Bm) == mpi_rank_) .and. (FUSS_block_owner(Bs) == mpi_rank_)
#else
    lr = .true.
    if (MOSE_domain%bc(i)%b < 0) Bm = 0
    if (MOSE_domain%bc(i)%bs < 0) Bs = 0
#endif
  end function is_MOSE_bc_local_local


  !===========================================================================
  ! Internal helpers
  !===========================================================================

  subroutine reset_level(L)
    implicit none
    type(level_sched_t), intent(inout) :: L
    if (allocated(L%F_recv_ent)) deallocate(L%F_recv_ent)
    if (allocated(L%F_send_ent)) deallocate(L%F_send_ent)
    if (allocated(L%M_recv_ent)) deallocate(L%M_recv_ent)
    if (allocated(L%M_send_ent)) deallocate(L%M_send_ent)
    if (allocated(L%F_remote_P)) deallocate(L%F_remote_P)
    if (allocated(L%M_remote_matID)) deallocate(L%M_remote_matID)
    if (allocated(L%partner_of_FUSS)) deallocate(L%partner_of_FUSS)
    if (allocated(L%partner_of_MOSE)) deallocate(L%partner_of_MOSE)
    if (allocated(L%own103_F)) deallocate(L%own103_F)
    if (allocated(L%own103_M)) deallocate(L%own103_M)
    if (allocated(L%rflux_F)) deallocate(L%rflux_F)
    if (allocated(L%rflux_M)) deallocate(L%rflux_M)
    L%n_own103_F = 0
    L%n_own103_M = 0
    if (allocated(L%RF_FUSS_recv_mosebc)) deallocate(L%RF_FUSS_recv_mosebc)
    if (allocated(L%RF_MOSE_recv_fussbc)) deallocate(L%RF_MOSE_recv_fussbc)
    call reset_dir(L%F_recv)
    call reset_dir(L%F_send)
    call reset_dir(L%M_recv)
    call reset_dir(L%M_send)
    call reset_dir(L%RF_FUSS_send)
    call reset_dir(L%RF_FUSS_recv)
    call reset_dir(L%RF_MOSE_send)
    call reset_dir(L%RF_MOSE_recv)
    call reset_dir(L%F_meta_recv)
    call reset_dir(L%F_meta_send)
    call reset_dir(L%M_meta_recv)
    call reset_dir(L%M_meta_send)
    L%built = .false.
  end subroutine reset_level

  subroutine reset_dir(d)
    implicit none
    type(dir_t), intent(inout) :: d
    if (allocated(d%groups)) deallocate(d%groups)
    if (allocated(d%buf)) deallocate(d%buf)
    if (allocated(d%reqs)) deallocate(d%reqs)
    d%n_entries = 0
    d%n_groups  = 0
    d%payload_per_entry = 0
  end subroutine reset_dir


#ifdef USE_MPI
  subroutine sort_recv_by_rank(arr, n)
    implicit none
    type(recv_entry_t), intent(inout) :: arr(:)
    integer, intent(in) :: n
    integer :: i, j
    type(recv_entry_t) :: key
    do i = 2, n
      key = arr(i)
      j = i - 1
      do while (j >= 1)
        if (arr(j)%remote_rank <= key%remote_rank) exit
        arr(j+1) = arr(j)
        j = j - 1
      end do
      arr(j+1) = key
    end do
  end subroutine sort_recv_by_rank

  subroutine sort_send_by_rank(arr, n)
    implicit none
    type(send_entry_t), intent(inout) :: arr(:)
    integer, intent(in) :: n
    integer :: i, j
    type(send_entry_t) :: key
    do i = 2, n
      key = arr(i)
      j = i - 1
      do while (j >= 1)
        if (arr(j)%remote_rank <= key%remote_rank) exit
        arr(j+1) = arr(j)
        j = j - 1
      end do
      arr(j+1) = key
    end do
  end subroutine sort_send_by_rank

  subroutine build_groups_recv(arr, n, d)
    implicit none
    type(recv_entry_t), intent(in) :: arr(:)
    integer, intent(in) :: n
    type(dir_t), intent(inout) :: d
    integer :: i, g, cur
    if (n <= 0) then
      d%n_groups = 0
      allocate(d%groups(0))
      return
    end if
    ! count
    g = 1
    do i = 2, n
      if (arr(i)%remote_rank /= arr(i-1)%remote_rank) g = g + 1
    end do
    d%n_groups = g
    allocate(d%groups(g))
    ! fill
    g = 1; cur = arr(1)%remote_rank
    d%groups(1)%rank = cur
    d%groups(1)%offset = 1
    d%groups(1)%count = 1
    do i = 2, n
      if (arr(i)%remote_rank /= cur) then
        g = g + 1
        cur = arr(i)%remote_rank
        d%groups(g)%rank = cur
        d%groups(g)%offset = i
        d%groups(g)%count = 1
      else
        d%groups(g)%count = d%groups(g)%count + 1
      end if
    end do
  end subroutine build_groups_recv

  subroutine build_groups_send(arr, n, d)
    implicit none
    type(send_entry_t), intent(in) :: arr(:)
    integer, intent(in) :: n
    type(dir_t), intent(inout) :: d
    integer :: i, g, cur
    if (n <= 0) then
      d%n_groups = 0
      allocate(d%groups(0))
      return
    end if
    g = 1
    do i = 2, n
      if (arr(i)%remote_rank /= arr(i-1)%remote_rank) g = g + 1
    end do
    d%n_groups = g
    allocate(d%groups(g))
    g = 1; cur = arr(1)%remote_rank
    d%groups(1)%rank = cur
    d%groups(1)%offset = 1
    d%groups(1)%count = 1
    do i = 2, n
      if (arr(i)%remote_rank /= cur) then
        g = g + 1
        cur = arr(i)%remote_rank
        d%groups(g)%rank = cur
        d%groups(g)%offset = i
        d%groups(g)%count = 1
      else
        d%groups(g)%count = d%groups(g)%count + 1
      end if
    end do
  end subroutine build_groups_send

  subroutine copy_groups(src, dst)
    implicit none
    type(dir_t), intent(in) :: src
    type(dir_t), intent(inout) :: dst
    integer :: g
    dst%n_groups = src%n_groups
    if (allocated(dst%groups)) deallocate(dst%groups)
    allocate(dst%groups(max(src%n_groups,0)))
    do g = 1, src%n_groups
      dst%groups(g) = src%groups(g)
    end do
  end subroutine copy_groups

  subroutine alloc_dir_buf(d)
    implicit none
    type(dir_t), intent(inout) :: d
    integer :: total
    total = max(d%n_entries * d%payload_per_entry, 1)
    if (allocated(d%buf)) deallocate(d%buf)
    allocate(d%buf(total)); d%buf = 0.0d0
    if (allocated(d%reqs)) deallocate(d%reqs)
    allocate(d%reqs(max(d%n_groups,1))); d%reqs = MPI_REQUEST_NULL
  end subroutine alloc_dir_buf

  subroutine post_recvs(d, tag)
    use hydra_Mod_Timers, only: timer_xchg_begin, timer_xchg_end
    implicit none
    type(dir_t), intent(inout) :: d
    integer, intent(in) :: tag
    integer :: g, ierr, off, cnt
    call timer_xchg_begin()
    do g = 1, d%n_groups
      off = (d%groups(g)%offset - 1) * d%payload_per_entry
      cnt = d%groups(g)%count * d%payload_per_entry
      call MPI_IRECV(d%buf(off+1), cnt, MPI_DOUBLE_PRECISION, &
                     d%groups(g)%rank, tag, MPI_COMM_WORLD, d%reqs(g), ierr)
      call check_mpi_error(ierr)
    end do
    call timer_xchg_end()
  end subroutine post_recvs

  subroutine post_sends(d, tag)
    use hydra_Mod_Timers, only: timer_xchg_begin, timer_xchg_end
    implicit none
    type(dir_t), intent(inout) :: d
    integer, intent(in) :: tag
    integer :: g, ierr, off, cnt
    call timer_xchg_begin()
    do g = 1, d%n_groups
      off = (d%groups(g)%offset - 1) * d%payload_per_entry
      cnt = d%groups(g)%count * d%payload_per_entry
      call MPI_ISEND(d%buf(off+1), cnt, MPI_DOUBLE_PRECISION, &
                     d%groups(g)%rank, tag, MPI_COMM_WORLD, d%reqs(g), ierr)
      call check_mpi_error(ierr)
    end do
    call timer_xchg_end()
  end subroutine post_sends

  subroutine wait_all(d)
    use hydra_Mod_Timers, only: timer_xchg_begin, timer_xchg_end
    implicit none
    type(dir_t), intent(inout) :: d
    integer :: ierr
    integer, allocatable :: statuses(:,:)
    if (d%n_groups <= 0) return
    call timer_xchg_begin()
    allocate(statuses(MPI_STATUS_SIZE, d%n_groups))
    call MPI_WAITALL(d%n_groups, d%reqs, statuses, ierr)
    call check_mpi_error(ierr)
    deallocate(statuses)
    call timer_xchg_end()
  end subroutine wait_all
#endif


  !---------------------------------------------------------------------------
  !> Accessors for the per-BC remote storage (used by Wall_Temperature routines
  !> when the BC's secondary block is owned by a different rank).
  function get_F_remote_P(level, i) result(ptr)
    integer, intent(in) :: level, i
    real(R8), pointer :: ptr(:)
    ptr => null()
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    if (.not. allocated(levels(level)%F_remote_P)) return
    if (i < 1 .or. i > size(levels(level)%F_remote_P, 2)) return
    ptr => levels(level)%F_remote_P(:, i)
  end function get_F_remote_P

  function get_M_remote_matID(level, i) result(v)
    integer, intent(in) :: level, i
    real(R8) :: v
    v = 0.0d0
    if (.not. allocated(levels)) return
    if (level < 1 .or. level > size(levels)) return
    if (.not. allocated(levels(level)%M_remote_matID)) return
    if (i < 1 .or. i > size(levels(level)%M_remote_matID)) return
    v = levels(level)%M_remote_matID(i)
  end function get_M_remote_matID

end module hydraMF_Mod_CouplingExchange
