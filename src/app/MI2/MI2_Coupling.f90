!> The exchange between MOSE and the two condensed solvers of hydra-MI2.
!
!  Within a coupling step:
!    MOSExIGLOO, MOSExICE : the gas state a condensed solver steps on
!    IGLOOxMOSE, ICExMOSE : the sources it hands back, into MOSE's residual at every Runge-Kutta stage
!  What a solver hands the gas in one cell is written once, in the cell kernels IGLOO_Mass_To_Gas, ICE_Cell and
!  Removal_Cell: the injection applies them, IGLOO_Handed, ICE_Handed and the removals' totals sum them for the
!  books and the report (hydraMI2_Report).
!
!  What lands outside a solver's step goes through a channel (hydraMI2_Channels), one per route:
!    removed(g)  ICE family g -> gas: what its burn-out removed (cond_phase%gone, drained after each ICE step), which
!                the gas has not received yet
!    iceMean(g)  ICE family g -> gas: a time-accurate ICE block of several steps hands the gas its steps' dt-weighted
!                mean source
!    smoke(c)    IGLOO material -> ICE family: the condensed product of IGLOO's burning, born in the family instead
!                of the gas
!  They go as rates when MOSE or ICE steps in pseudo-time (handRates: a converged steady field then holds per unit
!  time), as amounts when both are time-accurate.
!
!  When MOSE starts on its coarse grids, IGLOO and ICE still see its finest grid, filled from the level MOSE steps
!  (Fill_Finest_Grid), and that level takes their sources restricted to it from the finest grid (Before_MOSE_Iteration).
module hydraMI2_Mod_Coupling
  use iso_fortran_env,      only: R8 => real64
  use hydraMI2_Channels,    only: channel_t
  use IGLOO_data_block,     only: obj_sourceblock, obj_eulerblock
  use ICE_Advanced_Types_m, only: ICE_domain_type

  implicit none
  private

  public :: source_state_t, handed_t, removals_t, smoke_t, restricted_t
  public :: Setup_ICE_Coupling, Setup_IGLOO_Coupling, Set_IGLOO_Switches
  public :: Gas_Field_For_IGLOO, MOSExIGLOO, MOSExICE, IGLOOxMOSE, ICExMOSE
  public :: Check_IGLOO_Sources, Check_ICE_Sources
  public :: IGLOO_Handed, ICE_Handed, Removals_Handed
  public :: Removed_Store_Setup, Removals_Timed
  public :: Open_Removals, Close_Removals, After_ICE_Step, ICE_Block_Mean
  public :: IGLOO_Smoke_Accumulate, IGLOO_Smoke_Deliver
  public :: Before_MOSE_Iteration, Fill_Finest_Grid, Evaluate_Gas_Once, Restricted_Totals

  !> The sources a condensed solver last produced. They stay in force, and MOSE keeps reading them, until that
  !  solver's next solve: a cycle with X-iter = 0 runs on frozen sources. cycle = 0 for sources a restart read back.
  type :: source_state_t
    logical :: ready    = .false.   ! the solver has produced sources MOSE may read
    logical :: eulerian = .false.   ! IGLOO only: they were produced under the Eulerian coupling mode
    integer :: cycle    = 0         ! the (cycle, coupling step) that produced them
    integer :: step     = 0
  end type source_state_t

  !> What a condensed solver hands the gas, summed over the domain (a rate, or an amount): the mass of each gas
  !  species (allocated when a material converts), momentum, energy, and the part of the energy that is datum offsets.
  type :: handed_t
    real(R8), allocatable :: mass(:)
    real(R8) :: mom(3) = 0._R8, energy = 0._R8, datum = 0._R8
  end type handed_t

  !> What ICE's stored removals hand the gas over a MOSE block (Open_Removals ... Close_Removals): the rates in force
  !  over it (asRate), or the amounts its iterations handed (shared).
  type :: removals_t
    logical        :: asRate = .false., shared = .false.
    type(handed_t) :: h
  end type removals_t

  !> IGLOO's smoke a MOSE block made, withheld from the gas and given to ICE (IGLOO_Smoke_Deliver): its amounts (mass,
  !  momentum (3), energy on the gas's datum); under handRates the rates set as ICE's source (and the number); the
  !  cells of the block's iterations where a smoke was born with no parcel's state.
  type :: smoke_t
    real(R8) :: amount(5) = 0._R8, rate(6) = 0._R8
    integer  :: nNoState  = 0
  end type smoke_t

  !> .true.: what lands between solver steps goes as rates; .false.: as amounts. Set by Setup_ICE_Coupling.
  logical, protected, public :: handRates = .false.

  !> Number of (cell, species) extractions clamped by the consumed-species floor guard.
  integer, protected, public :: nMassLimited = 0

  !> Inflow mass flow MOSExIGLOO handed IGLOO at its last call [kg/s], and the number of face records carrying it.
  real(R8), protected, public :: mdotHanded = 0._R8
  integer,  protected, public :: nHanded    = 0

  !> family -> parcel material id, built once after IGLOO setup (getSourceTerms needs it).
  integer, allocatable :: famMatID(:)

  !> Under-relaxation store: the residual IGLOOxMOSE last blended, per block.
  type :: res_block_t
    real(R8), allocatable :: res(:,:,:,:)
  end type res_block_t
  type(res_block_t), allocatable :: oldResidual(:)

  !> A value per cell of a block (interior cells).
  type :: vol_t
    real(R8), allocatable :: v(:,:,:)
  end type vol_t

  !> IGLOO's smoke on the cells of a coarse level, per block: r(6, channel, cell) the sum over the cell's sub-cells of
  !  each smoke channel's rate (Smoke_Packet), nostate(cell) the sub-cells, once per material, whose smoke took no
  !  parcel's state.
  type :: smoke_rate_t
    real(R8), allocatable :: r(:,:,:,:,:)
    integer,  allocatable :: nostate(:,:,:)
  end type smoke_rate_t

  !> The condensed solvers' sources restricted from the finest grid to the coarser level MOSE steps (held, 0: none),
  !  from the sources in force iglooAt and iceAt: IGLOO's rates, its smoke's and its Eulerian fields, ICE's per-volume
  !  rates, and the removals' channels on the coarse cells. ratio, the sub-cells per coarse cell along i, j, k, and
  !  vol, their volume, are those of the level they were set up for (level).
  type :: coarse_src_t
    integer :: level = 0, held = 0
    integer,               allocatable :: ratio(:,:)
    type(vol_t),           allocatable :: vol(:)
    type(obj_sourceblock), allocatable :: src(:)
    type(smoke_rate_t),    allocatable :: smk(:)
    type(obj_eulerblock),  allocatable :: eul(:,:)
    type(ICE_domain_type)              :: ice
    type(channel_t),       allocatable :: removed(:)
    type(source_state_t)               :: iglooAt, iceAt
  end type coarse_src_t
  type(coarse_src_t) :: coarse
  !> The finest grid lags a coarser level MOSE stepped since its last fill.
  logical :: fineStale = .false.
  !> The finest grid's boundary records hold MOSE's evaluation of its gas (an iteration there, or a fill).
  logical :: gasEvaluated = .false.

  !> A total over block b of a quantity restricted to a coarse level by who (IGLOO or ICE): the finest grid's, the
  !  coarse level's, the sum of its magnitudes on the finest grid. The last restriction's, for the report.
  type :: restricted_t
    character(len=8)  :: who  = ''
    character(len=32) :: what = ''
    integer  :: b = 0
    real(R8) :: onFine = 0._R8, onCoarse = 0._R8, magnitude = 0._R8
  end type restricted_t
  type(restricted_t), allocatable :: totals(:)
  integer :: nTotals = 0

  !> removed(g): what the burn-out of ICE family g removed (mass, momentum (3), energy, per unit volume), from ICE to
  !  the gas; the store of every family, though only a converting one hands it on. Under rates timed: the store over
  !  the time it was made in is the rate in force, kept while ICE does not step, as its source is. Under amounts each
  !  MOSE iteration takes its share. The store, its time and its rate are restored by the hand-off file of a restart.
  type(channel_t), allocatable, public :: removed(:)

  !> iceMean(g): over a time-accurate ICE block, family g's dt*source, dt*heat, dt*birth (its components, in this
  !  order) and dt, cell by cell (After_ICE_Step), made into the block's rates after it (ICE_Block_Mean) when it made
  !  meanSteps > 1 steps.
  type(channel_t), allocatable :: iceMean(:)
  integer :: meanSteps = 0

  !> smoke(c): the smoke of an IGLOO material (from), born in an ICE family (to): the smoke of its combustion record
  !  (trigger SMOKE_OF_BURN) or of its cap oxide handed at burnout (SMOKE_OF_CAP); kg per kg, the smoke's datum
  !  offset, the material's d-smoke. Its store: what the current MOSE block made (mass, momentum (3), energy on ICE's
  !  datum, number). smokeIGOn: one at least. nSmokeNoState: the cells of the block's iterations where a smoke was
  !  born with no parcel's state.
  integer, parameter :: SMOKE_OF_BURN = 1, SMOKE_OF_CAP = 2
  type(channel_t), allocatable :: smoke(:)
  logical, protected, public   :: smokeIGOn = .false.
  integer :: nSmokeNoState = 0

  character(len=*), parameter :: E_ = ' [ERROR] MI2: '

contains

  ! ---------------------------------------------------------------------------------------------
  ! Setup
  ! ---------------------------------------------------------------------------------------------

  !> ICE's side of the coupler, after ICE%setup: how what lands between steps goes (handRates), material -> gas
  !  species of ICE's phase from the file ATLAS writes, the smoke families, the datum offsets, and the removals'
  !  channels. ICE credits a material its energy at T_ref on its own datum: its table's, or cp*T.
  subroutine Setup_ICE_Coupling ( ICE_domain )
    use MOSE_Config_Types_m,    only: MOSE_time_scheme => obj_time_scheme
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ICE_phase_prefix
    use ICE_Lib_Properties,     only: mat_e, mat_cp_const, mat_lv
    use ICE_Config_Types_m,     only: COMB_BECKSTEAD, ICE_condensed => obj_condensed, ICE_time_scheme => obj_time_scheme
    use hydraMI2_Channels,      only: As_Rate
    use hydraMI2_Stoichiometry, only: Load_Species_Map, Report_Stoichiometry, Build_Datum_Offsets, Smoke_Record, &
                                      Set_Smoke_Datum, MAP_ICE, T_datum_ref, productsICE
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    character(len=64), allocatable :: matnames(:)
    logical,           allocatable :: evap(:), burn(:), beck(:), Lvconst(:)
    real(R8),          allocatable :: eref(:), Lvmat(:), qmat(:), Thot(:), ehot(:), beta(:)
    character(len=64) :: sname, sphase
    real(R8) :: snu, sh, dh
    integer  :: m, nm, gs

    handRates = As_Rate( MOSE_time_scheme%time_accurate, ICE_time_scheme%time_accurate )

    nm = size(ICE_condensed)
    allocate(matnames(nm), evap(nm), burn(nm), beck(nm), Lvconst(nm))
    allocate(eref(nm), Lvmat(nm), qmat(nm), Thot(nm), ehot(nm), beta(nm))
    do m = 1, nm
      associate ( mat => ICE_condensed(m) )
        matnames(m) = mat%name
        evap(m)     = mat%evapSelect > 0
        burn(m)     = mat%combSelect > 0
        beck(m)     = mat%combSelect == COMB_BECKSTEAD
        beta(m)     = merge(mat%betaPart, 0._R8, beck(m))
        ! where the material converts: its boiling temperature, or the one it burns from (Beckstead: T-ign)
        Thot(m)     = merge(mat%Tign, mat%Tboil, beck(m))
        if ( mat%cs_varies ) then
          eref(m) = mat_e( mat, T_datum_ref )
          ehot(m) = mat_e( mat, Thot(m) )
        else
          eref(m) = mat_cp_const( mat ) * T_datum_ref
          ehot(m) = mat_cp_const( mat ) * Thot(m)
        endif
        Lvmat(m)    = mat_lv( mat, T_datum_ref )
        Lvconst(m)  = .not. mat%use_lv
        qmat(m)     = mat%qComb
      end associate
    enddo

    call Load_Species_Map( MAP_ICE, ICE_phase_prefix, matnames, evap, burn )

    ! the smoke of ICE's burning materials, which ICE itself makes (cond_phase%birth): its family, kg per kg, datum
    do m = 1, nm
      if ( .not. Smoke_Record( MAP_ICE, m, sname, snu, sh, sphase ) ) cycle
      call Smoke_Family( trim(matnames(m)), 'INPUT/'//trim(ICE_phase_prefix)//'products.txt', sname, sphase, sh, gs, dh )
      if ( .not. ICE_condensed(m)%dSmoke > 0._R8 ) then
        write(*,'(A)') E_//trim(matnames(m))//' burns into the smoke '//trim(sname)//' and its row of INPUT/'// &
          trim(ICE_phase_prefix)//'models.txt gives no d-smoke, the diameter each smoke drop is born with'
        error stop 1
      endif
      ICE_condensed(m)%smokeFam = gs
      ICE_condensed(m)%nuSmoke  = snu
      call Set_Smoke_Datum( m, dh )
      write(*,'(A,I0,A,ES15.7,A,ES12.4,A,ES15.7,A)') ' [MI2] ICE material '//trim(matnames(m))//' burns into the '// &
        'smoke '//trim(sname)//' (family ', gs, '): ', snu, ' kg per kg burnt, d-smoke ', ICE_condensed(m)%dSmoke, &
        ' m, datum offset ', dh, ' J/kg'
    enddo

    call Report_Stoichiometry( MAP_ICE, matnames, nm )
    call Build_Datum_Offsets( MAP_ICE, matnames, nm, beck, beta, eref, Lvmat, qmat, Thot, ehot, Lvconst )

    ! the removals reach the gas as the rate of their store, timed by ICE's steps, or as amounts
    call Removed_Store_Setup( ICE_domain, timed = handRates .and. productsICE%active )
    removed(:)%asRate = handRates

  end subroutine Setup_ICE_Coupling


  !> The channels of ICE's removals, one per family (a restart may have opened them first): each store, and when timed
  !  its time and rate too.
  subroutine Removed_Store_Setup ( ICE_domain, timed )
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: ngroups
    use hydraMI2_Channels,    only: Open_Store
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    logical,               intent(in) :: timed
    integer :: b, g, dims(3, ICE_domain%nb)

    if (.not. allocated(removed)) then
      allocate(removed(ngroups))
      do g = 1, ngroups
        removed(g)%name = 'ice-removed'
        removed(g)%from = g
      enddo
    endif
    do b = 1, ICE_domain%nb
      dims(:,b) = ICE_domain%blk(b)%dim(1:3)
    enddo
    do g = 1, ngroups
      call Open_Store( removed(g), dims, 5, timed=timed, rated=timed )
    enddo

  end subroutine Removed_Store_Setup


  !> The removals' channels are timed (their rates kept by the coupler).
  logical function Removals_Timed ()
    implicit none

    Removals_Timed = .false.
    if (.not. allocated(removed)) return
    if (size(removed) == 0) return
    if (.not. allocated(removed(1)%blk)) return
    Removals_Timed = allocated(removed(1)%blk(1)%r)

  end function Removals_Timed


  !> IGLOO's side of the coupler, after IGLOO%setup: the family -> material map, material -> gas species of IGLOO's
  !  phase from the file ATLAS writes (with the oxide cap of a material that burns with one), the smoke of its burning
  !  materials, the datum offsets, and the under-relaxation store (allocated whether or not relaxation is on).
  subroutine Setup_IGLOO_Coupling ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nprim
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nfam, IGLOO_phase_prefix
    use IGLOO_Lib_Combustion,   only: nuOx
    use IGLOO_Lib_Evaporation,  only: IGLOO_COMB_BECKSTEAD => COMB_BECKSTEAD
    use hydraMI2_Stoichiometry, only: Load_Species_Map, Report_Stoichiometry, Build_Datum_Offsets, MAP_IGLOO, &
                                      T_datum_ref
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    character(len=20), allocatable :: matnames(:)
    logical,           allocatable :: evap(:), burn(:), handsCap(:), beck(:), Lvconst(:)
    integer,           allocatable :: oxide(:)
    real(R8),          allocatable :: href(:), Lvmat(:), qmat(:), Thot(:), hhot(:), beta(:)
    integer :: b, m, g, fam, nmat

    ! family -> material, from the famID IGLOO%setup assigned
    if (allocated(famMatID)) deallocate(famMatID)
    allocate(famMatID(1:max(nfam,1)))
    famMatID = 1
    do m = 1, size(IGLOO%material)
      do g = 1, IGLOO%material(m)%ngroups
        fam = IGLOO%material(m)%group(g)%famID
        if (fam >= 1 .and. fam <= nfam) famMatID(fam) = m
      enddo
    enddo

    ! per material: what it converts by, and the enthalpy IGLOO credits it at T_ref and where it converts (boiling,
    ! or T-ign on Beckstead's law) on the datum of its table: cp*T + hOff for constant cp, hTab interpolated otherwise
    nmat = size(IGLOO%material)
    allocate(matnames(nmat), evap(nmat), burn(nmat), handsCap(nmat), oxide(nmat), beck(nmat), Lvconst(nmat))
    allocate(href(nmat), Lvmat(nmat), qmat(nmat), Thot(nmat), hhot(nmat), beta(nmat))
    do m = 1, nmat
      associate ( mat => IGLOO%material(m) )
        matnames(m) = mat%matName
        evap(m)     = mat%evapSelect > 0
        burn(m)     = mat%combSelect > 0
        ! a material that burns with an oxide cap names the material its oxide is, and whether the cap ends in the gas
        oxide(m)    = merge(mat%prodID, 0, burn(m))
        handsCap(m) = oxide(m) > 0 .and. .not. mat%resOn
        beck(m)     = mat%combSelect == IGLOO_COMB_BECKSTEAD
        beta(m)     = merge(mat%betaPart, 0._R8, beck(m))
        Thot(m)     = merge(mat%Tign, mat%Tboil, beck(m))
        Lvmat(m)    = mat%Lv
        Lvconst(m)  = .not. mat%LvVariable
        qmat(m)     = mat%qComb
        if ( mat%cpVariable ) then
          href(m) = Table_Enthalpy( mat%hTab, lbound(mat%hTab,1), T_datum_ref )
          hhot(m) = Table_Enthalpy( mat%hTab, lbound(mat%hTab,1), Thot(m) )
        else
          href(m) = mat%cp * T_datum_ref + mat%hOff
          hhot(m) = mat%cp * Thot(m) + mat%hOff
        endif
        write(*,'(A,I0,A,ES14.7,A)') ' [MI2] IGLOO material ', m, ': table datum '//trim(mat%hDatum)//', hOff = ', &
          mat%hOff, ' J/kg'
      end associate
    enddo

    call Load_Species_Map( MAP_IGLOO, IGLOO_phase_prefix, matnames, evap, burn, oxide, handsCap, nuOx )
    call Setup_IGLOO_Smoke( MOSE_domain, IGLOO, matnames )
    call Report_Stoichiometry( MAP_IGLOO, matnames, nmat )
    call Build_Datum_Offsets( MAP_IGLOO, matnames, nmat, beck, beta, href, Lvmat, qmat, Thot, hhot, Lvconst )

    if (allocated(oldResidual)) deallocate(oldResidual)
    allocate(oldResidual(1:size(MOSE_domain%blk)))
    do b = 1, size(oldResidual)
      allocate(oldResidual(b)%res(1:nprim, 1:MOSE_domain%blk(b)%dim(1), 1:MOSE_domain%blk(b)%dim(2), &
                                         1:MOSE_domain%blk(b)%dim(3)), source=0._R8)
    enddo

  contains

    !> IGLOO's enthalpy table ht(lb:), indexed by the temperature in K, interpolated linearly at T (clamped to its
    !  range).
    real(R8) function Table_Enthalpy ( ht, lb, T )
      integer,  intent(in) :: lb
      real(R8), intent(in) :: ht(lb:), T
      integer :: it
      it = max(lbound(ht,1), min(int(T), ubound(ht,1)-1))
      Table_Enthalpy = ht(it) + (ht(it+1)-ht(it))*(T - real(it,R8))
    end function Table_Enthalpy

  end subroutine Setup_IGLOO_Coupling


  !> The smoke channels of IGLOO's burning materials (the tail of a combustion record, and of the release record of a
  !  cap oxide handed at burnout): ICE runs and holds its family, the coupling is not relaxed, the material gives
  !  d-smoke. One channel per route, in the order of the materials, then of the records.
  subroutine Setup_IGLOO_Smoke ( MOSE_domain, IGLOO, matnames )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use IGLOO_module,            only: obj_IGLOO
    use IGLOO_variables,         only: IGLOO_phase_prefix
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydraMI2_Channels,       only: Open_Store
    use hydraMI2_Stoichiometry,  only: Smoke_Record, MAP_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    character(len=*),       intent(in) :: matnames(:)
    character(len=64)  :: sname, sphase
    character(len=256) :: file
    character(len=:), allocatable :: who
    real(R8) :: snu, sh, dh
    integer  :: m, p, gs, n, b, c, nc, dims(3, size(MOSE_domain%blk))
    logical  :: has
    type(channel_t), allocatable :: found(:)

    allocate(found(2*size(matnames)))
    nc = 0
    file = 'INPUT/'//trim(IGLOO_phase_prefix)//'products.txt'
    do m = 1, size(matnames)
      who = 'the IGLOO material '//trim(matnames(m))
      do n = SMOKE_OF_BURN, SMOKE_OF_CAP
        if ( n == SMOKE_OF_BURN ) then
          has = Smoke_Record( MAP_IGLOO, m, sname, snu, sh, sphase )
        else
          p   = IGLOO%material(m)%prodID
          has = .false.
          if ( IGLOO%material(m)%combSelect > 0 .and. p > 0 .and. .not. IGLOO%material(m)%resOn ) &
            has = Smoke_Record( MAP_IGLOO, p, sname, snu, sh, sphase, release=.true. )
        endif
        if ( .not. has ) cycle
        if ( .not. hydra_MI2_obj%run_ICE ) then
          write(*,'(A)') E_//trim(file)//': '//who//' burns into the smoke '//trim(sname)//', and the smoke has no '// &
            'family: ICE is not run. Run ICE with '//trim(sname)//' in its phase, or burn the material into gas species only'
          error stop 1
        endif
        if ( hydra_MI2_obj%relaxation ) then
          write(*,'(A)') E_//who//' burns into the smoke '//trim(sname)//': relaxation = true blends two coupling '// &
            'steps of the gas''s share, and the smoke''s is not blended; set relaxation = false'
          error stop 1
        endif
        if ( .not. IGLOO%material(m)%dSmoke > 0._R8 ) then
          write(*,'(A)') E_//who//' burns into the smoke '//trim(sname)//' and its row of INPUT/'// &
            trim(IGLOO_phase_prefix)//'models.txt gives no d-smoke, the diameter each smoke drop is born with'
          error stop 1
        endif
        call Smoke_Family( who, file, sname, sphase, sh, gs, dh )
        nc = nc + 1
        found(nc)%name    = 'igloo-smoke'
        found(nc)%from    = m
        found(nc)%to      = gs
        found(nc)%trigger = n
        found(nc)%nu      = snu
        found(nc)%dh      = dh
        found(nc)%d       = IGLOO%material(m)%dSmoke
        found(nc)%asRate  = handRates
        write(*,'(A,I0,A,ES15.7,A,ES12.4,A,ES15.7,A)') ' [MI2] IGLOO material '//trim(matnames(m))//' '// &
          trim(merge('burns into' // repeat(' ', 17), 'hands its cap at burnout to', n == SMOKE_OF_BURN))// &
          ' the smoke '//trim(sname)//' (ICE family ', gs, '): ', snu, ' kg per kg, d-smoke ', &
          IGLOO%material(m)%dSmoke, ' m, datum offset ', dh, ' J/kg'
      enddo
    enddo
    allocate(smoke(nc))
    do c = 1, nc
      smoke(c) = found(c)
    enddo
    smokeIGOn = nc > 0
    ! each channel's store: what a MOSE block's iterations make (mass, momentum (3), energy, number)
    do b = 1, size(MOSE_domain%blk)
      dims(:,b) = MOSE_domain%blk(b)%dim(1:3)
    enddo
    do c = 1, size(smoke)
      call Open_Store( smoke(c), dims, 6, timed=.false., rated=.false. )
    enddo
    if ( smokeIGOn .and. handRates ) &
      write(*,'(A)') ' [MI2] IGLOO''s smoke reaches ICE as rates: a source of each smoke family over ICE''s own steps'
    if ( smokeIGOn .and. .not. handRates ) &
      write(*,'(A)') ' [MI2] IGLOO''s smoke reaches ICE as amounts: what each MOSE block made, after the block'

  end subroutine Setup_IGLOO_Smoke


  !> The ICE family a smoke is born in, for ICE's materials and IGLOO's alike: the smoke material is a material of
  !  ICE's phase (the record's phase) with exactly one family, of the MK closure, that neither evaporates nor burns.
  !  Returns the family and the smoke's datum offset: its formation enthalpy sh less ICE's energy of it at T_datum_ref.
  subroutine Smoke_Family ( owner, file, sname, sphase, sh, gs, dh )
    use ICE_Global_m,           only: ICE_phase_prefix, mat_of, ngroups
    use ICE_Config_Types_m,     only: ICE_condensed => obj_condensed, ICE_time_scheme => obj_time_scheme
    use ICE_Lib_Properties,     only: mat_e, mat_cp_const
    use hydraMI2_Stoichiometry, only: T_datum_ref
    implicit none
    character(len=*), intent(in)  :: owner, file, sname, sphase
    real(R8),         intent(in)  :: sh
    integer,          intent(out) :: gs
    real(R8),         intent(out) :: dh
    integer :: g, ms, nfam

    ms = 0
    do g = 1, size(ICE_condensed)
      if ( trim(ICE_condensed(g)%name) == trim(sname) ) ms = g
    enddo
    if ( ms == 0 ) then
      write(*,'(A)') E_//trim(file)//': '//owner//' burns into the smoke '//trim(sname)//', which is not a '// &
        'material of ICE''s phase: declare it in INPUT/'//trim(ICE_phase_prefix)//'phase.txt with one population '// &
        'and its family with closure = MK'
      error stop 1
    endif
    if ( trim(sphase) /= trim(ICE_phase_prefix) ) then
      write(*,'(A)') E_//'the smoke of '//owner//' is born in the phase '//trim(sphase)//', and ICE runs the phase '// &
        trim(ICE_phase_prefix)//': run ATLAS GPB again'
      error stop 1
    endif
    nfam = count(mat_of(1:ngroups) == ms)
    gs   = findloc(mat_of(1:ngroups), ms, dim=1, back=.true.)
    if ( nfam /= 1 ) then
      write(*,'(A,I0,A)') E_//'the smoke '//trim(sname)//' of '//owner//' has ', nfam, &
        ' families in ICE: it is born in one, give its material one population'
      error stop 1
    endif
    if ( trim(ICE_time_scheme%model(gs)) /= 'MK' ) then
      write(*,'(A,I0,A)') E_//'the smoke '//trim(sname)//' is born in ICE family ', gs, ', whose closure is '// &
        trim(ICE_time_scheme%model(gs))//': set closure = MK (a smoke takes mass and number at a fixed diameter)'
      error stop 1
    endif
    if ( ICE_condensed(ms)%evapSelect > 0 .or. ICE_condensed(ms)%combSelect > 0 ) then
      write(*,'(A)') E_//'the smoke '//trim(sname)//' evaporates or burns in ICE: a smoke does neither'
      error stop 1
    endif
    if ( ICE_condensed(ms)%cs_varies ) then
      dh = sh - mat_e( ICE_condensed(ms), T_datum_ref )
    else
      dh = sh - mat_cp_const( ICE_condensed(ms) ) * T_datum_ref
    endif

  end subroutine Smoke_Family


  !> Select IGLOO's coupling mode for cycle n_cycle: IGLOO-iter = 2 -> Eulerian (eulSwitch), 1 -> source (srcSwitch).
  !  On the last cycle that runs IGLOO both switches are on for the writeout.
  subroutine Set_IGLOO_Switches ( IGLOO, n_cycle, eulerian )
    use IGLOO_module,            only: obj_IGLOO
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    type(obj_IGLOO), intent(inout) :: IGLOO
    integer,         intent(in)    :: n_cycle
    logical,         intent(out)   :: eulerian
    logical :: last

    eulerian        = hydra_MI2_obj%IGLOO_iter(n_cycle) == 2
    last            = all( hydra_MI2_obj%IGLOO_iter(n_cycle+1:) == 0 )
    IGLOO%eulSwitch = eulerian .or. last
    IGLOO%srcSwitch = .not. eulerian .or. last

  end subroutine Set_IGLOO_Switches


  ! ---------------------------------------------------------------------------------------------
  ! MOSE -> condensed solvers: the gas state
  ! ---------------------------------------------------------------------------------------------

  !> Derived gas state (rhoi, vel, rho, R, gam, T, mu, kappa) from one MOSE primitive column.
  subroutine Gas_State_From_Prim ( prim, rhoi, vel, rho, R, gam, T, mu, kappa )
    use MOSE_Global_m,           only: nsc, nu, nw, np
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    implicit none
    real(R8), intent(in)  :: prim(:)
    real(R8), intent(out) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa

    rhoi = prim(1:nsc)
    vel  = prim(nu:nw)
    call co_rotot_Rtot ( rhoi, rho, R )
    gam = f_gamma ( rhoi, prim(np), rho, R )
    T   = EOS ( prim(np), rho, R )
    call co_k_mi_lam_Wilke ( rhoi, rho, T, mu, kappa )

  end subroutine Gas_State_From_Prim


  !> Build the ORION gas field IGLOO ingests at setup from MOSE's domain: x, y, z, rho(1..nsc),
  !  u, v, w, T, mil, kl, gam, R (import_gas maps varname v to vars(v-3)).
  subroutine Gas_Field_For_IGLOO ( MOSE_domain, MOSE_IOfield, gas_for_igloo )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc
    use Lib_ORION_data,        only: orion_data, copyORION
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    type(orion_data),       intent(in)  :: MOSE_IOfield
    type(orion_data),       intent(out) :: gas_for_igloo
    integer  :: b, i, j, k, v, nv
    real(R8) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa

    ! start from MOSE's field: mesh, block names and format options
    call copyORION ( MOSE_IOfield, gas_for_igloo )

    nv = nsc + 8   ! rho(1..nsc), u, v, w, T, mil, kl, gam, R

    if (allocated(gas_for_igloo%varnames)) deallocate(gas_for_igloo%varnames)
    allocate(gas_for_igloo%varnames(3+nv))
    gas_for_igloo%varnames(1:3) = [ 'x', 'y', 'z' ]
    do v = 1, nsc
      write(gas_for_igloo%varnames(3+v),'(A,I0,A)') 'rho(', v, ')'
    enddo
    gas_for_igloo%varnames(4+nsc:) = [character(len=3) :: 'u', 'v', 'w', 'T', 'mil', 'kl', 'gam', 'R']

    do b = 1, size(MOSE_domain%blk)
      if (allocated(gas_for_igloo%block(b)%vars)) deallocate(gas_for_igloo%block(b)%vars)
      allocate(gas_for_igloo%block(b)%vars(1:nv, 1:MOSE_domain%blk(b)%dim(1), 1:MOSE_domain%blk(b)%dim(2), &
                                                 1:MOSE_domain%blk(b)%dim(3)))

      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) PRIVATE(rhoi,vel,rho,kappa,mu,R,gam,T,i,j,k)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), rhoi, vel, rho, R, gam, T, mu, kappa )
        gas_for_igloo%block(b)%vars(1:nsc  ,i,j,k) = rhoi
        gas_for_igloo%block(b)%vars(nsc+1: ,i,j,k) = [vel, T, mu, kappa, gam, R]
      enddo ; enddo ; enddo
      !$OMP END PARALLEL DO
    enddo

    write(*,'(A,I0,A,I0,A)') ' [MI2] built gas field for IGLOO: ', nsc, &
      ' species, ', nv, ' variables (rho(i), u, v, w, T, mil, kl, gam, R)'

  end subroutine Gas_Field_For_IGLOO


  !> Push the MOSE gas state (interior only: IGLOO rebuilds its own ghost ring in solve) into IGLOO's gasblock and
  !  re-impose the gas inflow mass flow mdotGas [kg/s] on the inlet faces, every sweep; mdotHanded and nHanded say
  !  how much and on how many face records. Opens its own parallel regions.
  subroutine MOSExIGLOO ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc, gc, MOSE_phase_prefix
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(obj_IGLOO),        intent(inout) :: IGLOO
    integer  :: b, i, j, k, f, m, n, lower, upper, nSum, noFlow
    real(R8) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa
    real(R8) :: mSum, mdot

    upper = 0
    mdotHanded = 0._R8
    nHanded    = 0
    do b = 1, size(MOSE_domain%blk)

      associate ( gb => IGLOO%gasblock(b) )
      ! the gasblock, ghost ring included under gas-order = 2, lies within MOSE's ghosted extent
      if ( any(lbound(gb%temperature) < 1-gc) .or. any(ubound(gb%temperature) > MOSE_domain%blk(b)%dim(1:3)+gc) ) then
        write(*,'(A,I0,A)') E_//'IGLOO gasblock ', b, ' extends beyond MOSE''s ghosted range: the two meshes disagree'
        error stop 1
      endif

      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) PRIVATE(rhoi,vel,rho,kappa,mu,R,gam,T,i,j,k)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), rhoi, vel, rho, R, gam, T, mu, kappa )
        gb%density    (1:nsc,i,j,k) = rhoi
        gb%velocity   (1:3  ,i,j,k) = vel
        gb%temperature(      i,j,k) = T
        gb%gam        (      i,j,k) = gam
        gb%R          (      i,j,k) = R
        gb%mil        (      i,j,k) = mu
        gb%kl         (      i,j,k) = kappa
      enddo ; enddo ; enddo
      !$OMP END PARALLEL DO
      end associate

      ! the gas inflow mass flow on the inlet faces (MOSE_domain%bc is flat over blocks and faces, in the order of
      ! MOSE's file); other records carry no mass flow and leave IGLOO's value alone. A solid-propellant grain in its
      ! inflow mode (ATLAS 502, grain-wall = F: a 404 with srm_inflow) hands the gas mass flux MOSE last computed at
      ! the face, bc%g_srm = mdot_ref (p/p_ref)^n at the boundary cell's pressure (its mdot_ref already carries
      ! 1 - krho_tot, ATLAS): IGLOO then injects krho/(1 - krho_tot) times it.
      do f = 1, 6
        lower = upper + 1
        upper = upper + MOSE_domain%n_bf(b,f)
        do i = lower, upper
          if (MOSE_domain%bc(i)%b == b .and. MOSE_domain%bc(i)%f == f) cycle
          write(*,'(A,I0,3A,I0,A,I0,A,I0,A,I0,A)') E_//'record ', i, ' of INPUT/', trim(MOSE_phase_prefix), &
            'bc.txt is of block ', MOSE_domain%bc(i)%b, ' face ', MOSE_domain%bc(i)%f, ' where the coupler expects block ', &
            b, ' face ', f, ': the gas boundary records must be grouped by block, then by face (1 to 6), in that order'
          error stop 1
        enddo
        mSum = 0._R8; nSum = 0; noFlow = 0
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,m,n,mdot) REDUCTION(+:mSum,nSum) REDUCTION(max:noFlow)
        do i = lower, upper
          call IGLOO%geoblock(b)%ijk2fmn ( MOSE_domain%bc(i)%i, MOSE_domain%bc(i)%j, &
                                           MOSE_domain%bc(i)%k, f, m, n )
          if (m<=0 .or. n<=0) cycle
          select case (MOSE_domain%bc(i)%type)
          case (401, 402, 405)             ! mdot is the face mass flow MOSE computes [kg/s] along the solver normal:
            mdot = abs(MOSE_domain%bc(i)%mdot)   ! negative for an inflow on faces 2, 4 and 6
          case (406, 407)
            mdot = MOSE_domain%bc(i)%mdot
          case (403, 404, 420)             ! mdot is the prescribed mass flux [kg/(s m2)], an inflow whatever its sign,
            if (MOSE_domain%bc(i)%srm_inflow) then   ! or for a grain the burn law's at the cell pressure (g_srm >= 0)
              mdot = abs(MOSE_domain%bc(i)%g_srm) * IGLOO%geoblock(b)%face(f)%cell(m,n)%area
            else
              mdot = abs(MOSE_domain%bc(i)%mdot) * IGLOO%geoblock(b)%face(f)%cell(m,n)%area
            end if
          case default                     ! no mass flow in the record: a particle inlet that scales with it is refused below
            if (IGLOO%geoblock(b)%face(f)%cell(m,n)%bcdef == 401) noFlow = max(noFlow, i)
            cycle
          end select
          IGLOO%geoblock(b)%face(f)%cell(m,n)%mdotGas = mdot
          mSum = mSum + mdot
          nSum = nSum + 1
        enddo
        !$OMP END PARALLEL DO
        if (noFlow > 0) then
          write(*,'(A,I0,A,I0,A,I0,A)') E_//'block ', b, ' face ', f, ': particles injected in proportion to '// &
            'the gas mass flow (particle BC 401) on a gas boundary of kind ', MOSE_domain%bc(noFlow)%type, &
            ', whose record holds no mass flow'
          error stop 1
        endif
        mdotHanded = mdotHanded + mSum
        nHanded    = nHanded + nSum
      enddo

    enddo

  end subroutine MOSExIGLOO


  !> Push the MOSE gas state into ICE's gas_phase: prim = (rho, u, v, w, T), interior only, plus
  !  R, gam, mu, k. Opens its own parallel region.
  subroutine MOSExICE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    integer  :: b, i, j, k
    real(R8) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa

    do b = 1, size(MOSE_domain%blk)
      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) PRIVATE(rhoi,vel,rho,kappa,mu,R,gam,T,i,j,k)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), rhoi, vel, rho, R, gam, T, mu, kappa )
        ICE_domain%blk(b)%gas_phase%prim(1:5,i,j,k) = [rho, vel, T]
        ICE_domain%blk(b)%gas_phase%R   (    i,j,k) = R
        ICE_domain%blk(b)%gas_phase%gam (    i,j,k) = gam
        ICE_domain%blk(b)%gas_phase%mu  (    i,j,k) = mu
        ICE_domain%blk(b)%gas_phase%k   (    i,j,k) = kappa
      enddo ; enddo ; enddo
      !$OMP END PARALLEL DO
    enddo

  end subroutine MOSExICE


  !> After MOSE stepped a coarser level (a start on its coarse grids): the finest grid, which IGLOO and ICE see, filled
  !  from it, its state prolongated level by level, then evaluated. Nothing when the finest grid does not lag.
  subroutine Fill_Finest_Grid ( MOSE_simulation )
    use MOSE_Advanced_Types_m, only: MOSE_simulation_type
    use MOSE_Config_Types_m,   only: MOSE_multigrid => obj_multigrid
    use MOSE_Mod_Multigrid,    only: Prolongation
    implicit none
    type(MOSE_simulation_type), intent(inout) :: MOSE_simulation
    integer :: m

    if ( .not. fineStale ) return
    do m = MOSE_multigrid%MG_level, 2, -1
      call Prolongation( Fine=MOSE_simulation%domain(m-1), Coarse=MOSE_simulation%domain(m) )
    enddo
    call Evaluate_Finest_Grid( MOSE_simulation%domain(1) )
    fineStale = .false.

  end subroutine Fill_Finest_Grid


  !> Before IGLOO's hand-off: the finest grid evaluated once when neither a MOSE iteration there nor a fill has done
  !  it in this run (IGLOO solves before MOSE's first iteration): its boundary records hold the inflow mass flow
  !  MOSExIGLOO hands IGLOO only from MOSE's evaluation of the gas.
  subroutine Evaluate_Gas_Once ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain

    if ( gasEvaluated ) return
    call Evaluate_Finest_Grid( MOSE_domain )

  end subroutine Evaluate_Gas_Once


  !> What a MOSE iteration computes on its finest grid, MOSE_domain, before its update (derived variables, dtlocal,
  !  ghosts, fluxes, boundary records: the inflow mass flow IGLOO takes), without the condensed sources: no update, no
  !  time step. Opens its own parallel regions.
  subroutine Evaluate_Finest_Grid ( MOSE_domain )
    use MOSE_Advanced_Types_m,    only: MOSE_domain_type
    use MOSE_Config_Types_m,      only: MOSE_time_scheme => obj_time_scheme
    use MOSE_Lib_Preconditioning, only: Update_Derived_Variables
    use MOSE_Mod_dt,              only: Compute_dt
    use MOSE_Mod_Explicit,        only: Assemble_Residual
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain

    !$omp parallel
    call Update_Derived_Variables( MOSE_domain )
    !$omp end parallel
    MOSE_domain%dtglobal = 1d5
    call Compute_dt( MOSE_domain, MOSE_time_scheme%cfl, MOSE_time_scheme%vnn, MOSE_time_scheme%rampa_cfl_iter )
    !$omp parallel
    call Assemble_Residual( MOSE_domain, No_Source )
    !$omp end parallel
    gasEvaluated = .true.

  end subroutine Evaluate_Finest_Grid


  !> The external routine of the finest grid's residual when the coupler evaluates it: no source.
  subroutine No_Source ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
  end subroutine No_Source


  ! ---------------------------------------------------------------------------------------------
  ! Condensed solvers -> MOSE: the sources, from MOSE's per-RK-stage callback. Inside MOSE's parallel region:
  ! orphaned !$OMP DO only, no ASSOCIATE, and cell kernels that write only their arguments.
  ! ---------------------------------------------------------------------------------------------

  !> Inject IGLOO's reaction on the gas into MOSE's residual (R = R - S): on MOSE's finest grid from IGLOO's own
  !  arrays, on a coarser level from those restricted to it. iter of iMax is the MOSE iteration within its block and
  !  relCoeff the cycle's under-relaxation coefficient (1: the source is injected as it is).
  subroutine IGLOOxMOSE ( MOSE_domain, IGLOO, eulerian, iter, iMax, relCoeff )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: eulerian
    integer,                intent(in)    :: iter, iMax
    real(R8),               intent(in)    :: relCoeff

    if (MOSE_domain%mg_level == 1) then
      call IGLOO_Residual( MOSE_domain, IGLOO, IGLOO%source, IGLOO%euler, eulerian, iter, iMax, relCoeff )
    else
      call IGLOO_Residual( MOSE_domain, IGLOO, coarse%src, coarse%eul, eulerian, iter, iMax, relCoeff, coarse%vol, &
                           coarse%smk )
    endif

  end subroutine IGLOOxMOSE


  !> IGLOO's reaction on the gas into MOSE's residual, from the arrays given on the cells of MOSE_domain: src the
  !  source mode's rates, eul the Eulerian fields; when given, volE the volume the Eulerian branch's per-volume source
  !  is taken over (else the cell's own) and smk the smoke's rates (else those of the cell's column).
  !    source branch : source{Mass,Mom,En} are absolute rates; their mass through IGLOO_Mass_To_Gas, under the
  !                    consumed-species floor guard
  !    euler  branch : getSourceTerms is per unit volume (*vol); momentum and energy only
  subroutine IGLOO_Residual ( MOSE_domain, IGLOO, src, eul, eulerian, iter, iMax, relCoeff, volE, smk )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nu, nw, np, nprim
    use MOSE_Config_Types_m,    only: obj_chemistry
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nfam
    use hydraMI2_Stoichiometry, only: productsIGLOO
    implicit none
    type(MOSE_domain_type),             intent(inout) :: MOSE_domain
    type(obj_IGLOO),                    intent(in)    :: IGLOO
    type(obj_sourceblock), allocatable, intent(in)    :: src(:)
    type(obj_eulerblock),  allocatable, intent(in)    :: eul(:,:)
    logical,                            intent(in)    :: eulerian
    integer,                            intent(in)    :: iter, iMax
    real(R8),                           intent(in)    :: relCoeff
    type(vol_t),          optional,     intent(in)    :: volE(:)
    type(smoke_rate_t),   optional,     intent(in)    :: smk(:)
    logical  :: blend
    integer  :: b, i, j, k, m
    ! res mirrors blk%R(1:nprim): species 1:nsc, momentum nu:nw, energy np, then soot/passive/RANS (left at zero)
    real(R8) :: res(nprim), sources(4), datum, vE
    real(R8) :: rhoi(nsc), rho, u(3), T, R, gam, kappa, mu, coeff

    call Relaxation_Weight ( iter, iMax, relCoeff, coeff, blend )

    do b = 1, size(MOSE_domain%blk)

      !$OMP DO COLLAPSE(3) PRIVATE(res,sources,datum,rhoi,rho,u,T,R,gam,kappa,mu,i,j,k,m,vE)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)

        res = 0._R8
        if (eulerian) then
          call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), rhoi, u, rho, R, gam, T, mu, kappa )
          sources = 0._R8
          do m = 1, nfam
            sources = sources + IGLOO%getSourceTerms( [rho,u,T,gam,R,mu,kappa],      &
                                    eul(b,m)%velocity(:,i,j,k),                      &
                                    eul(b,m)%temperature(i,j,k),                     &
                                    eul(b,m)%density(i,j,k),                         &
                                    eul(b,m)%np(i,j,k), famMatID(m) )
          enddo
          res(nu:nw) = - sources(1:3)
          res(np)    = - sources(4)
          vE = MOSE_domain%blk(b)%vol(i,j,k)
          if (present(volE)) vE = volE(b)%v(i,j,k)
          res = res * vE
        else
          res(nu:nw) = src(b)%sourceMom(:,i,j,k)
          res(np)    = src(b)%sourceEn(i,j,k)
          if (productsIGLOO%active) then
            if (present(smk) .and. smokeIGOn) then
              call IGLOO_Mass_To_Gas ( src(b)%sourceMass(:,i,j,k), res, datum, smk(b)%r(:,:,i,j,k) )
            else
              call IGLOO_Mass_To_Gas ( src(b)%sourceMass(:,i,j,k), res, datum )
            endif
            call Floor_Guard ( res(1:nsc), MOSE_domain%blk(b)%P(1:nsc,i,j,k), MOSE_domain%blk(b)%vol(i,j,k), &
                               MOSE_domain%blk(b)%dtlocal(i,j,k) * obj_chemistry%strangcoeff )
          endif
        endif

        if (blend) then
          res(1:np) = coeff*res(1:np) + (1._R8-coeff)*oldResidual(b)%res(1:np,i,j,k)
          oldResidual(b)%res(:,i,j,k) = res
        endif

        MOSE_domain%blk(b)%R(:,i,j,k) = MOSE_domain%blk(b)%R(:,i,j,k) - res

      enddo ; enddo ; enddo
      !$OMP END DO

    enddo

  end subroutine IGLOO_Residual


  !> Inject ICE's reaction on the gas into MOSE's residual: on MOSE's finest grid from ICE's own arrays and the
  !  removals' channels, on a coarser level from those restricted to it.
  subroutine ICExMOSE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(ICE_domain_type),  intent(in)    :: ICE_domain

    if (MOSE_domain%mg_level == 1) then
      call ICE_Residual( MOSE_domain, ICE_domain, removed )
    else
      call ICE_Residual( MOSE_domain, coarse%ice, coarse%removed )
    endif

  end subroutine ICExMOSE


  !> ICE's reaction on the gas into MOSE's residual, from the arrays given on the cells of MOSE_domain: what each
  !  family's source in ICE_domain hands the gas (ICE_Cell) and the removals of the channels rem (Removal_Cell: a
  !  steady run, the rate of the store; a time-accurate one, this MOSE iteration's share of it as a rate over its step,
  !  which each stage and both Strang halves integrate to rate*dtlocal), under the consumed-species floor guard, times
  !  the cell volume.
  subroutine ICE_Residual ( MOSE_domain, ICE_domain, rem )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nprim
    use MOSE_Config_Types_m,    only: obj_chemistry
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(ICE_domain_type),  intent(in)    :: ICE_domain
    type(channel_t),        intent(in)    :: rem(:)
    integer  :: b, g, i, j, k, imat, nshare
    real(R8) :: r(nprim), dtIter, datum, gr(5)
    logical  :: mapped, consumes

    ! the families whose material consumes a gas species share what a cell holds of it
    nshare = count([(Consumes_Gas(mat_of(g)), g = 1, ngroups)])

    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups

        imat     = mat_of(g)
        mapped   = Converts(imat)
        consumes = Consumes_Gas(imat)

        !$OMP DO COLLAPSE(3) PRIVATE(r,i,j,k,dtIter,datum,gr)
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)

          r = 0._R8
          call ICE_Cell ( imat, mapped, ICE_domain%blk(b)%cond_phase(g)%source(:,i,j,k), nbase(g)-1, &
                          ICE_domain%blk(b)%cond_phase(g)%birth(:,i,j,k), r, datum )
          if (mapped) then
            if (rem(g)%asRate) then
              if (allocated(rem(g)%blk(b)%r)) then
                gr = rem(g)%blk(b)%r(:,i,j,k)
                call Removal_Cell ( imat, gr, r, datum )
              endif
            elseif (rem(g)%share > 0._R8) then
              dtIter = MOSE_domain%blk(b)%dtlocal(i,j,k)
              if (dtIter > 0._R8) then
                gr = rem(g)%blk(b)%a(:,i,j,k)*rem(g)%share/dtIter
                call Removal_Cell ( imat, gr, r, datum )
              endif
            endif
            ! a consumed species leaves at most this family's share of what the cell holds (per unit volume)
            if (consumes) &
              call Floor_Guard ( r(1:nsc), MOSE_domain%blk(b)%P(1:nsc,i,j,k), 1._R8,          &
                                 MOSE_domain%blk(b)%dtlocal(i,j,k) * obj_chemistry%strangcoeff * real(nshare, R8) )
          endif

          MOSE_domain%blk(b)%R(:,i,j,k) = MOSE_domain%blk(b)%R(:,i,j,k) - r * MOSE_domain%blk(b)%vol(i,j,k)

        enddo ; enddo ; enddo
        !$OMP END DO

      enddo
    enddo

  contains

    !> ICE material imat converts and consumes a gas species.
    logical function Consumes_Gas ( imat )
      integer, intent(in) :: imat
      Consumes_Gas = .false.
      if (Converts(imat)) Consumes_Gas = any(productsICE%nu(imat,1:nsc) < 0._R8)
    end function Consumes_Gas

  end subroutine ICE_Residual


  !> ICE material imat converts: it has a record in the products file, whose mass the gas receives.
  logical function Converts ( imat )
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    integer, intent(in) :: imat

    Converts = .false.
    if (productsICE%active) Converts = productsICE%mapped(imat)

  end function Converts


  !> What IGLOO's mass source in one cell (col, IGLOO%source(b)%sourceMass(:,i,j,k)) hands the gas, added to r (the
  !  residual's layout) as absolute rates: each material through its row (a capped material through its three, with
  !  Cap_Rates) with its datum offset, less the smoke born in ICE: the rates of each smoke channel smk(:, channel) when
  !  given (a coarse cell's, summed over its sub-cells), else Smoke_Packet's of col. datum: the offsets' part of what
  !  it adds to the energy. Thread-safe.
  subroutine IGLOO_Mass_To_Gas ( col, r, datum, smk )
    use MOSE_Global_m,          only: nsc, nu, nw, np
    use IGLOO_variables,        only: nm
    use hydraMI2_Stoichiometry, only: productsIGLOO, Cap_Rates
    implicit none
    real(R8),           intent(in)    :: col(:)
    real(R8),           intent(inout) :: r(:)
    real(R8),           intent(out)   :: datum
    real(R8), optional, intent(in)    :: smk(:,:)
    integer  :: imat, isp, c
    real(R8) :: dm, cr(4), sr(6)

    datum = 0._R8
    do imat = 1, nm
      dm = col(imat)
      ! an oxide cap: the material burnt into gas species, the oxide the cap gained, the cap handed to the gas
      if (productsIGLOO%capped(imat)) then
        cr = Cap_Rates( imat, dm, col(nm+2*imat-1), col(nm+2*imat) )
        do isp = 1, nsc
          r(isp) = r(isp) + productsIGLOO%nu(imat,isp)*cr(2) + productsIGLOO%nuCap(imat,isp)*cr(3) &
                          + productsIGLOO%nuRel(imat,isp)*cr(4)
        enddo
        r(np) = r(np) + productsIGLOO%dh(imat)*cr(1) - productsIGLOO%dhOxide(imat)*(cr(3) - cr(4))
        datum = datum + productsIGLOO%dh(imat)*cr(1) - productsIGLOO%dhOxide(imat)*(cr(3) - cr(4))
        cycle
      endif
      if (dm == 0._R8) cycle
      do isp = 1, nsc
        if (productsIGLOO%nu(imat,isp) == 0._R8) cycle
        r(isp) = r(isp) + productsIGLOO%nu(imat,isp)*dm
      enddo
      r(np) = r(np) + productsIGLOO%dh(imat)*dm
      datum = datum + productsIGLOO%dh(imat)*dm
    enddo

    ! the smoke of a material that makes one is born in its ICE family: the gas receives IGLOO's books less it
    if (.not. smokeIGOn) return
    do c = 1, size(smoke)
      if (present(smk)) then
        sr = smk(:,c)
      else
        call Smoke_Packet( smoke(c), col, sr )
      endif
      r(nu:nw) = r(nu:nw) - sr(2:4)
      r(np)    = r(np) - sr(5) - smoke(c)%dh*sr(1)
      datum    = datum - smoke(c)%dh*sr(1)
    enddo

  end subroutine IGLOO_Mass_To_Gas


  !> What the source of one ICE family in one cell hands the gas, per unit volume, into r (the residual's layout):
  !  -source (src = cond_phase%source(:,i,j,k)), the momentum src(2:4) and the energy src(ie); for a material that
  !  converts (mapped) also the mass -src(1) through its row with its datum offset, less the smoke the family makes
  !  (birth = cond_phase%birth(:,i,j,k)). datum: the offsets' part of the energy. Thread-safe.
  subroutine ICE_Cell ( imat, mapped, src, ie, birth, r, datum )
    use MOSE_Global_m,          only: nsc, nu, nw, np
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    integer,  intent(in)    :: imat, ie
    logical,  intent(in)    :: mapped
    real(R8), intent(in)    :: src(:), birth(:)
    real(R8), intent(inout) :: r(:)
    real(R8), intent(out)   :: datum
    real(R8) :: dm

    r(nu:nw) = - src(2:4)
    r(np)    = - src(ie)
    datum    = 0._R8
    if (.not. mapped) return

    dm = - src(1)
    r(1:nsc) = productsICE%nu(imat,1:nsc)*dm
    r(np)    = r(np) + productsICE%dh(imat)*dm
    ! the smoke the family makes is not the gas's: its momentum, its energy and its datum offset stay out
    r(nu:nw) = r(nu:nw) - birth(2:4)
    r(np)    = r(np) - birth(5) - productsICE%dhSmoke(imat)*birth(1)
    datum    = productsICE%dh(imat)*dm - productsICE%dhSmoke(imat)*birth(1)

  end subroutine ICE_Cell


  !> What a rate g of ICE's removals of material imat (mass, momentum (3), energy, per unit volume) hands the gas,
  !  added to r (the residual's layout): the mass through its row with the removal's datum offset (the material's,
  !  less that of the smoke it made), which is added to datum too. Thread-safe.
  subroutine Removal_Cell ( imat, g, r, datum )
    use MOSE_Global_m,          only: nsc, nu, nw, np
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    integer,  intent(in)    :: imat
    real(R8), intent(in)    :: g(5)
    real(R8), intent(inout) :: r(:), datum
    integer  :: isp
    real(R8) :: dh

    dh = productsICE%dh(imat) - productsICE%nuSmoke(imat)*productsICE%dhSmoke(imat)
    do isp = 1, nsc
      r(isp) = r(isp) + productsICE%nu(imat,isp)*g(1)
    enddo
    r(nu:nw) = r(nu:nw) + g(2:4)
    r(np)    = r(np) + g(5) + dh*g(1)
    datum    = datum + dh*g(1)

  end subroutine Removal_Cell


  !> Under-relaxation of IGLOO's source over a MOSE block of iMax iterations: the weight coeff of the new source
  !  against the stored one, and whether to blend at all. relCoeff outside (0,1) is no relaxation. A
  !  single-iteration block blends with relCoeff; a longer one ramps the weight up over its first iMax*(1-relCoeff)
  !  iterations and injects the source as it is afterwards.
  subroutine Relaxation_Weight ( iter, iMax, relCoeff, coeff, blend )
    implicit none
    integer,  intent(in)  :: iter, iMax
    real(R8), intent(in)  :: relCoeff
    real(R8), intent(out) :: coeff
    logical,  intent(out) :: blend
    real(R8), parameter :: toll = 1.e-15_R8

    blend = ( relCoeff > toll .and. relCoeff < (1._R8-toll) )
    if (iter < iMax*(1._R8-relCoeff) .and. blend) then
      coeff = real(iter,R8) / real(iMax,R8)
      coeff = coeff*(1._R8/(1._R8-relCoeff) - exp(-(coeff/relCoeff))**2)
      coeff = max(0._R8, min(1._R8, coeff))
    elseif ((iMax==1) .and. blend) then
      coeff = relCoeff
    else
      coeff = 1._R8
      blend = .false.
    endif

  end subroutine Relaxation_Weight


  !> Consumed-species floor guard: cap the extraction RATE of every species at the rate that empties the cell over
  !  dtStage, the step of this stage MOSE's Newstate applies (dtlocal*strangcoeff, divided among the sources that
  !  share the cell). rate is per `per`: the cell volume for an absolute rate, 1 for a rate per unit volume; rhoi is
  !  the species density in the cell. Counts what it clamps in nMassLimited.
  subroutine Floor_Guard ( rate, rhoi, per, dtStage )
    implicit none
    real(R8), intent(inout) :: rate(:)
    real(R8), intent(in)    :: rhoi(:), per, dtStage
    integer  :: s
    real(R8) :: avail

    if (dtStage <= 0._R8) return
    do s = 1, size(rate)
      if (rate(s) >= 0._R8) cycle
      avail = rhoi(s) * per / dtStage
      if (-rate(s) > avail) then
        rate(s) = -avail
        !$OMP ATOMIC UPDATE
        nMassLimited = nMassLimited + 1
      endif
    enddo

  end subroutine Floor_Guard


  ! ---------------------------------------------------------------------------------------------
  ! After a condensed solver's solve, outside any parallel region: what it hands the gas is refused before MOSE
  ! applies it when it is not finite (MOSE's conservative update keeps the old temperature of a cell whose energy
  ! residual is NaN, so such a value would act unseen), or when ICE loses mass no gas species receives.
  ! ---------------------------------------------------------------------------------------------

  !> IGLOO's sources: under the source mode the rates IGLOOxMOSE reads (sourceMom, sourceEn, and sourceMass when a
  !  material converts); under the Eulerian mode the fields it turns into sources (velocity, temperature, density,
  !  np of every family).
  subroutine Check_IGLOO_Sources ( MOSE_domain, IGLOO, eulerian )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nsm, nfam
    use hydraMI2_Stoichiometry, only: productsIGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    logical,                intent(in) :: eulerian
    character(len=*), parameter :: SRC = 'source', EUL = 'Eulerian field'
    integer :: b, m, i, j, k

    do b = 1, size(MOSE_domain%blk)
      i = MOSE_domain%blk(b)%dim(1); j = MOSE_domain%blk(b)%dim(2); k = MOSE_domain%blk(b)%dim(3)
      if (eulerian) then
        do m = 1, nfam
          associate ( e => IGLOO%euler(b,m) )
          if (allocated(e%velocity))    call Refuse_4( 'IGLOO', EUL, b, m, 'velocity', 1, e%velocity(:,1:i,1:j,1:k) )
          if (allocated(e%temperature)) call Refuse_3( 'IGLOO', EUL, b, m, 'temperature', e%temperature(1:i,1:j,1:k) )
          if (allocated(e%density))     call Refuse_3( 'IGLOO', EUL, b, m, 'density', e%density(1:i,1:j,1:k) )
          if (allocated(e%np))          call Refuse_3( 'IGLOO', EUL, b, m, 'np', e%np(1:i,1:j,1:k) )
          end associate
        enddo
      else
        associate ( s => IGLOO%source(b) )
        if (allocated(s%sourceMom)) call Refuse_4( 'IGLOO', SRC, b, 0, 'sourceMom', 1, s%sourceMom(:,1:i,1:j,1:k) )
        if (allocated(s%sourceEn))  call Refuse_3( 'IGLOO', SRC, b, 0, 'sourceEn', s%sourceEn(1:i,1:j,1:k) )
        if (productsIGLOO%active .and. allocated(s%sourceMass)) &
          call Refuse_4( 'IGLOO', SRC, b, 0, 'sourceMass', 1, s%sourceMass(1:nsm(),1:i,1:j,1:k) )
        end associate
      endif
    enddo

  end subroutine Check_IGLOO_Sources


  !> ICE's sources, the slots ICExMOSE reads: the momentum source(2:4) and the energy source(nbase-1) of every family,
  !  the mass source(1) of a material that converts. And a family whose material has no record in the products file
  !  that loses mass: ICExMOSE would not hand that mass to the gas.
  subroutine Check_ICE_Sources ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use ICE_Global_m,          only: ngroups, nbase, mat_of, ICE_phase_prefix
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(ICE_domain_type),  intent(in) :: ICE_domain
    character(len=16) :: ie
    integer  :: b, g, i, j, k
    real(R8) :: peak

    do g = 1, ngroups
      write(ie,'(A,I0,A)') 'source(', nbase(g)-1, ')'
      peak = 0._R8
      do b = 1, size(MOSE_domain%blk)
        if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%source)) cycle
        i = MOSE_domain%blk(b)%dim(1); j = MOSE_domain%blk(b)%dim(2); k = MOSE_domain%blk(b)%dim(3)
        associate ( s => ICE_domain%blk(b)%cond_phase(g)%source )
        call Refuse_4( 'ICE', 'source', b, g, 'source', 2, s(2:4,1:i,1:j,1:k) )
        call Refuse_3( 'ICE', 'source', b, g, trim(ie)//', the energy', s(nbase(g)-1,1:i,1:j,1:k) )
        if (Converts(mat_of(g))) then
          call Refuse_3( 'ICE', 'source', b, g, 'source(1), the mass', s(1,1:i,1:j,1:k) )
        else
          peak = max(peak, maxval(abs(s(1,:,:,:))))
        endif
        end associate
      enddo
      if (peak > 0._R8) then
        write(*,'(A,I0,A,I0,3A,ES12.4,A)') E_//'ICE family ', g, ' (material ', mat_of(g), &
          ') loses mass and INPUT/', trim(ICE_phase_prefix), 'products.txt has no record for its material (peak of '// &
          'cond_phase%source(1) ', peak, ' kg/(m3 s)): the gas would not receive it'
        error stop 1
      endif
    enddo

  end subroutine Check_ICE_Sources


  !> Stop when a(slot, i, j, k), what solver who hands the gas in block b (family or group g; 0: none), holds a value
  !  that is not finite, naming the first one: name(first + slot - 1), its cell, its value.
  subroutine Refuse_4 ( who, kind, b, g, name, first, a )
    use MOSE_Finite_m, only: Not_Finite
    implicit none
    character(len=*), intent(in) :: who, kind, name
    integer,          intent(in) :: b, g, first
    real(R8),         intent(in) :: a(:,:,:,:)
    character(len=64) :: label
    integer :: at(4)

    if (.not. any(Not_Finite(a))) return
    at = findloc(Not_Finite(a), .true.)
    write(label,'(2A,I0,A)') name, '(', first + at(1) - 1, ')'
    call Refuse_At( who, kind, b, g, trim(label), at(2:4), a(at(1),at(2),at(3),at(4)) )

  end subroutine Refuse_4


  !> As Refuse_4, for a field of one value per cell, named label.
  subroutine Refuse_3 ( who, kind, b, g, label, a )
    use MOSE_Finite_m, only: Not_Finite
    implicit none
    character(len=*), intent(in) :: who, kind, label
    integer,          intent(in) :: b, g
    real(R8),         intent(in) :: a(:,:,:)
    integer :: at(3)

    if (.not. any(Not_Finite(a))) return
    at = findloc(Not_Finite(a), .true.)
    call Refuse_At( who, kind, b, g, label, at, a(at(1),at(2),at(3)) )

  end subroutine Refuse_3


  subroutine Refuse_At ( who, kind, b, g, label, ijk, v )
    implicit none
    character(len=*), intent(in) :: who, kind, label
    integer,          intent(in) :: b, g, ijk(3)
    real(R8),         intent(in) :: v
    character(len=16) :: val, gw
    character(len=48) :: cell

    write(val,'(ES12.5)') v
    gw = ''
    if (g > 0) write(gw,'(A,I0)') trim(merge(', group  ', ', family ', who == 'ICE'))//' ', g
    write(cell,'(A,3(1X,I0))') ', cell', ijk
    write(*,'(A,I0,A)') E_//who//' handed a non-finite '//kind//' (block ', b, trim(gw)//trim(cell)//', '//label// &
      ') = '//trim(adjustl(val))
    error stop 1

  end subroutine Refuse_At


  ! ---------------------------------------------------------------------------------------------
  ! What the sources in force hand the gas, summed over the interior cells (for the books and the report). Outside
  ! any parallel region.
  ! ---------------------------------------------------------------------------------------------

  !> IGLOO's source-mode rates, before under-relaxation and the floor guard [per s].
  function IGLOO_Handed ( MOSE_domain, IGLOO ) result( h )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nu, nw, np, nprim
    use IGLOO_module,           only: obj_IGLOO
    use hydraMI2_Stoichiometry, only: productsIGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    type(handed_t) :: h
    real(R8), allocatable :: w(:,:,:,:)
    integer  :: b, i, j, k
    real(R8) :: r(nprim)

    if (productsIGLOO%active) allocate(h%mass(nsc), source=0._R8)
    do b = 1, size(MOSE_domain%blk)
      if (.not. allocated(IGLOO%source(b)%sourceMom)) cycle
      call Cells( MOSE_domain%blk(b)%dim, w )
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        r = 0._R8
        r(nu:nw) = IGLOO%source(b)%sourceMom(:,i,j,k)
        r(np)    = IGLOO%source(b)%sourceEn(i,j,k)
        if (productsIGLOO%active) call IGLOO_Mass_To_Gas( IGLOO%source(b)%sourceMass(:,i,j,k), r, w(0,i,j,k) )
        w(1:np,i,j,k) = r(1:np)
      enddo; enddo; enddo
      call Add_Cells( h, w )
    enddo

  end function IGLOO_Handed


  !> ICE's sources of every family, as ICExMOSE applies them without the removals and the floor guard [per s].
  function ICE_Handed ( MOSE_domain, ICE_domain ) result( h )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, np, nprim
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(ICE_domain_type),  intent(in) :: ICE_domain
    type(handed_t) :: h
    real(R8), allocatable :: w(:,:,:,:)
    integer  :: b, g, i, j, k
    real(R8) :: r(nprim), datum, vol

    if (productsICE%active) allocate(h%mass(nsc), source=0._R8)
    do b = 1, size(MOSE_domain%blk)
      call Cells( MOSE_domain%blk(b)%dim, w )
      do g = 1, ngroups
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          r = 0._R8
          call ICE_Cell( mat_of(g), Converts(mat_of(g)), ICE_domain%blk(b)%cond_phase(g)%source(:,i,j,k), nbase(g)-1, &
                         ICE_domain%blk(b)%cond_phase(g)%birth(:,i,j,k), r, datum )
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          w(:,i,j,k) = w(:,i,j,k) + [datum, r(1:np)]*vol
        enddo; enddo; enddo
      enddo
      call Add_Cells( h, w )
    enddo

  end function ICE_Handed


  !> What the removals' stores times share (1: what they still hold), or (rate) the rates in force, hand the gas:
  !  Removal_Cell times the cell volume.
  function Removals_Handed ( MOSE_domain, share, rate ) result( h )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, np, nprim
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    real(R8),               intent(in) :: share
    logical, optional,      intent(in) :: rate
    type(handed_t) :: h
    real(R8), allocatable :: w(:,:,:,:)
    integer  :: b, g, i, j, k
    real(R8) :: r(nprim), datum, vol, gr(5)
    logical  :: asRate

    if (.not. productsICE%active) return
    asRate = .false.
    if (present(rate)) asRate = rate
    allocate(h%mass(nsc), source=0._R8)
    do b = 1, size(MOSE_domain%blk)
      call Cells( MOSE_domain%blk(b)%dim, w )
      do g = 1, ngroups
        if (.not. Converts(mat_of(g))) cycle
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          if (asRate) then
            gr = removed(g)%blk(b)%r(:,i,j,k)*vol
          else
            gr = removed(g)%blk(b)%a(:,i,j,k)*share*vol
          endif
          r = 0._R8; datum = 0._R8
          call Removal_Cell( mat_of(g), gr, r, datum )
          w(:,i,j,k) = w(:,i,j,k) + [datum, r(1:np)]
        enddo; enddo; enddo
      enddo
      call Add_Cells( h, w )
    enddo

  end function Removals_Handed


  !> w(0:np, cell): what a block's cells hand the gas (0 the datum part, 1:np the residual's layout), allocated to the
  !  block's interior and zero.
  subroutine Cells ( dim, w )
    use MOSE_Global_m, only: np
    implicit none
    integer,               intent(in)    :: dim(:)
    real(R8), allocatable, intent(inout) :: w(:,:,:,:)

    if (allocated(w)) deallocate(w)
    allocate(w(0:np, dim(1), dim(2), dim(3)), source=0._R8)

  end subroutine Cells


  !> h plus the sums of w over the cells (each with SUM, whose accuracy the closed-domain balances rely on).
  subroutine Add_Cells ( h, w )
    use MOSE_Global_m, only: nsc, nu, np
    implicit none
    type(handed_t), intent(inout) :: h
    real(R8),       intent(in)    :: w(0:,:,:,:)
    integer :: s

    if (allocated(h%mass)) then
      do s = 1, nsc
        h%mass(s) = h%mass(s) + sum(w(s,:,:,:))
      enddo
    endif
    do s = 1, 3
      h%mom(s) = h%mom(s) + sum(w(nu+s-1,:,:,:))
    enddo
    h%energy = h%energy + sum(w(np,:,:,:))
    h%datum  = h%datum  + sum(w(0,:,:,:))

  end subroutine Add_Cells


  ! ---------------------------------------------------------------------------------------------
  ! ICE's removals, around a MOSE block. Outside any parallel region.
  ! ---------------------------------------------------------------------------------------------

  !> At a MOSE block's start: how the removals' channels reach the gas over its nIter iterations (ICExMOSE reads
  !  them). Rates: where ICE stepped since the last block, the store over the time it was made in becomes the rate
  !  in force; rm holds the rates' totals. Amounts, once ICE has produced sources (ready): each iteration hands
  !  1/nIter of the store.
  subroutine Open_Removals ( MOSE_domain, ready, nIter, rm )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Channels,      only: Store_To_Rate, Open_Share
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    logical,                intent(in)  :: ready
    integer,                intent(in)  :: nIter
    type(removals_t),       intent(out) :: rm
    integer :: g

    if (.not. handRates) then
      if (ready) then
        do g = 1, ngroups
          call Open_Share( removed(g), nIter )
        enddo
      endif
      return
    endif

    rm%asRate = .true.
    if (.not. productsICE%active) return
    do g = 1, ngroups
      if (Converts(mat_of(g))) call Store_To_Rate( removed(g) )
    enddo
    rm%h = Removals_Handed( MOSE_domain, 1._R8, rate=.true. )

  end subroutine Open_Removals


  !> After a MOSE block that made nRun of the iterations Open_Removals shared the stores over (amounts): what they
  !  handed (rm%h), taken out of the stores; the shares closed.
  subroutine Close_Removals ( MOSE_domain, nRun, rm )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Global_m,          only: ngroups, mat_of
    use hydraMI2_Channels,     only: Close_Share, Keep
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    integer,                intent(in)    :: nRun
    type(removals_t),       intent(inout) :: rm
    integer  :: g
    real(R8) :: frac

    ! every family's share was opened alike
    frac = 0._R8
    do g = 1, ngroups
      frac = Close_Share( removed(g), nRun )
    enddo
    if (.not. frac > 0._R8) return
    rm%shared = .true.
    rm%h      = Removals_Handed( MOSE_domain, frac )
    do g = 1, ngroups
      if (Converts(mat_of(g))) call Keep( removed(g), 1._R8 - frac )
    enddo

  end subroutine Close_Removals


  ! ---------------------------------------------------------------------------------------------
  ! ICE's channels, after each of its steps and after its block. Outside any parallel region.
  ! ---------------------------------------------------------------------------------------------

  !> After an ICE iteration: what each family's burn-out removed (cond_phase%gone) moved into its channel, and under
  !  rates each converting family's dt added, cell by cell, to the time its store is being made over. A
  !  time-accurate ICE: each family's dt*source, dt*heat, dt*birth and dt added to the block's mean channel.
  subroutine After_ICE_Step ( ICE_domain )
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use ICE_Config_Types_m,     only: ICE_time_scheme => obj_time_scheme
    use hydraMI2_Channels,      only: Open_Store, Clock
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(ICE_domain_type), intent(inout) :: ICE_domain
    integer :: b, g, i, j, k, ns, dims(3, ICE_domain%nb)

    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        removed(g)%blk(b)%a = removed(g)%blk(b)%a + ICE_domain%blk(b)%cond_phase(g)%gone
        ICE_domain%blk(b)%cond_phase(g)%gone = 0._R8
        if (handRates .and. productsICE%active .and. Converts(mat_of(g))) &
          call Clock( removed(g), b, ICE_domain%blk(b)%cond_phase(g)%dt )
      enddo
    enddo

    if (.not. ICE_time_scheme%time_accurate) return
    if (.not. allocated(iceMean)) then
      allocate(iceMean(ngroups))
      do b = 1, ICE_domain%nb
        dims(:,b) = ICE_domain%blk(b)%dim(1:3)
      enddo
      do g = 1, ngroups
        iceMean(g)%name = 'ice-mean'
        iceMean(g)%from = g
        call Open_Store( iceMean(g), dims, Mean_Size(g), timed=.true., rated=.false. )
      enddo
    endif
    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        ns = size(ICE_domain%blk(b)%cond_phase(g)%source, 1)
        do k = 1, ICE_domain%blk(b)%dim(3)
        do j = 1, ICE_domain%blk(b)%dim(2)
        do i = 1, ICE_domain%blk(b)%dim(1)
          associate (a => iceMean(g)%blk(b)%a, dt => ICE_domain%blk(b)%cond_phase(g)%dt(i,j,k))
            a(1:ns,i,j,k)       = a(1:ns,i,j,k)       + dt*ICE_domain%blk(b)%cond_phase(g)%source(:,i,j,k)
            a(ns+1,i,j,k)       = a(ns+1,i,j,k)       + dt*ICE_domain%blk(b)%cond_phase(g)%heat(i,j,k)
            a(ns+2:ns+6,i,j,k)  = a(ns+2:ns+6,i,j,k)  + dt*ICE_domain%blk(b)%cond_phase(g)%birth(:,i,j,k)
          end associate
        enddo; enddo; enddo
        call Clock( iceMean(g), b, ICE_domain%blk(b)%cond_phase(g)%dt )
      enddo
    enddo
    meanSteps = meanSteps + 1

  contains

    !> The components of family g's mean channel: its source, its heat, its birth (5).
    integer function Mean_Size ( g )
      integer, intent(in) :: g
      Mean_Size = size(ICE_domain%blk(1)%cond_phase(g)%source, 1) + 6
    end function Mean_Size

  end subroutine After_ICE_Step


  !> After a time-accurate ICE block: where it made more than one step, each family's source, heat and birth become
  !  the block's rates, its mean channel's store over its time; the channels emptied.
  subroutine ICE_Block_Mean ( ICE_domain )
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: ngroups
    use hydraMI2_Channels,    only: Empty
    implicit none
    type(ICE_domain_type), intent(inout) :: ICE_domain
    integer  :: b, g, i, j, k, ns
    real(R8) :: t

    if (.not. allocated(iceMean)) return
    do g = 1, ngroups
      if (meanSteps > 1) then
        do b = 1, ICE_domain%nb
          ns = size(ICE_domain%blk(b)%cond_phase(g)%source, 1)
          do k = 1, ICE_domain%blk(b)%dim(3)
          do j = 1, ICE_domain%blk(b)%dim(2)
          do i = 1, ICE_domain%blk(b)%dim(1)
            t = iceMean(g)%blk(b)%t(i,j,k)
            if (.not. t > 0._R8) cycle
            ICE_domain%blk(b)%cond_phase(g)%source(:,i,j,k) = iceMean(g)%blk(b)%a(1:ns,i,j,k)/t
            ICE_domain%blk(b)%cond_phase(g)%heat(i,j,k)     = iceMean(g)%blk(b)%a(ns+1,i,j,k)/t
            ICE_domain%blk(b)%cond_phase(g)%birth(:,i,j,k)  = iceMean(g)%blk(b)%a(ns+2:ns+6,i,j,k)/t
          enddo; enddo; enddo
        enddo
      endif
      call Empty( iceMean(g) )
    enddo
    meanSteps = 0

  end subroutine ICE_Block_Mean


  ! ---------------------------------------------------------------------------------------------
  ! IGLOO's smoke, from the gas's share to ICE
  ! ---------------------------------------------------------------------------------------------

  !> The smoke channel ch makes in a cell, as a rate [per s, the cell's]: mass, momentum (3), energy on ICE's datum
  !  with the kinetic part, number (each drop at the channel's d); born with the parcels' velocity and temperature
  !  weighted by what they burnt there (IGLOO's columns). col: the cell's column of IGLOO's mass source. nostate: a
  !  smoke born where no parcel's columns gave its state. Thread-safe.
  subroutine Smoke_Packet ( ch, col, rate, nostate )
    use IGLOO_variables,        only: nm, smokeRow
    use ICE_Config_Types_m,     only: ICE_condensed => obj_condensed
    use ICE_Global_m,           only: mat_of
    use ICE_Lib_Properties,     only: mat_rho, mat_e, mat_cp_const
    use hydraMI2_Stoichiometry, only: productsIGLOO, Cap_Rates, T_datum_ref
    implicit none
    type(channel_t),   intent(in)  :: ch
    real(R8),          intent(in)  :: col(:)
    real(R8),          intent(out) :: rate(6)
    logical, optional, intent(out) :: nostate
    real(R8), parameter :: pi = 3.14159265358979323846_R8
    real(R8) :: md, w, u(3), T, es, cr(4)
    integer  :: imat, ms

    rate = 0._R8
    if (present(nostate)) nostate = .false.
    imat = ch%from
    ! the mass that makes it: the material burnt into its combustion record, or its cap oxide handed at burnout
    if (productsIGLOO%capped(imat)) then
      cr = Cap_Rates( imat, col(imat), col(nm+2*imat-1), col(nm+2*imat) )
      md = ch%nu*merge(cr(2), cr(4), ch%trigger == SMOKE_OF_BURN)
    elseif (ch%trigger == SMOKE_OF_BURN) then
      md = ch%nu*col(imat)
    else
      md = 0._R8
    endif
    if (md == 0._R8) return
    w = col(smokeRow(imat,1))
    if (w > 0._R8) then
      u = col(smokeRow(imat,2):smokeRow(imat,4))/w
      T = col(smokeRow(imat,5))/w
    else
      u = 0._R8
      T = T_datum_ref
      if (present(nostate)) nostate = .true.
    endif
    ms = mat_of(ch%to)
    if (ICE_condensed(ms)%cs_varies) then
      es = mat_e(ICE_condensed(ms), T)
    else
      es = mat_cp_const(ICE_condensed(ms))*T
    endif
    rate(1)   = md
    rate(2:4) = md*u
    rate(5)   = md*(es + 0.5_R8*dot_product(u, u))
    rate(6)   = md/(pi/6._R8*mat_rho(ICE_condensed(ms), T)*ch%d**3)

  end subroutine Smoke_Packet


  !> After a MOSE iteration on the level it stepped: each smoke channel's rate times the iteration's dtlocal, added cell
  !  by cell to its store (the gas received the same rates less, over the same step). A cell where a material's smoke
  !  had no parcel's state counts once. Outside any parallel region.
  subroutine IGLOO_Smoke_Accumulate ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    integer  :: b, i, j, k, c, mNoState
    real(R8) :: sr(6), dt
    logical  :: nostate

    if (MOSE_domain%mg_level > 1) then
      call Smoke_Accumulate_Coarse( MOSE_domain )
      return
    endif

    do b = 1, size(MOSE_domain%blk)
      if (.not. allocated(IGLOO%source(b)%sourceMass)) cycle
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        dt = MOSE_domain%blk(b)%dtlocal(i,j,k)
        mNoState = 0
        do c = 1, size(smoke)
          call Smoke_Packet( smoke(c), IGLOO%source(b)%sourceMass(:,i,j,k), sr, nostate )
          ! the channels of a material are consecutive
          if (nostate .and. smoke(c)%from /= mNoState) then
            nSmokeNoState = nSmokeNoState + 1
            mNoState      = smoke(c)%from
          endif
          if (sr(1) == 0._R8) cycle
          smoke(c)%blk(b)%a(:,i,j,k) = smoke(c)%blk(b)%a(:,i,j,k) + sr*dt
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine IGLOO_Smoke_Accumulate


  !> IGLOO_Smoke_Accumulate on a coarse level, from the smoke's rates restricted to it (each sub-cell's, summed); a
  !  coarse cell's amounts go to the store of its first sub-cell (only totals are read in a steady run).
  subroutine Smoke_Accumulate_Coarse ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer  :: b, i, j, k, c, f(3)
    real(R8) :: dt

    do b = 1, size(MOSE_domain%blk)
      if (.not. allocated(coarse%smk(b)%r)) cycle
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        f  = ([i, j, k] - 1)*coarse%ratio(:,b) + 1
        dt = MOSE_domain%blk(b)%dtlocal(i,j,k)
        nSmokeNoState = nSmokeNoState + coarse%smk(b)%nostate(i,j,k)
        do c = 1, size(smoke)
          if (coarse%smk(b)%r(1,c,i,j,k) == 0._R8) cycle
          smoke(c)%blk(b)%a(:,f(1),f(2),f(3)) = smoke(c)%blk(b)%a(:,f(1),f(2),f(3)) + coarse%smk(b)%r(:,c,i,j,k)*dt
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine Smoke_Accumulate_Coarse


  !> After a MOSE block: the smoke its iterations made (IGLOO_Smoke_Accumulate), which the gas did not get, handed to
  !  each ICE family a smoke channel feeds, in the channels' form (every smoke channel steps on MOSE's clock into
  !  ICE's). Amounts: added to ICE's state cell by cell (per volume). Rates: the rates of IGLOO's sources as they stand
  !  (srcOn; else 0), per volume, set as the family's source (cond_phase%ext). sm: the totals, energies on the gas's
  !  formation scale; the stores emptied. Outside any parallel region.
  subroutine IGLOO_Smoke_Deliver ( MOSE_domain, ICE_domain, IGLOO, srcOn, sm )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    use ICE_Global_m,          only: ngroups, mat_of, nbase, ncond
    use ICE_Config_Types_m,    only: ICE_condensed => obj_condensed
    use ICE_Lib_Model,         only: prim_2_cons, cons_2_prim, assign_prim_2_cons, assign_cons_2_prim
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: srcOn
    type(smoke_t),          intent(out)   :: sm
    integer  :: b, g, i, j, k, ns, c
    real(R8) :: vol, ac(6), sr(6), pk(6), dhg
    real(R8), allocatable :: q(:)
    logical  :: asRate

    sm%nNoState   = nSmokeNoState
    nSmokeNoState = 0
    asRate = smoke(1)%asRate
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        if (.not. any(smoke(:)%to == g)) cycle
        ! the smoke family's datum offset: its material's, whichever channel feeds it
        do c = 1, size(smoke)
          if (smoke(c)%to == g) dhg = smoke(c)%dh
        enddo
        ! the conversions of family g's closure (ICE points them at the family it last stepped)
        call assign_prim_2_cons(g)
        call assign_cons_2_prim(g)
        ns = nbase(g)
        if (allocated(q)) deallocate(q)
        allocate(q(ncond(g)))
        if (asRate) then
          if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%ext)) &
            allocate(ICE_domain%blk(b)%cond_phase(g)%ext(ncond(g), ICE_domain%blk(b)%dim(1), &
                     ICE_domain%blk(b)%dim(2), ICE_domain%blk(b)%dim(3)))
          ICE_domain%blk(b)%cond_phase(g)%ext = 0._R8
        endif
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          ! rates: the source of family g from IGLOO's sources as they stand
          if (asRate .and. srcOn) then
            if (allocated(IGLOO%source(b)%sourceMass)) then
              sr = 0._R8
              do c = 1, size(smoke)
                if (smoke(c)%to /= g) cycle
                call Smoke_Packet( smoke(c), IGLOO%source(b)%sourceMass(:,i,j,k), pk )
                sr = sr + pk
              enddo
              q = 0._R8
              q(1:4)  = sr(1:4)/vol
              q(ns-1) = sr(5)/vol
              q(ns)   = sr(6)/vol
              ICE_domain%blk(b)%cond_phase(g)%ext(:,i,j,k) = q
              sm%rate(1:4) = sm%rate(1:4) + sr(1:4)
              sm%rate(5)   = sm%rate(5) + sr(5) + dhg*sr(1)
              sm%rate(6)   = sm%rate(6) + sr(6)
            endif
          endif
          ! the amounts of the block: out of the gas's share, and into ICE's state when they go as amounts
          ac = 0._R8
          do c = 1, size(smoke)
            if (smoke(c)%to == g) ac = ac + smoke(c)%blk(b)%a(:,i,j,k)
          enddo
          if (all(ac == 0._R8)) cycle
          if (.not. asRate) then
            q = prim_2_cons(ICE_domain%blk(b)%cond_phase(g)%prim(:,i,j,k), ICE_condensed(mat_of(g)))
            q(1:4)  = q(1:4)  + ac(1:4)/vol
            q(ns-1) = q(ns-1) + ac(5)/vol
            q(ns)   = q(ns)   + ac(6)/vol
            ICE_domain%blk(b)%cond_phase(g)%prim(:,i,j,k) = cons_2_prim(q, ICE_condensed(mat_of(g)))
          endif
          sm%amount(1:4) = sm%amount(1:4) + ac(1:4)
          sm%amount(5)   = sm%amount(5) + ac(5) + dhg*ac(1)
          do c = 1, size(smoke)
            if (smoke(c)%to == g) smoke(c)%blk(b)%a(:,i,j,k) = 0._R8
          enddo
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine IGLOO_Smoke_Deliver


  ! ---------------------------------------------------------------------------------------------
  ! MOSE's coarse levels: the condensed sources restricted from the finest grid to the level MOSE steps. Outside any
  ! parallel region.
  ! ---------------------------------------------------------------------------------------------

  !> Before a MOSE iteration: on a coarser level than the finest grid, the condensed sources in force (IGLOO_src,
  !  ICE_src) restricted to it when that level or those sources changed since the last restriction (rebuilt); the
  !  finest grid then lags the level until Fill_Finest_Grid.
  subroutine Before_MOSE_Iteration ( MOSE_simulation, IGLOO, ICE_simulation, IGLOO_src, ICE_src, rebuilt )
    use MOSE_Advanced_Types_m, only: MOSE_simulation_type
    use MOSE_Config_Types_m,   only: MOSE_multigrid => obj_multigrid
    use ICE_Advanced_Types_m,  only: ICE_simulation_type
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_simulation_type), intent(in)  :: MOSE_simulation
    type(obj_IGLOO),            intent(in)  :: IGLOO
    type(ICE_simulation_type),  intent(in)  :: ICE_simulation
    type(source_state_t),       intent(in)  :: IGLOO_src, ICE_src
    logical,                    intent(out) :: rebuilt
    integer :: lev

    lev       = MOSE_multigrid%MG_level
    rebuilt   = .false.
    fineStale = lev > 1
    ! an iteration on the finest grid evaluates its gas
    if ( lev == 1 ) gasEvaluated = .true.
    if ( lev == 1 ) return
    if ( lev == coarse%held .and. Same_State(IGLOO_src, coarse%iglooAt) .and. Same_State(ICE_src, coarse%iceAt) ) return

    nTotals = 0
    if ( IGLOO_src%ready ) call Restrict_IGLOO_Sources( MOSE_simulation%domain(1), MOSE_simulation%domain(lev), IGLOO )
    if ( ICE_src%ready ) call Restrict_ICE_Sources( MOSE_simulation%domain(1), MOSE_simulation%domain(lev), &
                                                    ICE_simulation%domain(1) )
    coarse%held    = lev
    coarse%iglooAt = IGLOO_src
    coarse%iceAt   = ICE_src
    rebuilt        = .true.

  contains

    !> The same sources in force: produced by the same solve.
    logical function Same_State ( a, b )
      type(source_state_t), intent(in) :: a, b
      Same_State = (a%ready .eqv. b%ready) .and. (a%eulerian .eqv. b%eulerian) .and. a%cycle == b%cycle .and. &
                   a%step == b%step
    end function Same_State

  end subroutine Before_MOSE_Iteration


  !> The totals of the last restriction (Before_MOSE_Iteration), for the report.
  function Restricted_Totals () result( t )
    implicit none
    type(restricted_t), allocatable :: t(:)

    allocate(t(nTotals))
    if (nTotals > 0) t = totals(1:nTotals)

  end function Restricted_Totals


  !> One total of a restriction recorded: who, what, block b, the finest grid's f, the coarse level's c, magnitude a.
  subroutine Add_Total ( who, what, b, f, c, a )
    implicit none
    character(len=*), intent(in) :: who, what
    integer,          intent(in) :: b
    real(R8),         intent(in) :: f, c, a
    type(restricted_t), allocatable :: more(:)

    if (.not. allocated(totals)) allocate(totals(64))
    if (nTotals == size(totals)) then
      allocate(more(2*size(totals)))
      more(1:nTotals) = totals
      call move_alloc(more, totals)
    endif
    nTotals = nTotals + 1
    totals(nTotals) = restricted_t( who=who, what=what, b=b, onFine=f, onCoarse=c, magnitude=a )

  end subroutine Add_Total


  !> The coarse level crs set up: per block the sub-cells per coarse cell (finest-grid cells over coarse cells) and
  !  their volume. Stops on a level whose cells do not hold a whole number of the finest grid's.
  subroutine Coarse_Level_Setup ( fine, crs )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    type(MOSE_domain_type), intent(in) :: fine, crs
    integer :: b, n(3), r(3), lo(3), hi(3), ci, cj, ck

    if (coarse%level == crs%mg_level) return
    if (allocated(coarse%ratio)) deallocate(coarse%ratio)
    if (allocated(coarse%vol))   deallocate(coarse%vol)
    allocate(coarse%ratio(3, size(crs%blk)), coarse%vol(size(crs%blk)))
    do b = 1, size(crs%blk)
      n = crs%blk(b)%dim
      if (any(mod(fine%blk(b)%dim, n) /= 0)) then
        write(*,'(A,I0,A,I0,A,3(1X,I0),A,3(1X,I0),A)') E_//'block ', b, ' of MOSE''s grid level ', crs%mg_level, &
          ' has', n, ' cells against', fine%blk(b)%dim, ' on the finest grid: not a whole number of sub-cells per cell'
        error stop 1
      endif
      r = fine%blk(b)%dim / n
      coarse%ratio(:,b) = r
      allocate(coarse%vol(b)%v(n(1), n(2), n(3)))
      do ck = 1, n(3)
      do cj = 1, n(2)
      do ci = 1, n(1)
        lo = ([ci, cj, ck] - 1)*r + 1
        hi = [ci, cj, ck]*r
        coarse%vol(b)%v(ci,cj,ck) = sum(fine%blk(b)%vol(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
      enddo; enddo; enddo
    enddo
    coarse%level = crs%mg_level

  end subroutine Coarse_Level_Setup


  !> Stop when an array of a solver does not hold the cells of block b of MOSE's finest grid (fine): ext, its extents.
  subroutine Check_Cells ( what, ext, fine, b )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    character(len=*),       intent(in) :: what
    integer,                intent(in) :: ext(3), b
    type(MOSE_domain_type), intent(in) :: fine

    if (all(ext == fine%blk(b)%dim)) return
    write(*,'(A,I0,A,3(1X,I0),A,3(1X,I0))') E_//what//' of block ', b, ' has', ext, ' cells, MOSE''s finest grid', &
      fine%blk(b)%dim
    error stop 1

  end subroutine Check_Cells


  !> IGLOO's sources restricted to the coarse level crs: rates summed over the sub-cells (the smoke's too, each
  !  sub-cell's: they are not linear in IGLOO's columns), Eulerian fields weighted as IGLOO's normalizeEUL defines them
  !  (density, number by volume; velocity, temperature by mass).
  subroutine Restrict_IGLOO_Sources ( fine, crs, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    use IGLOO_variables,       only: toll
    implicit none
    type(MOSE_domain_type), intent(in) :: fine, crs
    type(obj_IGLOO),        intent(in) :: IGLOO
    character(len=6), parameter :: EUL_NAME(6) = [character(len=6) :: 'mass', 'number', 'rhou', 'rhov', 'rhow', 'rhoT']
    integer  :: b, m, c, n(3), r(3), i, j, k, ci, cj, ck, mNoState
    real(R8) :: v, sv, sm, sn, smu(3), smT, ft(6), ct(6), at(6), sr(6)
    logical  :: nostate
    character(len=32) :: nm

    call Coarse_Level_Setup( fine, crs )
    if (allocated(coarse%src)) deallocate(coarse%src)
    if (allocated(coarse%smk)) deallocate(coarse%smk)
    if (allocated(coarse%eul)) deallocate(coarse%eul)
    allocate(coarse%src(size(crs%blk)), coarse%smk(size(crs%blk)))
    if (allocated(IGLOO%euler)) allocate(coarse%eul(size(crs%blk), size(IGLOO%euler, 2)))
    do b = 1, size(crs%blk)
      n = crs%blk(b)%dim
      r = coarse%ratio(:,b)
      ! the source mode's rates: sums
      if (allocated(IGLOO%source(b)%sourceMass)) then
        call Check_Cells( 'IGLOO''s sourceMass', shape3(IGLOO%source(b)%sourceMass), fine, b )
        allocate(coarse%src(b)%sourceMass(size(IGLOO%source(b)%sourceMass, 1), n(1), n(2), n(3)))
        call Sum_Sub_Cells( IGLOO%source(b)%sourceMass, coarse%src(b)%sourceMass )
        do c = 1, size(IGLOO%source(b)%sourceMass, 1)
          write(nm,'(A,I0,A)') 'igloo-sourceMass(', c, ')'
          call Add_Total( 'IGLOO', nm, b, sum(IGLOO%source(b)%sourceMass(c,:,:,:)), &
                          sum(coarse%src(b)%sourceMass(c,:,:,:)), sum(abs(IGLOO%source(b)%sourceMass(c,:,:,:))) )
        enddo
        ! the smoke's rates: each sub-cell's, summed; a sub-cell without a parcel's state counted once per material
        if (smokeIGOn) then
          allocate(coarse%smk(b)%r(6, size(smoke), n(1), n(2), n(3)), source=0._R8)
          allocate(coarse%smk(b)%nostate(n(1), n(2), n(3)), source=0)
          do ck = 1, n(3)
          do cj = 1, n(2)
          do ci = 1, n(1)
            do k = (ck-1)*r(3)+1, ck*r(3)
            do j = (cj-1)*r(2)+1, cj*r(2)
            do i = (ci-1)*r(1)+1, ci*r(1)
              mNoState = 0
              do c = 1, size(smoke)
                call Smoke_Packet( smoke(c), IGLOO%source(b)%sourceMass(:,i,j,k), sr, nostate )
                coarse%smk(b)%r(:,c,ci,cj,ck) = coarse%smk(b)%r(:,c,ci,cj,ck) + sr
                if (nostate .and. smoke(c)%from /= mNoState) then
                  coarse%smk(b)%nostate(ci,cj,ck) = coarse%smk(b)%nostate(ci,cj,ck) + 1
                  mNoState = smoke(c)%from
                endif
              enddo
            enddo; enddo; enddo
          enddo; enddo; enddo
        endif
      endif
      if (allocated(IGLOO%source(b)%sourceMom)) then
        call Check_Cells( 'IGLOO''s sourceMom', shape3(IGLOO%source(b)%sourceMom), fine, b )
        allocate(coarse%src(b)%sourceMom(3, n(1), n(2), n(3)))
        call Sum_Sub_Cells( IGLOO%source(b)%sourceMom, coarse%src(b)%sourceMom )
        do c = 1, 3
          write(nm,'(A,I0,A)') 'igloo-sourceMom(', c, ')'
          call Add_Total( 'IGLOO', nm, b, sum(IGLOO%source(b)%sourceMom(c,:,:,:)), &
                          sum(coarse%src(b)%sourceMom(c,:,:,:)), sum(abs(IGLOO%source(b)%sourceMom(c,:,:,:))) )
        enddo
      endif
      if (allocated(IGLOO%source(b)%sourceEn)) then
        call Check_Cells( 'IGLOO''s sourceEn', shape(IGLOO%source(b)%sourceEn), fine, b )
        allocate(coarse%src(b)%sourceEn(n(1), n(2), n(3)))
        do ck = 1, n(3)
        do cj = 1, n(2)
        do ci = 1, n(1)
          coarse%src(b)%sourceEn(ci,cj,ck) = sum(IGLOO%source(b)%sourceEn((ci-1)*r(1)+1:ci*r(1), &
                                                 (cj-1)*r(2)+1:cj*r(2), (ck-1)*r(3)+1:ck*r(3)))
        enddo; enddo; enddo
        call Add_Total( 'IGLOO', 'igloo-sourceEn', b, sum(IGLOO%source(b)%sourceEn), sum(coarse%src(b)%sourceEn), &
                        sum(abs(IGLOO%source(b)%sourceEn)) )
      endif
      ! the Eulerian fields: density and number by volume, velocity and temperature by mass
      if (.not. allocated(coarse%eul)) cycle
      do m = 1, size(coarse%eul, 2)
        if (.not. allocated(IGLOO%euler(b,m)%density)) cycle
        associate (fe => IGLOO%euler(b,m))
        call Check_Cells( 'IGLOO''s Eulerian density', shape(fe%density), fine, b )
        allocate(coarse%eul(b,m)%density(n(1), n(2), n(3)), coarse%eul(b,m)%np(n(1), n(2), n(3)), &
                 coarse%eul(b,m)%velocity(3, n(1), n(2), n(3)), coarse%eul(b,m)%temperature(n(1), n(2), n(3)))
        ! the totals of the finest grid, cell by cell
        ft = 0._R8; ct = 0._R8; at = 0._R8
        do k = 1, fine%blk(b)%dim(3)
        do j = 1, fine%blk(b)%dim(2)
        do i = 1, fine%blk(b)%dim(1)
          v  = fine%blk(b)%vol(i,j,k)
          ft = ft + [fe%density(i,j,k), fe%np(i,j,k), fe%density(i,j,k)*fe%velocity(:,i,j,k), &
                     fe%density(i,j,k)*fe%temperature(i,j,k)]*v
          at = at + abs([fe%density(i,j,k), fe%np(i,j,k), fe%density(i,j,k)*fe%velocity(:,i,j,k), &
                         fe%density(i,j,k)*fe%temperature(i,j,k)])*v
        enddo; enddo; enddo
        do ck = 1, n(3)
        do cj = 1, n(2)
        do ci = 1, n(1)
          sm = 0._R8; sn = 0._R8; smu = 0._R8; smT = 0._R8
          do k = (ck-1)*r(3)+1, ck*r(3)
          do j = (cj-1)*r(2)+1, cj*r(2)
          do i = (ci-1)*r(1)+1, ci*r(1)
            v   = fine%blk(b)%vol(i,j,k)
            sm  = sm  + fe%density(i,j,k)*v
            sn  = sn  + fe%np(i,j,k)*v
            smu = smu + fe%density(i,j,k)*fe%velocity(:,i,j,k)*v
            smT = smT + fe%density(i,j,k)*fe%temperature(i,j,k)*v
          enddo; enddo; enddo
          sv = coarse%vol(b)%v(ci,cj,ck)
          coarse%eul(b,m)%density(ci,cj,ck) = sm/sv
          coarse%eul(b,m)%np(ci,cj,ck)      = sn/sv
          if (coarse%eul(b,m)%density(ci,cj,ck) > toll) then
            coarse%eul(b,m)%velocity(:,ci,cj,ck)  = smu/sm
            coarse%eul(b,m)%temperature(ci,cj,ck) = smT/sm
          else
            coarse%eul(b,m)%velocity(:,ci,cj,ck)  = 0._R8
            coarse%eul(b,m)%temperature(ci,cj,ck) = 0._R8
          endif
          ct = ct + [coarse%eul(b,m)%density(ci,cj,ck), coarse%eul(b,m)%np(ci,cj,ck), &
                     coarse%eul(b,m)%density(ci,cj,ck)*coarse%eul(b,m)%velocity(:,ci,cj,ck), &
                     coarse%eul(b,m)%density(ci,cj,ck)*coarse%eul(b,m)%temperature(ci,cj,ck)]*sv
        enddo; enddo; enddo
        end associate
        do c = 1, 6
          write(nm,'(A,I0,2A)') 'igloo-euler(', m, ')-', trim(EUL_NAME(c))
          call Add_Total( 'IGLOO', nm, b, ft(c), ct(c), at(c) )
        enddo
      enddo
    enddo

  contains

    !> y, the sum over the sub-cells of each coarse cell of block b of x (component first)
    subroutine Sum_Sub_Cells ( x, y )
      real(R8), intent(in)  :: x(:,:,:,:)
      real(R8), intent(out) :: y(:,:,:,:)
      y = 0._R8
      do ck = 1, n(3)
      do cj = 1, n(2)
      do ci = 1, n(1)
        do k = (ck-1)*r(3)+1, ck*r(3)
        do j = (cj-1)*r(2)+1, cj*r(2)
        do i = (ci-1)*r(1)+1, ci*r(1)
          y(:,ci,cj,ck) = y(:,ci,cj,ck) + x(:,i,j,k)
        enddo; enddo; enddo
      enddo; enddo; enddo
    end subroutine Sum_Sub_Cells

  end subroutine Restrict_IGLOO_Sources


  !> ICE's per-volume rates restricted to the coarse level crs, and the removals' rates in force: per coarse cell the
  !  sub-cells' rate times volume summed, over the coarse cell's own volume, which ICE_Residual multiplies back (on a
  !  curved mesh it is not the sub-cells' sum). A start on the coarse grids steps MOSE in pseudo-time, so the removals'
  !  channels go as rates (handRates): their rates alone are restricted.
  subroutine Restrict_ICE_Sources ( fine, crs, ICE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Global_m,          only: ngroups
    implicit none
    type(MOSE_domain_type), intent(in) :: fine, crs
    type(ICE_domain_type),  intent(in) :: ICE_domain
    integer :: b, g, n(3), r(3), i, j, k, ci, cj, ck

    call Coarse_Level_Setup( fine, crs )
    if (allocated(coarse%ice%blk)) deallocate(coarse%ice%blk)
    allocate(coarse%ice%blk(size(crs%blk)))
    if (allocated(coarse%removed)) deallocate(coarse%removed)
    allocate(coarse%removed(size(removed)))
    do g = 1, size(removed)
      coarse%removed(g)%name   = removed(g)%name
      coarse%removed(g)%from   = removed(g)%from
      coarse%removed(g)%asRate = removed(g)%asRate
      allocate(coarse%removed(g)%blk(size(crs%blk)))
    enddo
    do b = 1, size(crs%blk)
      n = crs%blk(b)%dim
      r = coarse%ratio(:,b)
      allocate(coarse%ice%blk(b)%cond_phase(ngroups))
      do g = 1, ngroups
        associate (fc => ICE_domain%blk(b)%cond_phase(g))
        if (allocated(fc%source)) then
          call Check_Cells( 'ICE''s source', shape3(fc%source), fine, b )
          allocate(coarse%ice%blk(b)%cond_phase(g)%source(size(fc%source, 1), n(1), n(2), n(3)))
          call Per_Coarse_Volume( fc%source, coarse%ice%blk(b)%cond_phase(g)%source, 'ice-source' )
        endif
        if (allocated(fc%birth)) then
          call Check_Cells( 'ICE''s birth', shape3(fc%birth), fine, b )
          allocate(coarse%ice%blk(b)%cond_phase(g)%birth(size(fc%birth, 1), n(1), n(2), n(3)))
          call Per_Coarse_Volume( fc%birth, coarse%ice%blk(b)%cond_phase(g)%birth, 'ice-birth' )
        endif
        end associate
        if (allocated(removed(g)%blk(b)%r)) then
          call Check_Cells( 'the rate of ICE''s removals', shape3(removed(g)%blk(b)%r), fine, b )
          allocate(coarse%removed(g)%blk(b)%r(size(removed(g)%blk(b)%r, 1), n(1), n(2), n(3)))
          call Per_Coarse_Volume( removed(g)%blk(b)%r, coarse%removed(g)%blk(b)%r, 'ice-removed-rate' )
        endif
      enddo
    enddo

  contains

    !> y, per coarse cell of block b, the sum over its sub-cells of x times their volume, over the coarse cell's
    !  volume; one total per component
    subroutine Per_Coarse_Volume ( x, y, what )
      real(R8),         intent(in)  :: x(:,:,:,:)
      real(R8),         intent(out) :: y(:,:,:,:)
      character(len=*), intent(in)  :: what
      real(R8) :: ft(size(x,1)), ct(size(x,1)), at(size(x,1))
      integer  :: c
      character(len=32) :: nm
      ! the totals of the finest grid, cell by cell
      ft = 0._R8; ct = 0._R8; at = 0._R8
      do k = 1, fine%blk(b)%dim(3)
      do j = 1, fine%blk(b)%dim(2)
      do i = 1, fine%blk(b)%dim(1)
        ft = ft + x(:,i,j,k)*fine%blk(b)%vol(i,j,k)
        at = at + abs(x(:,i,j,k))*fine%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
      do ck = 1, n(3)
      do cj = 1, n(2)
      do ci = 1, n(1)
        y(:,ci,cj,ck) = 0._R8
        do k = (ck-1)*r(3)+1, ck*r(3)
        do j = (cj-1)*r(2)+1, cj*r(2)
        do i = (ci-1)*r(1)+1, ci*r(1)
          y(:,ci,cj,ck) = y(:,ci,cj,ck) + x(:,i,j,k)*fine%blk(b)%vol(i,j,k)
        enddo; enddo; enddo
        y(:,ci,cj,ck) = y(:,ci,cj,ck)/crs%blk(b)%vol(ci,cj,ck)
        ct = ct + y(:,ci,cj,ck)*crs%blk(b)%vol(ci,cj,ck)
      enddo; enddo; enddo
      do c = 1, size(x,1)
        write(nm,'(A,I0,A,I0,A)') what//'(', g, ',', c, ')'
        call Add_Total( 'ICE', nm, b, ft(c), ct(c), at(c) )
      enddo
    end subroutine Per_Coarse_Volume

  end subroutine Restrict_ICE_Sources


  !> The last three extents of a rank-4 array: its cells.
  pure function shape3 ( x ) result( e )
    implicit none
    real(R8), intent(in) :: x(:,:,:,:)
    integer :: e(3)

    e = [size(x,2), size(x,3), size(x,4)]

  end function shape3

end module hydraMI2_Mod_Coupling
