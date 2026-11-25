!> The exchange between MOSE and the two condensed solvers of hydra-MI2.
!
!  The gas is the ledger of every exchange. A condensed mass that leaves its solver enters the gas's books at the
!  solver's energy plus its datum offset; a record of the products file turns it into gas species; a condensed product
!  of that record born in another solver leaves the gas's books at the receiver's energy plus its own offset. So a
!  route between two solvers needs no code of its own: it is a record, read from the products file ATLAS writes.
!
!  Within a coupling step:
!    MOSExIGLOO, MOSExICE : the gas state a condensed solver steps on
!    IGLOOxMOSE, ICExMOSE : what it hands the gas, into MOSE's residual at every Runge-Kutta stage
!  What one cell hands the gas is written once, in the cell kernels IGLOO_Cell, ICE_Cell and Removal_Cell; the
!  injection applies them, IGLOO_Handed, ICE_Handed and Removals_Handed sum them for the books (hydraMI2_Report).
!
!  What lands outside a solver's step goes through a channel (hydraMI2_Channels):
!    removed(g)  ICE family g -> gas: what its burn-out removed (take_removed after each ICE step)
!    iceMean(g)  a time-accurate ICE block of several steps hands the gas its steps' dt-weighted mean
!    gen(c)      IGLOO material -> ICE family: the mass a record of the material generates, born in the family
!  They go as rates when MOSE or ICE steps in pseudo-time (handRates), as amounts when both are time-accurate.
module hydraMI2_Mod_Coupling
  use iso_fortran_env,   only: R8 => real64
  use hydraMI2_Channels, only: channel_t

  implicit none
  private

  public :: source_state_t, handed_t, removals_t, gen_t
  public :: Setup_ICE_Coupling, Setup_IGLOO_Coupling, Set_IGLOO_Switches
  public :: Gas_Field_For_IGLOO, MOSExIGLOO, MOSExICE, IGLOOxMOSE, ICExMOSE
  public :: Check_IGLOO_Sources, Check_ICE_Sources
  public :: IGLOO_Handed, ICE_Handed, Removals_Handed
  public :: Removed_Store_Setup, Removals_Timed
  public :: Open_Removals, Close_Removals, After_ICE_Step, ICE_Block_Mean
  public :: IGLOO_Gen_Accumulate, IGLOO_Gen_Deliver

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

  !> IGLOO's generated mass a MOSE block made, withheld from the gas and given to ICE (IGLOO_Gen_Deliver): its amounts
  !  (mass, momentum (3), energy on the gas's datum); under handRates the rates set in ICE (and the number); the cells
  !  of the block's iterations where a generated mass was born with no parcel's state.
  type :: gen_t
    real(R8) :: amount(5) = 0._R8, rate(6) = 0._R8
    integer  :: nNoState  = 0
  end type gen_t

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

  !> removed(g): what the burn-out of ICE family g removed and the gas has not received yet, per unit volume: the
  !  gas's mass, momentum (3), energy, then the generated mass (take_removed). Under rates timed: the store over the
  !  time it was made in is the rate in force. Restored by the hand-off file of a restart.
  type(channel_t), allocatable, public :: removed(:)

  !> iceMean(g): over a time-accurate ICE block, family g's dt*source, dt*heat, dt*gen_rate and dt, cell by cell,
  !  made into the block's rates after it (ICE_Block_Mean) when it made meanSteps > 1 steps.
  type(channel_t), allocatable :: iceMean(:)
  integer :: meanSteps = 0

  !> gen(c): the mass a record of an IGLOO material (from) generates, born in an ICE family (to): GEN_OF_PROCESS, the
  !  material's process record; GEN_OF_RELEASE, the release record of the material it carries, handed at burnout.
  !  d: the material's d-gen; dh: the generated material's datum offset in ICE. Its store: what the current MOSE
  !  block made (mass, momentum (3), energy on ICE's datum, number). genIGOn: one at least.
  integer, parameter :: GEN_OF_PROCESS = 1, GEN_OF_RELEASE = 2
  type(channel_t), allocatable :: gen(:)
  logical, protected, public   :: genIGOn = .false.
  integer :: nGenNoState = 0

  character(len=*), parameter :: E_ = ' [ERROR] MI2: '

contains

  ! ---------------------------------------------------------------------------------------------
  ! Setup
  ! ---------------------------------------------------------------------------------------------

  !> ICE's side of the coupler, after ICE%setup: how what lands between steps goes (handRates), the chemistry of
  !  ICE's materials, the mass they generate (declared to ICE), the datum offsets, the removals' channels. ICE credits
  !  a material its energy at T_ref on its own datum: its table's, or cp*T.
  subroutine Setup_ICE_Coupling ( ICE_domain )
    use MOSE_Config_Types_m,    only: MOSE_time_scheme => obj_time_scheme
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ICE_phase_prefix
    use ICE_Lib_Properties,     only: mat_lv
    use ICE_Config_Types_m,     only: ICE_condensed => obj_condensed, ICE_time_scheme => obj_time_scheme
    use ICE_Mod_Exchange,       only: set_generation
    use hydraMI2_Channels,      only: As_Rate
    use hydraMI2_Chemistry,     only: chemICE, IS_ICE, T_datum_ref, Load_Chemistry, Set_Datum_Offsets, &
                                      Set_Product_Datum, Report_Chemistry
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    character(len=64), allocatable :: names(:)
    real(R8),          allocatable :: eref(:), Lv(:)
    logical,           allocatable :: none(:)
    character(len=:),  allocatable :: why
    integer  :: m, nm, gs
    real(R8) :: dh

    handRates = As_Rate( MOSE_time_scheme%time_accurate, ICE_time_scheme%time_accurate )

    nm = size(ICE_condensed)
    allocate(names(nm), eref(nm), Lv(nm), none(nm))
    none = .false.
    do m = 1, nm
      names(m) = ICE_condensed(m)%name
      eref(m)  = ICE_Energy( m, T_datum_ref )
      Lv(m)    = mat_lv( ICE_condensed(m), T_datum_ref )
    enddo
    call Load_Chemistry( IS_ICE, ICE_phase_prefix, names, ICE_condensed(:)%evapSelect > 0, &
                         ICE_condensed(:)%combSelect > 0, none, none )
    call Set_Datum_Offsets( IS_ICE, eref, Lv, eref )

    ! the mass ICE's burning materials generate, which ICE makes itself, in the family the record names
    do m = 1, nm
      if ( .not. chemICE%mat(m)%process%on .or. len_trim(chemICE%mat(m)%process%phase) == 0 ) cycle
      call Gen_Family( trim(names(m)), chemICE%mat(m)%process%product, chemICE%mat(m)%process%phase, &
                       chemICE%mat(m)%process%hProduct, gs, dh )
      call set_generation( m, chemICE%mat(m)%process%product, chemICE%mat(m)%process%nuProduct, why )
      if ( len(why) > 0 ) then
        write(*,'(A)') E_//'ICE material '//trim(names(m))//' cannot generate '// &
          trim(chemICE%mat(m)%process%product)//': '//why
        error stop 1
      endif
      call Set_Product_Datum( IS_ICE, m, .false., dh )
      write(*,'(A,I0,A,ES15.7,A,ES12.4,A)') ' [MI2] ICE material '//trim(names(m))//' generates '// &
        trim(chemICE%mat(m)%process%product)//' (family ', gs, '): ', chemICE%mat(m)%process%nuProduct, &
        ' kg per kg burnt, d-gen ', ICE_condensed(m)%dGen, ' m'
    enddo
    call Report_Chemistry( IS_ICE )

    ! the removals reach the gas as the rate of their store, timed by ICE's steps, or as amounts
    call Removed_Store_Setup( ICE_domain, timed = handRates .and. chemICE%active )
    removed(:)%asRate = handRates

  end subroutine Setup_ICE_Coupling


  !> ICE's energy of material m at T on its own datum: its table's, or cp*T.
  real(R8) function ICE_Energy ( m, T )
    use ICE_Config_Types_m, only: ICE_condensed => obj_condensed
    use ICE_Lib_Properties, only: mat_e, mat_cp_const
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T

    if ( ICE_condensed(m)%cs_varies ) then
      ICE_Energy = mat_e( ICE_condensed(m), T )
    else
      ICE_Energy = mat_cp_const( ICE_condensed(m) ) * T
    endif

  end function ICE_Energy


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
      call Open_Store( removed(g), dims, 6, timed=timed, rated=timed )
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


  !> IGLOO's side of the coupler, after IGLOO%setup: the family -> material map, the chemistry of IGLOO's materials
  !  (what each carries, and the carried material it hands the gas at burnout), the mass they generate in ICE, the
  !  datum offsets, the under-relaxation store (allocated whether or not relaxation is on).
  subroutine Setup_IGLOO_Coupling ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nprim
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nfam, IGLOO_phase_prefix
    use hydraMI2_Chemistry,     only: chemIGLOO, IS_IGLOO, T_datum_ref, Load_Chemistry, Set_Datum_Offsets, &
                                      Report_Chemistry
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(obj_IGLOO),        intent(inout) :: IGLOO
    character(len=64), allocatable :: names(:)
    logical,           allocatable :: carries(:), released(:)
    real(R8),          allocatable :: href(:), Lv(:), hcar(:)
    integer  :: b, m, g, fam, nm
    real(R8) :: yield

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

    ! per material: its name, what it carries and hands the gas, its enthalpy at T_ref on the datum of its table, its
    ! latent heat, and the enthalpy of the material it carries at T_ref
    nm = size(IGLOO%material)
    allocate(names(nm), carries(nm), released(nm), href(nm), Lv(nm), hcar(nm))
    released = .false.
    do m = 1, nm
      associate ( mat => IGLOO%material(m) )
        names(m)   = mat%matName
        carries(m) = mat%carryShare > 0._R8 .or. mat%carryLaw > 0
        if ( carries(m) .and. .not. mat%resOn ) released(mat%carriedID) = .true.
        if ( mat%cpVariable ) then
          href(m) = Table_Enthalpy( mat%hTab, lbound(mat%hTab,1), T_datum_ref )
        else
          href(m) = mat%cp * T_datum_ref + mat%hOff
        endif
        Lv(m)   = mat%Lv
        hcar(m) = mat%cpCar * T_datum_ref + mat%hOffCar
      end associate
    enddo
    call Load_Chemistry( IS_IGLOO, IGLOO_phase_prefix, names, IGLOO%material(:)%evapSelect > 0, &
                         IGLOO%material(:)%combSelect > 0, carries, released )
    call Set_Datum_Offsets( IS_IGLOO, href, Lv, hcar )

    ! the yield IGLOO's carry law holds is the one of the carry record
    do m = 1, nm
      if ( .not. carries(m) ) cycle
      yield = IGLOO%material(m)%carryYield
      ! 1e-7: one formula with two tables of atomic weights that differ in their last digit stays within it
      if ( abs(chemIGLOO%mat(m)%carry%nuProduct - yield) > 1.e-7_R8*yield ) then
        write(*,'(A,ES15.7,A,ES15.7,A)') E_//'the carry record of IGLOO material '//trim(names(m))//' forms ', &
          chemIGLOO%mat(m)%carry%nuProduct, ' kg per kg, and its carry row gives yield = ', yield, &
          ': run ATLAS GPB again'
        error stop 1
      endif
    enddo

    call Setup_IGLOO_Gen( MOSE_domain, IGLOO, names )
    call Report_Chemistry( IS_IGLOO )

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


  !> The generation channels of IGLOO's materials: a process record whose product is born in an ICE phase, and the
  !  release record of a carried material handed at burnout that names one. ICE runs and holds the family, the
  !  coupling is not relaxed, the material gives d-gen, and IGLOO tracks the state its processed mass leaves at.
  subroutine Setup_IGLOO_Gen ( MOSE_domain, IGLOO, names )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use IGLOO_module,            only: obj_IGLOO
    use IGLOO_Mod_Exchange,      only: gen_diameter, track_processed
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydraMI2_Channels,       only: Open_Store
    use hydraMI2_Chemistry,      only: chemIGLOO, record_t
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(obj_IGLOO),        intent(inout) :: IGLOO
    character(len=*),       intent(in)    :: names(:)
    type(channel_t), allocatable :: found(:)
    type(record_t) :: rec
    character(len=:), allocatable :: who, why
    real(R8) :: dh
    integer  :: m, n, gs, b, c, nc, dims(3, size(MOSE_domain%blk))

    allocate(found(2*size(names)))
    nc = 0
    do m = 1, size(names)
      who = 'IGLOO material '//trim(names(m))
      do n = GEN_OF_PROCESS, GEN_OF_RELEASE
        if ( n == GEN_OF_PROCESS ) then
          rec = chemIGLOO%mat(m)%process
        elseif ( chemIGLOO%mat(m)%carry%on .and. .not. IGLOO%material(m)%resOn ) then
          rec = chemIGLOO%mat(IGLOO%material(m)%carriedID)%release
        else
          cycle
        endif
        if ( .not. rec%on .or. len_trim(rec%phase) == 0 ) cycle
        if ( .not. hydra_MI2_obj%run_ICE ) then
          write(*,'(A)') E_//who//' generates '//trim(rec%product)//', born in ICE, and ICE is not run'
          error stop 1
        endif
        if ( hydra_MI2_obj%relaxation ) then
          write(*,'(A)') E_//who//' generates '//trim(rec%product)//': relaxation = true blends two coupling '// &
            'steps of the gas''s share, and the generated mass is not blended; set relaxation = false'
          error stop 1
        endif
        if ( .not. gen_diameter(IGLOO, m) > 0._R8 ) then
          write(*,'(A)') E_//who//' generates '//trim(rec%product)//' and its combustion row gives no d-gen, the '// &
            'diameter each generated drop is born with'
          error stop 1
        endif
        call track_processed( IGLOO, m, why )
        if ( len(why) > 0 ) then
          write(*,'(A)') E_//who//' generates '//trim(rec%product)//': '//why
          error stop 1
        endif
        call Gen_Family( who, rec%product, rec%phase, rec%hProduct, gs, dh )
        nc = nc + 1
        found(nc)%name    = 'igloo-gen'
        found(nc)%from    = m
        found(nc)%to      = gs
        found(nc)%trigger = n
        found(nc)%nu      = rec%nuProduct
        found(nc)%dh      = dh
        found(nc)%d       = gen_diameter(IGLOO, m)
        found(nc)%asRate  = handRates
        write(*,'(A,I0,A,ES15.7,A,ES12.4,A)') ' [MI2] '//who//' generates '//trim(rec%product)// &
          trim(merge('                        ', ' from what it carries   ', n == GEN_OF_PROCESS))//' (ICE family ', &
          gs, '): ', rec%nuProduct, ' kg per kg, d-gen ', found(nc)%d, ' m'
      enddo
    enddo
    allocate(gen(nc))
    do c = 1, nc
      gen(c) = found(c)
    enddo
    genIGOn = nc > 0
    ! each channel's store: what a MOSE block's iterations make (mass, momentum (3), energy, number)
    do b = 1, size(MOSE_domain%blk)
      dims(:,b) = MOSE_domain%blk(b)%dim(1:3)
    enddo
    do c = 1, size(gen)
      call Open_Store( gen(c), dims, 6, timed=.false., rated=.false. )
    enddo
    if ( genIGOn .and. handRates ) &
      write(*,'(A)') ' [MI2] IGLOO''s generated mass reaches ICE as rates: births over ICE''s own steps'
    if ( genIGOn .and. .not. handRates ) &
      write(*,'(A)') ' [MI2] IGLOO''s generated mass reaches ICE as amounts: what each MOSE block made, after the block'

  end subroutine Setup_IGLOO_Gen


  !> The ICE family a generated mass of material `product` is born in (ICE's birth_family: one family, MK, inert) and
  !  the generated material's datum offset there: its formation enthalpy hProduct less ICE's energy of it at T_ref.
  !  phase: the ICE phase the record names.
  subroutine Gen_Family ( owner, product, phase, hProduct, gs, dh )
    use ICE_Global_m,       only: ICE_phase_prefix, mat_of
    use ICE_Mod_Exchange,   only: birth_family
    use hydraMI2_Chemistry, only: T_datum_ref
    implicit none
    character(len=*), intent(in)  :: owner, product, phase
    real(R8),         intent(in)  :: hProduct
    integer,          intent(out) :: gs
    real(R8),         intent(out) :: dh
    character(len=:), allocatable :: why

    if ( trim(phase) /= trim(ICE_phase_prefix) ) then
      write(*,'(A)') E_//'the mass '//owner//' generates is born in the phase '//trim(phase)//', and ICE runs '// &
        'the phase '//trim(ICE_phase_prefix)//': run ATLAS GPB again'
      error stop 1
    endif
    gs = birth_family( product, why )
    if ( gs == 0 ) then
      write(*,'(A)') E_//'the mass '//owner//' generates cannot be born in ICE: '//why
      error stop 1
    endif
    dh = hProduct - ICE_Energy( mat_of(gs), T_datum_ref )

  end subroutine Gen_Family


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


  !> Push the MOSE gas state into ICE (set_gas), interior cells. Opens its own parallel region.
  subroutine MOSExICE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use ICE_Mod_Exchange,      only: set_gas
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
        call set_gas( ICE_domain%blk(b), i, j, k, rho, vel, T, R, gam, mu, kappa )
      enddo ; enddo ; enddo
      !$OMP END PARALLEL DO
    enddo

  end subroutine MOSExICE


  ! ---------------------------------------------------------------------------------------------
  ! Condensed solvers -> MOSE: the sources, from MOSE's per-RK-stage callback. Inside MOSE's parallel region:
  ! orphaned !$OMP DO only, no ASSOCIATE, and cell kernels that write only their arguments.
  ! ---------------------------------------------------------------------------------------------

  !> Inject IGLOO's reaction on the gas into MOSE's residual (R = R - S).
  !    source branch : IGLOO_Cell (absolute rates), under the consumed-species floor guard
  !    euler  branch : getSourceTerms is per unit volume (*vol); momentum and energy only
  !  iter of iMax is the MOSE iteration within its block and relCoeff the cycle's under-relaxation coefficient
  !  (1: the source is injected as it is).
  subroutine IGLOOxMOSE ( MOSE_domain, IGLOO, eulerian, iter, iMax, relCoeff )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nu, nw, np, nprim
    use MOSE_Config_Types_m,    only: obj_chemistry
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nfam
    use hydraMI2_Chemistry,     only: chemIGLOO
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: eulerian
    integer,                intent(in)    :: iter, iMax
    real(R8),               intent(in)    :: relCoeff
    logical  :: blend
    integer  :: b, i, j, k, m
    ! res mirrors blk%R(1:nprim): species 1:nsc, momentum nu:nw, energy np, then soot/passive/RANS (left at zero)
    real(R8) :: res(nprim), sources(4), datum
    real(R8) :: rhoi(nsc), rho, u(3), T, R, gam, kappa, mu, coeff

    call Relaxation_Weight ( iter, iMax, relCoeff, coeff, blend )

    do b = 1, size(MOSE_domain%blk)

      !$OMP DO COLLAPSE(3) PRIVATE(res,sources,datum,rhoi,rho,u,T,R,gam,kappa,mu,i,j,k,m)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)

        res = 0._R8
        if (eulerian) then
          call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), rhoi, u, rho, R, gam, T, mu, kappa )
          sources = 0._R8
          do m = 1, nfam
            sources = sources + IGLOO%getSourceTerms( [rho,u,T,gam,R,mu,kappa],      &
                                    IGLOO%euler(b,m)%velocity(:,i,j,k),              &
                                    IGLOO%euler(b,m)%temperature(i,j,k),             &
                                    IGLOO%euler(b,m)%density(i,j,k),                 &
                                    IGLOO%euler(b,m)%np(i,j,k), famMatID(m) )
          enddo
          res(nu:nw) = - sources(1:3)
          res(np)    = - sources(4)
          res = res * MOSE_domain%blk(b)%vol(i,j,k)
        else
          call IGLOO_Cell ( IGLOO, b, i, j, k, res, datum )
          if (chemIGLOO%active) &
            call Floor_Guard ( res(1:nsc), MOSE_domain%blk(b)%P(1:nsc,i,j,k), MOSE_domain%blk(b)%vol(i,j,k), &
                               MOSE_domain%blk(b)%dtlocal(i,j,k) * obj_chemistry%strangcoeff )
        endif

        if (blend) then
          res(1:np) = coeff*res(1:np) + (1._R8-coeff)*oldResidual(b)%res(1:np,i,j,k)
          oldResidual(b)%res(:,i,j,k) = res
        endif

        MOSE_domain%blk(b)%R(:,i,j,k) = MOSE_domain%blk(b)%R(:,i,j,k) - res

      enddo ; enddo ; enddo
      !$OMP END DO

    enddo

  end subroutine IGLOOxMOSE


  !> Inject ICE's reaction on the gas into MOSE's residual: what each family gives the gas over its last step
  !  (ICE_Cell) and the removals ICE's burn-out made (Removal_Cell: a steady run, the rate of the store; a
  !  time-accurate one, this MOSE iteration's share of it as a rate over its step), under the consumed-species floor
  !  guard, times the cell volume.
  subroutine ICExMOSE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nprim
    use MOSE_Config_Types_m,    only: obj_chemistry
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use ICE_Mod_Exchange,       only: exchange_t, exchange_rate
    use hydraMI2_Chemistry,     only: chemICE
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(ICE_domain_type),  intent(in)    :: ICE_domain
    type(exchange_t) :: x
    integer  :: b, g, i, j, k, imat, nshare
    real(R8) :: r(nprim), dtIter, datum, gr(6)
    logical  :: mapped, consumes

    ! the families whose material consumes a gas species share what a cell holds of it
    nshare = count([(Consumes_Gas(mat_of(g)), g = 1, ngroups)])

    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups

        imat     = mat_of(g)
        mapped   = Converts(imat)
        consumes = Consumes_Gas(imat)

        !$OMP DO COLLAPSE(3) PRIVATE(x,r,i,j,k,dtIter,datum,gr)
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)

          r = 0._R8
          x = exchange_rate( ICE_domain%blk(b), g, i, j, k )
          call ICE_Cell ( imat, mapped, x, r, datum )
          if (mapped) then
            if (removed(g)%asRate) then
              if (allocated(removed(g)%blk(b)%r)) then
                gr = removed(g)%blk(b)%r(:,i,j,k)
                call Removal_Cell ( imat, gr, r, datum )
              endif
            elseif (removed(g)%share > 0._R8) then
              dtIter = MOSE_domain%blk(b)%dtlocal(i,j,k)
              if (dtIter > 0._R8) then
                gr = removed(g)%blk(b)%a(:,i,j,k)*removed(g)%share/dtIter
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
      if (Converts(imat)) Consumes_Gas = any(chemICE%mat(imat)%process%nu(1:nsc) < 0._R8)
    end function Consumes_Gas

  end subroutine ICExMOSE


  !> ICE material imat converts: it has a process record, through which the gas receives its mass.
  logical function Converts ( imat )
    use hydraMI2_Chemistry, only: chemICE
    implicit none
    integer, intent(in) :: imat

    Converts = .false.
    if (allocated(chemICE%mat)) Converts = chemICE%mat(imat)%process%on

  end function Converts


  !> The mass IGLOO material m's parcels consume in a cell (x: its exchange_rate) by each of its records [kg/s]: F(1)
  !  through its process record, F(2) through its carry record (the part that forms what they carry), F(3) the carried
  !  material they hand the gas at burnout, through that material's release record. The parcels lose F(1) + F(2) of
  !  their own material, gain the carried material deposited and lose the one handed, so the net mass the gas gets is
  !  x%mass. Thread-safe.
  pure function Material_Flows ( m, x ) result(F)
    use IGLOO_Mod_Exchange, only: exchange_t
    use hydraMI2_Chemistry, only: chemIGLOO
    implicit none
    integer,          intent(in) :: m
    type(exchange_t), intent(in) :: x
    real(R8) :: F(3)

    F = [x%mass, 0._R8, 0._R8]
    if (.not. chemIGLOO%mat(m)%carry%on) return
    F(2) = x%carried_in/chemIGLOO%mat(m)%carry%nuProduct
    F(1) = x%mass + x%carried_in - x%carried_out - F(2)
    F(3) = x%carried_out

  end function Material_Flows


  !> What IGLOO's parcels hand the gas in geometry cell (b, i, j, k), added to r (the residual's layout) as absolute
  !  rates: the momentum and energy of every parcel (gas_rate); for each material that converts, its gas species
  !  through its records with its datum offsets; less what it generates in ICE, born there. datum: the offsets' part of
  !  the energy. Thread-safe.
  subroutine IGLOO_Cell ( IGLOO, b, i, j, k, r, datum )
    use MOSE_Global_m,      only: nu, nw, np
    use IGLOO_module,       only: obj_IGLOO
    use IGLOO_Mod_Exchange, only: exchange_t, exchange_rate, gas_rate
    use hydraMI2_Chemistry, only: chemIGLOO, Gas_Of
    implicit none
    type(obj_IGLOO), intent(in)    :: IGLOO
    integer,         intent(in)    :: b, i, j, k
    real(R8),        intent(inout) :: r(:)
    real(R8),        intent(out)   :: datum
    type(exchange_t) :: x
    real(R8) :: q(4), F(3), sr(6)
    integer  :: m, ic, c

    q = gas_rate( IGLOO, b, i, j, k )
    r(nu:nw) = r(nu:nw) + q(1:3)
    r(np)    = r(np) + q(4)
    datum    = 0._R8
    if (.not. chemIGLOO%active) return

    do m = 1, size(chemIGLOO%mat)
      x = exchange_rate( IGLOO, m, b, i, j, k )
      F = Material_Flows( m, x )
      call Gas_Of( chemIGLOO%mat(m)%process, F(1), r )
      call Gas_Of( chemIGLOO%mat(m)%carry,   F(2), r )
      datum = datum + chemIGLOO%mat(m)%dh*(F(1) + F(2))
      if (chemIGLOO%mat(m)%carry%on) then
        ic = IGLOO%material(m)%carriedID
        call Gas_Of( chemIGLOO%mat(ic)%release, F(3), r )
        datum = datum - chemIGLOO%mat(m)%carry%dhProduct*(x%carried_in - x%carried_out)
      endif
    enddo

    ! the mass a material generates is born in its ICE family: it leaves the gas's books
    do c = 1, size(gen)
      call Gen_Packet( gen(c), IGLOO, b, i, j, k, sr )
      r(nu:nw) = r(nu:nw) - sr(2:4)
      r(np)    = r(np) - sr(5)
      datum    = datum - gen(c)%dh*sr(1)
    enddo
    r(np) = r(np) + datum

  end subroutine IGLOO_Cell


  !> What one ICE family gives the gas in a cell, per unit volume (x: its exchange_rate, or a rate of its removals),
  !  into r (the residual's layout): the momentum and energy; for a material that converts (mapped) the mass it loses,
  !  mass + generated, through its process record with its datum offset, less the generated mass's offset (it stays
  !  condensed). datum: the offsets' part of the energy. Thread-safe.
  subroutine ICE_Cell ( imat, mapped, x, r, datum )
    use MOSE_Global_m,      only: nu, nw, np
    use ICE_Mod_Exchange,   only: exchange_t
    use hydraMI2_Chemistry, only: chemICE, Gas_Of
    implicit none
    integer,          intent(in)    :: imat
    logical,          intent(in)    :: mapped
    type(exchange_t), intent(in)    :: x
    real(R8),         intent(inout) :: r(:)
    real(R8),         intent(out)   :: datum
    real(R8) :: F

    r(nu:nw) = r(nu:nw) + x%momentum
    r(np)    = r(np) + x%energy
    datum    = 0._R8
    if (.not. mapped) return

    F = x%mass + x%generated
    call Gas_Of( chemICE%mat(imat)%process, F, r )
    datum = chemICE%mat(imat)%dh*F - chemICE%mat(imat)%process%dhProduct*x%generated
    r(np) = r(np) + datum

  end subroutine ICE_Cell


  !> A rate g of ICE's removals of material imat, per unit volume ([gas mass, momentum (3), energy, generated mass]),
  !  added to r as ICE_Cell adds an exchange, its offsets added to datum too. Thread-safe.
  subroutine Removal_Cell ( imat, g, r, datum )
    use ICE_Mod_Exchange, only: exchange_t
    implicit none
    integer,  intent(in)    :: imat
    real(R8), intent(in)    :: g(6)
    real(R8), intent(inout) :: r(:), datum
    type(exchange_t) :: x
    real(R8) :: d

    x%mass      = g(1)
    x%momentum  = g(2:4)
    x%energy    = g(5)
    x%generated = g(6)
    call ICE_Cell ( imat, .true., x, r, d )
    datum = datum + d

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

  !> IGLOO's sources: under the source mode what IGLOO_Cell reads (gas_rate, and every material's exchange_rate when
  !  one converts); under the Eulerian mode the fields it turns into sources (velocity, temperature, density, np).
  subroutine Check_IGLOO_Sources ( MOSE_domain, IGLOO, eulerian )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nfam
    use IGLOO_Mod_Exchange,     only: exchange_t, exchange_rate, gas_rate
    use MOSE_Finite_m,          only: Not_Finite
    use hydraMI2_Chemistry,     only: chemIGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    logical,                intent(in) :: eulerian
    character(len=*), parameter :: EUL = 'Eulerian field'
    type(exchange_t) :: x
    character(len=32) :: what
    integer :: b, m, i, j, k, n(3)

    do b = 1, size(MOSE_domain%blk)
      n = MOSE_domain%blk(b)%dim(1:3)
      if (eulerian) then
        do m = 1, nfam
          associate ( e => IGLOO%euler(b,m) )
          if (allocated(e%velocity))    call Refuse_4( 'IGLOO', EUL, b, m, 'velocity', 1, e%velocity(:,1:n(1),1:n(2),1:n(3)) )
          if (allocated(e%temperature)) call Refuse_3( 'IGLOO', EUL, b, m, 'temperature', e%temperature(1:n(1),1:n(2),1:n(3)) )
          if (allocated(e%density))     call Refuse_3( 'IGLOO', EUL, b, m, 'density', e%density(1:n(1),1:n(2),1:n(3)) )
          if (allocated(e%np))          call Refuse_3( 'IGLOO', EUL, b, m, 'np', e%np(1:n(1),1:n(2),1:n(3)) )
          end associate
        enddo
        cycle
      endif
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)
        if (any(Not_Finite(gas_rate(IGLOO, b, i, j, k)))) &
          call Refuse_At( 'IGLOO', 'source', b, 0, 'momentum or energy', [i, j, k], 0._R8 )
        if (.not. chemIGLOO%active) cycle
        do m = 1, size(IGLOO%material)
          x = exchange_rate( IGLOO, m, b, i, j, k )
          if (.not. any(Not_Finite([x%mass, x%carried_in, x%carried_out, x%processed, x%u, x%T]))) cycle
          write(what,'(A,I0)') 'mass of material ', m
          call Refuse_At( 'IGLOO', 'source', b, 0, trim(what), [i, j, k], x%mass )
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine Check_IGLOO_Sources


  !> ICE's sources, the slots exchange_rate reads: the momentum source(2:4) and the energy source(nbase-1) of every
  !  family, the mass source(1) of a material that converts. And a family whose material has no record in the products
  !  file that loses mass: ICExMOSE would not hand that mass to the gas.
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
    use MOSE_Global_m,          only: nsc, np, nprim
    use IGLOO_module,           only: obj_IGLOO
    use hydraMI2_Chemistry,     only: chemIGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    type(handed_t) :: h
    real(R8), allocatable :: w(:,:,:,:)
    integer  :: b, i, j, k
    real(R8) :: r(nprim)

    if (chemIGLOO%active) allocate(h%mass(nsc), source=0._R8)
    do b = 1, size(MOSE_domain%blk)
      if (.not. allocated(IGLOO%source(b)%sourceMom)) cycle
      call Cells( MOSE_domain%blk(b)%dim, w )
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        r = 0._R8
        call IGLOO_Cell( IGLOO, b, i, j, k, r, w(0,i,j,k) )
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
    use ICE_Global_m,           only: ngroups, mat_of
    use ICE_Mod_Exchange,       only: exchange_rate
    use hydraMI2_Chemistry,     only: chemICE
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(ICE_domain_type),  intent(in) :: ICE_domain
    type(handed_t) :: h
    real(R8), allocatable :: w(:,:,:,:)
    integer  :: b, g, i, j, k
    real(R8) :: r(nprim), datum, vol

    if (chemICE%active) allocate(h%mass(nsc), source=0._R8)
    do b = 1, size(MOSE_domain%blk)
      call Cells( MOSE_domain%blk(b)%dim, w )
      do g = 1, ngroups
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          r = 0._R8
          call ICE_Cell( mat_of(g), Converts(mat_of(g)), exchange_rate(ICE_domain%blk(b), g, i, j, k), r, datum )
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
    use hydraMI2_Chemistry,     only: chemICE
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    real(R8),               intent(in) :: share
    logical, optional,      intent(in) :: rate
    type(handed_t) :: h
    real(R8), allocatable :: w(:,:,:,:)
    integer  :: b, g, i, j, k
    real(R8) :: r(nprim), datum, vol, gr(6)
    logical  :: asRate

    if (.not. chemICE%active) return
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
    use hydraMI2_Chemistry,     only: chemICE
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
    if (.not. chemICE%active) return
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

  !> After an ICE iteration: what each family's burn-out removed (take_removed) moved into its channel, and under
  !  rates each converting family's dt added, cell by cell, to the time its store is being made over. A
  !  time-accurate ICE: each family's dt*source, dt*heat, dt*gen_rate and dt added to the block's mean channel.
  subroutine After_ICE_Step ( ICE_domain )
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use ICE_Config_Types_m,     only: ICE_time_scheme => obj_time_scheme
    use ICE_Mod_Exchange,       only: exchange_t, take_removed
    use hydraMI2_Channels,      only: Open_Store, Clock
    use hydraMI2_Chemistry,     only: chemICE
    implicit none
    type(ICE_domain_type), intent(inout) :: ICE_domain
    type(exchange_t) :: x
    integer :: b, g, i, j, k, ns, dims(3, ICE_domain%nb)

    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        do k = 1, ICE_domain%blk(b)%dim(3)
        do j = 1, ICE_domain%blk(b)%dim(2)
        do i = 1, ICE_domain%blk(b)%dim(1)
          call take_removed( ICE_domain%blk(b), g, i, j, k, x )
          removed(g)%blk(b)%a(:,i,j,k) = removed(g)%blk(b)%a(:,i,j,k) + [x%mass, x%momentum, x%energy, x%generated]
        enddo; enddo; enddo
        if (handRates .and. chemICE%active .and. Converts(mat_of(g))) &
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
        call Open_Store( iceMean(g), dims, size(ICE_domain%blk(1)%cond_phase(g)%source, 1) + 6, timed=.true., &
                         rated=.false. )
      enddo
    endif
    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        ns = size(ICE_domain%blk(b)%cond_phase(g)%source, 1)
        do k = 1, ICE_domain%blk(b)%dim(3)
        do j = 1, ICE_domain%blk(b)%dim(2)
        do i = 1, ICE_domain%blk(b)%dim(1)
          associate (a => iceMean(g)%blk(b)%a, dt => ICE_domain%blk(b)%cond_phase(g)%dt(i,j,k))
            a(1:ns,i,j,k)      = a(1:ns,i,j,k)      + dt*ICE_domain%blk(b)%cond_phase(g)%source(:,i,j,k)
            a(ns+1,i,j,k)      = a(ns+1,i,j,k)      + dt*ICE_domain%blk(b)%cond_phase(g)%heat(i,j,k)
            a(ns+2:ns+6,i,j,k) = a(ns+2:ns+6,i,j,k) + dt*ICE_domain%blk(b)%cond_phase(g)%gen_rate(:,i,j,k)
          end associate
        enddo; enddo; enddo
        call Clock( iceMean(g), b, ICE_domain%blk(b)%cond_phase(g)%dt )
      enddo
    enddo
    meanSteps = meanSteps + 1

  end subroutine After_ICE_Step


  !> After a time-accurate ICE block: where it made more than one step, each family's source, heat and gen_rate become
  !  the block's rates, its mean channel's store over its time; the channels emptied. (exchange_rate then reads the
  !  block's mean.)
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
            ICE_domain%blk(b)%cond_phase(g)%source(:,i,j,k)   = iceMean(g)%blk(b)%a(1:ns,i,j,k)/t
            ICE_domain%blk(b)%cond_phase(g)%heat(i,j,k)       = iceMean(g)%blk(b)%a(ns+1,i,j,k)/t
            ICE_domain%blk(b)%cond_phase(g)%gen_rate(:,i,j,k) = iceMean(g)%blk(b)%a(ns+2:ns+6,i,j,k)/t
          enddo; enddo; enddo
        enddo
      endif
      call Empty( iceMean(g) )
    enddo
    meanSteps = 0

  end subroutine ICE_Block_Mean


  ! ---------------------------------------------------------------------------------------------
  ! IGLOO's generated mass, from the gas's books to ICE
  ! ---------------------------------------------------------------------------------------------

  !> The mass generation channel ch carries in geometry cell (b, i, j, k), as a rate [per s, the cell's]: ICE's birth
  !  packet (mass, momentum (3), energy on ICE's datum with the kinetic part, number at the channel's d), born with the
  !  velocity and temperature the processed mass of the material left at. nostate: born where IGLOO tracked no state
  !  (at rest, at T_datum_ref). Thread-safe.
  subroutine Gen_Packet ( ch, IGLOO, b, i, j, k, rate, nostate )
    use IGLOO_module,       only: obj_IGLOO
    use IGLOO_Mod_Exchange, only: exchange_t, exchange_rate
    use ICE_Mod_Exchange,   only: birth_packet
    use hydraMI2_Chemistry, only: T_datum_ref
    implicit none
    type(channel_t),   intent(in)  :: ch
    type(obj_IGLOO),   intent(in)  :: IGLOO
    integer,           intent(in)  :: b, i, j, k
    real(R8),          intent(out) :: rate(6)
    logical, optional, intent(out) :: nostate
    type(exchange_t) :: x
    real(R8) :: F(3), md, u(3), T

    rate = 0._R8
    if (present(nostate)) nostate = .false.
    x  = exchange_rate( IGLOO, ch%from, b, i, j, k )
    F  = Material_Flows( ch%from, x )
    md = ch%nu*merge(F(1), F(3), ch%trigger == GEN_OF_PROCESS)
    if (md == 0._R8) return
    if (x%tracked) then
      u = x%u
      T = x%T
    else
      u = 0._R8
      T = T_datum_ref
      if (present(nostate)) nostate = .true.
    endif
    rate = birth_packet( ch%to, md, u, T, ch%d )

  end subroutine Gen_Packet


  !> After a MOSE iteration: each generation channel's rate times the iteration's dtlocal, added cell by cell to its
  !  store (the gas received the same rates less, over the same step). A cell where a material's generated mass had no
  !  parcel's state counts once. Outside any parallel region.
  subroutine IGLOO_Gen_Accumulate ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    integer  :: b, i, j, k, c, mNoState
    real(R8) :: sr(6), dt
    logical  :: nostate

    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        dt = MOSE_domain%blk(b)%dtlocal(i,j,k)
        mNoState = 0
        do c = 1, size(gen)
          call Gen_Packet( gen(c), IGLOO, b, i, j, k, sr, nostate )
          ! the channels of a material are consecutive
          if (nostate .and. gen(c)%from /= mNoState) then
            nGenNoState = nGenNoState + 1
            mNoState    = gen(c)%from
          endif
          if (sr(1) == 0._R8) cycle
          gen(c)%blk(b)%a(:,i,j,k) = gen(c)%blk(b)%a(:,i,j,k) + sr*dt
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine IGLOO_Gen_Accumulate


  !> After a MOSE block: the generated mass its iterations made (IGLOO_Gen_Accumulate), which the gas did not get,
  !  handed to each ICE family a generation channel feeds. Amounts: born at once (add_birth, per volume). Rates: the
  !  rates of IGLOO's sources as they stand (srcOn; else 0), per volume, the family's births in force (add_birth_rate).
  !  sm: the totals, energies on the gas's formation scale; the stores emptied. Outside any parallel region.
  subroutine IGLOO_Gen_Deliver ( MOSE_domain, ICE_domain, IGLOO, srcOn, sm )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    use ICE_Global_m,          only: ngroups
    use ICE_Mod_Exchange,      only: add_birth, add_birth_rate, clear_import_rates
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: srcOn
    type(gen_t),            intent(out)   :: sm
    integer  :: b, g, i, j, k, c
    real(R8) :: vol, ac(6), sr(6), pk(6), dhg
    logical  :: asRate

    sm%nNoState = nGenNoState
    nGenNoState = 0
    asRate = gen(1)%asRate
    do g = 1, ngroups
      if (.not. any(gen(:)%to == g)) cycle
      if (asRate) call clear_import_rates( ICE_domain, g )
      ! the generated material's datum offset, whichever channel feeds it
      do c = 1, size(gen)
        if (gen(c)%to == g) dhg = gen(c)%dh
      enddo
      do b = 1, size(MOSE_domain%blk)
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          ! rates: the births of family g from IGLOO's sources as they stand
          if (asRate .and. srcOn) then
            sr = 0._R8
            do c = 1, size(gen)
              if (gen(c)%to /= g) cycle
              call Gen_Packet( gen(c), IGLOO, b, i, j, k, pk )
              sr = sr + pk
            enddo
            if (sr(1) > 0._R8) call add_birth_rate( ICE_domain%blk(b), g, i, j, k, sr/vol )
            sm%rate(1:4) = sm%rate(1:4) + sr(1:4)
            sm%rate(5)   = sm%rate(5) + sr(5) + dhg*sr(1)
            sm%rate(6)   = sm%rate(6) + sr(6)
          endif
          ! the amounts of the block: out of the gas's books, and into ICE's state when they go as amounts
          ac = 0._R8
          do c = 1, size(gen)
            if (gen(c)%to == g) ac = ac + gen(c)%blk(b)%a(:,i,j,k)
          enddo
          if (all(ac == 0._R8)) cycle
          if (.not. asRate) call add_birth( ICE_domain%blk(b), g, i, j, k, ac/vol )
          sm%amount(1:4) = sm%amount(1:4) + ac(1:4)
          sm%amount(5)   = sm%amount(5) + ac(5) + dhg*ac(1)
          do c = 1, size(gen)
            if (gen(c)%to == g) gen(c)%blk(b)%a(:,i,j,k) = 0._R8
          enddo
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine IGLOO_Gen_Deliver

end module hydraMI2_Mod_Coupling
