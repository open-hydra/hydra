!> The exchange between MOSE and the two condensed solvers of hydra-MI2.
!
!  Within a coupling step:
!    MOSExIGLOO, MOSExICE : the gas state a condensed solver steps on
!    IGLOOxMOSE, ICExMOSE : the sources it hands back, into MOSE's residual at every Runge-Kutta stage
!  What a solver hands the gas in one cell is written once, in the cell kernels IGLOO_Mass_To_Gas, ICE_Cell and
!  Removal_Cell: the injection applies them, IGLOO_Handed, ICE_Handed and the removals' totals sum them for the
!  books and the report (hydraMI2_Report).
!
!  Between two MOSE blocks, what lands outside a solver's step:
!    ICE's removals    (cond_phase%gone): what a converting material lost and the gas has not received yet
!    ICE's mean source : a time-accurate ICE block of several steps hands the gas its steps' dt-weighted mean
!    IGLOO's smoke     : the condensed product of IGLOO's burning, born in an ICE family instead of the gas
!  They go as rates when MOSE or ICE steps in pseudo-time (handRates: a converged steady field then holds per unit
!  time), as amounts when both are time-accurate.
module hydraMI2_Mod_Coupling
  use iso_fortran_env, only: R8 => real64

  implicit none
  private

  public :: source_state_t, handed_t, removals_t, smoke_t
  public :: Setup_ICE_Coupling, Setup_IGLOO_Coupling, Set_IGLOO_Switches
  public :: Gas_Field_For_IGLOO, MOSExIGLOO, MOSExICE, IGLOOxMOSE, ICExMOSE
  public :: Check_IGLOO_Sources, Check_ICE_Sources
  public :: IGLOO_Handed, ICE_Handed, Removals_Handed
  public :: Open_Removals, Close_Removals, Gone_Rate_Setup, ICE_Step_Clocks, ICE_Block_Mean
  public :: IGLOO_Smoke_Accumulate, IGLOO_Smoke_Deliver

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

  !> Amounts: the share of ICE's stored removals each MOSE iteration of the current block hands the gas, 1/MOSE-iter
  !  while a block runs, 0 outside one.
  real(R8) :: goneShare = 0._R8

  !> Rates, per block: a(5, ICE family, cell) the store ICE made over the time t(ICE family, cell) it was made in, as
  !  a rate per unit volume (mass, momentum (3), energy), kept while ICE does not step, as its source is. Restored by
  !  the hand-off file of a restart.
  type, public :: gone_rate_t
    real(R8), allocatable :: a(:,:,:,:,:), t(:,:,:,:)
  end type gone_rate_t
  type(gone_rate_t), allocatable, public :: goneRate(:)

  !> A time-accurate ICE block of more than one step: per block and family, the sums over its steps of dt*source,
  !  dt*heat, dt*birth and dt, cell by cell (ICE_Step_Clocks), turned into the block's rates after it (ICE_Block_Mean).
  type :: src_sum_t
    real(R8), allocatable :: s(:,:,:,:), h(:,:,:), t(:,:,:), b(:,:,:,:)
  end type src_sum_t
  type(src_sum_t), allocatable :: srcSum(:,:)
  integer :: srcSteps = 0

  !> The smoke of IGLOO's burning materials: per material, for the smoke of its combustion record (1) and the smoke
  !  its cap oxide becomes when handed at burnout (2), the ICE family it is born in (0: none), kg per kg, the smoke's
  !  datum offset; the material's d-smoke. smokeIGOn: one at least.
  type :: igloo_smoke_t
    integer  :: fam(2) = 0
    real(R8) :: nu(2) = 0._R8, dh(2) = 0._R8, dSmoke = 0._R8
  end type igloo_smoke_t
  type(igloo_smoke_t), allocatable :: smokeIG(:)
  logical, protected, public :: smokeIGOn = .false.

  !> IGLOO's smoke the current MOSE block made, per block (6, ICE family, cell): mass, momentum (3), energy, number;
  !  each smoke family's datum offset; the cells where a smoke was born with no parcel's state.
  type :: smoke_block_t
    real(R8), allocatable :: a(:,:,:,:,:)
  end type smoke_block_t
  type(smoke_block_t), allocatable :: smokeAcc(:)
  real(R8),            allocatable :: dhSmokeFam(:)
  integer :: nSmokeNoState = 0

  character(len=*), parameter :: E_ = ' [ERROR] MI2: '

contains

  ! ---------------------------------------------------------------------------------------------
  ! Setup
  ! ---------------------------------------------------------------------------------------------

  !> ICE's side of the coupler, after ICE%setup: how what lands between steps goes (handRates), material -> gas
  !  species of ICE's phase from the file ATLAS writes, the smoke families, and the datum offsets. ICE credits a
  !  material its energy at T_ref on its own datum: its table's, or cp*T.
  subroutine Setup_ICE_Coupling ()
    use MOSE_Config_Types_m,    only: MOSE_time_scheme => obj_time_scheme
    use ICE_Global_m,           only: ICE_phase_prefix
    use ICE_Lib_Properties,     only: mat_e, mat_cp_const, mat_lv
    use ICE_Config_Types_m,     only: COMB_BECKSTEAD, ICE_condensed => obj_condensed, ICE_time_scheme => obj_time_scheme
    use hydraMI2_Stoichiometry, only: Load_Species_Map, Report_Stoichiometry, Build_Datum_Offsets, Smoke_Record, &
                                      Set_Smoke_Datum, MAP_ICE, T_datum_ref
    implicit none
    character(len=64), allocatable :: matnames(:)
    logical,           allocatable :: evap(:), burn(:), beck(:), Lvconst(:)
    real(R8),          allocatable :: eref(:), Lvmat(:), qmat(:), Thot(:), ehot(:), beta(:)
    character(len=64) :: sname, sphase
    real(R8) :: snu, sh, dh
    integer  :: m, nm, gs

    handRates = .not. (MOSE_time_scheme%time_accurate .and. ICE_time_scheme%time_accurate)

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

  end subroutine Setup_ICE_Coupling


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
    call Setup_IGLOO_Smoke( IGLOO, matnames )
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


  !> The smoke of IGLOO's burning materials (the tail of a combustion record, and of the release record of a cap oxide
  !  handed at burnout): ICE runs and holds its family, the coupling is not relaxed, the material gives d-smoke.
  subroutine Setup_IGLOO_Smoke ( IGLOO, matnames )
    use IGLOO_module,            only: obj_IGLOO
    use IGLOO_variables,         only: IGLOO_phase_prefix
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydraMI2_Stoichiometry,  only: Smoke_Record, MAP_IGLOO
    implicit none
    type(obj_IGLOO),  intent(in) :: IGLOO
    character(len=*), intent(in) :: matnames(:)
    character(len=64)  :: sname, sphase
    character(len=256) :: file
    character(len=:), allocatable :: who
    real(R8) :: snu, sh, dh
    integer  :: m, p, gs, n
    logical  :: has

    allocate(smokeIG(max(size(matnames),1)))
    file = 'INPUT/'//trim(IGLOO_phase_prefix)//'products.txt'
    do m = 1, size(matnames)
      who = 'the IGLOO material '//trim(matnames(m))
      ! n = 1: the smoke of its combustion record; n = 2: the smoke of its cap oxide, handed at burnout
      do n = 1, 2
        if ( n == 1 ) then
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
        smokeIG(m)%fam(n) = gs
        smokeIG(m)%nu(n)  = snu
        smokeIG(m)%dh(n)  = dh
        smokeIG(m)%dSmoke = IGLOO%material(m)%dSmoke
        smokeIGOn = .true.
        write(*,'(A,I0,A,ES15.7,A,ES12.4,A,ES15.7,A)') ' [MI2] IGLOO material '//trim(matnames(m))//' '// &
          trim(merge('burns into' // repeat(' ', 17), 'hands its cap at burnout to', n == 1))//' the smoke '// &
          trim(sname)//' (ICE family ', gs, '): ', snu, ' kg per kg, d-smoke ', smokeIG(m)%dSmoke, &
          ' m, datum offset ', dh, ' J/kg'
      enddo
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


  ! ---------------------------------------------------------------------------------------------
  ! Condensed solvers -> MOSE: the sources, from MOSE's per-RK-stage callback. Inside MOSE's parallel region:
  ! orphaned !$OMP DO only, no ASSOCIATE, and cell kernels that write only their arguments.
  ! ---------------------------------------------------------------------------------------------

  !> Inject IGLOO's reaction on the gas into MOSE's residual (R = R - S).
  !    source branch : source{Mass,Mom,En} are absolute rates; their mass through IGLOO_Mass_To_Gas, under the
  !                    consumed-species floor guard
  !    euler  branch : getSourceTerms is per unit volume (*vol); momentum and energy only
  !  iter of iMax is the MOSE iteration within its block and relCoeff the cycle's under-relaxation coefficient
  !  (1: the source is injected as it is).
  subroutine IGLOOxMOSE ( MOSE_domain, IGLOO, eulerian, iter, iMax, relCoeff )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nu, nw, np, nprim
    use MOSE_Config_Types_m,    only: obj_chemistry
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nfam
    use hydraMI2_Stoichiometry, only: productsIGLOO
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
          res(nu:nw) = IGLOO%source(b)%sourceMom(:,i,j,k)
          res(np)    = IGLOO%source(b)%sourceEn(i,j,k)
          if (productsIGLOO%active) then
            call IGLOO_Mass_To_Gas ( IGLOO%source(b)%sourceMass(:,i,j,k), res, datum )
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

  end subroutine IGLOOxMOSE


  !> Inject ICE's reaction on the gas into MOSE's residual: what each family's source hands the gas (ICE_Cell) and the
  !  removals ICE stores (Removal_Cell: a steady run, the rate of the store; a time-accurate one, this MOSE
  !  iteration's share of it as a rate over its step, which each stage and both Strang halves integrate to
  !  rate*dtlocal), under the consumed-species floor guard, times the cell volume.
  subroutine ICExMOSE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nprim
    use MOSE_Config_Types_m,    only: obj_chemistry
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(ICE_domain_type),  intent(in)    :: ICE_domain
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
            if (handRates) then
              if (allocated(goneRate)) then
                gr = goneRate(b)%a(:,g,i,j,k)
                call Removal_Cell ( imat, gr, r, datum )
              endif
            elseif (goneShare > 0._R8) then
              dtIter = MOSE_domain%blk(b)%dtlocal(i,j,k)
              if (dtIter > 0._R8) then
                gr = ICE_domain%blk(b)%cond_phase(g)%gone(:,i,j,k)*goneShare/dtIter
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

  end subroutine ICExMOSE


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
  !  Cap_Rates) with its datum offset, less the smoke born in ICE. datum: the offsets' part of what it adds to the
  !  energy. Thread-safe.
  subroutine IGLOO_Mass_To_Gas ( col, r, datum )
    use MOSE_Global_m,          only: nsc, nu, nw, np
    use IGLOO_variables,        only: nm
    use hydraMI2_Stoichiometry, only: productsIGLOO, Cap_Rates
    implicit none
    real(R8), intent(in)    :: col(:)
    real(R8), intent(inout) :: r(:)
    real(R8), intent(out)   :: datum
    integer  :: imat, isp
    real(R8) :: dm, cr(4), srate(6,2)

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
    do imat = 1, nm
      call IGLOO_Smoke_Rates( col, imat, srate )
      r(nu:nw) = r(nu:nw) - srate(2:4,1) - srate(2:4,2)
      r(np)    = r(np) - srate(5,1) - srate(5,2) - smokeIG(imat)%dh(1)*srate(1,1) - smokeIG(imat)%dh(2)*srate(1,2)
      datum    = datum - smokeIG(imat)%dh(1)*srate(1,1) - smokeIG(imat)%dh(2)*srate(1,2)
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


  !> What the removals ICE stores times share (1: what they still hold), or (rate) the rates in force (goneRate),
  !  hand the gas: Removal_Cell times the cell volume.
  function Removals_Handed ( MOSE_domain, ICE_domain, share, rate ) result( h )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, np, nprim
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(ICE_domain_type),  intent(in) :: ICE_domain
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
            gr = goneRate(b)%a(:,g,i,j,k)*vol
          else
            gr = ICE_domain%blk(b)%cond_phase(g)%gone(:,i,j,k)*share*vol
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

  !> At a MOSE block's start: how the removals ICE stores reach the gas over its nIter iterations (ICExMOSE reads
  !  them). Rates: where ICE stepped since the last block, the store over the time it was made in becomes the rate
  !  in force (goneRate), the store and its time emptied; rm holds the rates' totals. Amounts, once ICE has produced
  !  sources (ready): each iteration hands 1/nIter of the store.
  subroutine Open_Removals ( MOSE_domain, ICE_domain, ready, nIter, rm )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    logical,                intent(in)    :: ready
    integer,                intent(in)    :: nIter
    type(removals_t),       intent(out)   :: rm
    integer  :: b, g, i, j, k
    real(R8) :: tg

    if (.not. handRates) then
      if (ready .and. nIter > 0) goneShare = 1._R8 / real(nIter, R8)
      return
    endif

    rm%asRate = .true.
    if (.not. productsICE%active) return
    call Gone_Rate_Setup( ICE_domain )
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        if (.not. Converts(mat_of(g))) cycle
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          tg = goneRate(b)%t(g,i,j,k)
          if (.not. tg > 0._R8) cycle
          goneRate(b)%a(:,g,i,j,k) = ICE_domain%blk(b)%cond_phase(g)%gone(:,i,j,k)/tg
          ICE_domain%blk(b)%cond_phase(g)%gone(:,i,j,k) = 0._R8
          goneRate(b)%t(g,i,j,k) = 0._R8
        enddo; enddo; enddo
      enddo
    enddo
    rm%h = Removals_Handed( MOSE_domain, ICE_domain, 1._R8, rate=.true. )

  end subroutine Open_Removals


  !> After a MOSE block that made nRun of the iterations Open_Removals shared the store over (amounts): what they
  !  handed (rm%h), taken out of the store; the share closed.
  subroutine Close_Removals ( MOSE_domain, ICE_domain, nRun, rm )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use ICE_Global_m,          only: ngroups, mat_of
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    integer,                intent(in)    :: nRun
    type(removals_t),       intent(inout) :: rm
    integer  :: b, g
    real(R8) :: frac

    if (.not. goneShare > 0._R8) return
    frac      = real(nRun, R8)*goneShare
    goneShare = 0._R8
    rm%shared = .true.
    rm%h      = Removals_Handed( MOSE_domain, ICE_domain, frac )
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        if (Converts(mat_of(g))) ICE_domain%blk(b)%cond_phase(g)%gone = ICE_domain%blk(b)%cond_phase(g)%gone*(1._R8 - frac)
      enddo
    enddo

  end subroutine Close_Removals


  !> The arrays of goneRate, allocated at the first call and zero.
  subroutine Gone_Rate_Setup ( ICE_domain )
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: ngroups
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer :: b, n(3)

    if (allocated(goneRate)) return
    allocate(goneRate(ICE_domain%nb))
    do b = 1, ICE_domain%nb
      n = ICE_domain%blk(b)%dim(1:3)
      allocate(goneRate(b)%a(5, ngroups, n(1), n(2), n(3)), goneRate(b)%t(ngroups, n(1), n(2), n(3)), source=0._R8)
    enddo

  end subroutine Gone_Rate_Setup


  ! ---------------------------------------------------------------------------------------------
  ! ICE's clocks, after each of its steps and after its block. Outside any parallel region.
  ! ---------------------------------------------------------------------------------------------

  !> After an ICE iteration. A steady run (handRates): each converting family's dt added, cell by cell, to the time
  !  its store of removals is being made over. A time-accurate ICE: each family's dt*source, dt*heat, dt*birth (a
  !  family that makes a smoke) and dt added to the block's sums.
  subroutine ICE_Step_Clocks ( ICE_domain )
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use ICE_Config_Types_m,     only: ICE_time_scheme => obj_time_scheme
    use hydraMI2_Stoichiometry, only: productsICE
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer :: b, g, i, j, k

    if (handRates .and. productsICE%active) then
      call Gone_Rate_Setup( ICE_domain )
      do b = 1, ICE_domain%nb
        do g = 1, ngroups
          if (Converts(mat_of(g))) goneRate(b)%t(g,:,:,:) = goneRate(b)%t(g,:,:,:) + ICE_domain%blk(b)%cond_phase(g)%dt
        enddo
      enddo
    endif

    if (.not. ICE_time_scheme%time_accurate) return
    if (.not. allocated(srcSum)) allocate(srcSum(ICE_domain%nb, ngroups))
    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        associate (cp => ICE_domain%blk(b)%cond_phase(g), sm => srcSum(b,g))
          if (.not. allocated(sm%s)) then
            allocate(sm%s, mold=cp%source); allocate(sm%h, mold=cp%heat); allocate(sm%t, mold=cp%heat)
            allocate(sm%b, mold=cp%birth)
            sm%s = 0._R8; sm%h = 0._R8; sm%t = 0._R8; sm%b = 0._R8
          endif
          do k = 1, ICE_domain%blk(b)%dim(3)
          do j = 1, ICE_domain%blk(b)%dim(2)
          do i = 1, ICE_domain%blk(b)%dim(1)
            sm%s(:,i,j,k) = sm%s(:,i,j,k) + cp%dt(i,j,k)*cp%source(:,i,j,k)
            sm%h(i,j,k)   = sm%h(i,j,k)   + cp%dt(i,j,k)*cp%heat(i,j,k)
            sm%t(i,j,k)   = sm%t(i,j,k)   + cp%dt(i,j,k)
            sm%b(:,i,j,k) = sm%b(:,i,j,k) + cp%dt(i,j,k)*cp%birth(:,i,j,k)
          enddo; enddo; enddo
        end associate
      enddo
    enddo
    srcSteps = srcSteps + 1

  end subroutine ICE_Step_Clocks


  !> After a time-accurate ICE block: where it made more than one step, each family's source, heat and birth become
  !  the block's rates, the sums of ICE_Step_Clocks divided by the time; the sums cleared.
  subroutine ICE_Block_Mean ( ICE_domain )
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: ngroups
    implicit none
    type(ICE_domain_type), intent(inout) :: ICE_domain
    integer :: b, g, i, j, k

    if (.not. allocated(srcSum)) return
    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        associate (cp => ICE_domain%blk(b)%cond_phase(g), sm => srcSum(b,g))
          if (srcSteps > 1) then
            do k = 1, ICE_domain%blk(b)%dim(3)
            do j = 1, ICE_domain%blk(b)%dim(2)
            do i = 1, ICE_domain%blk(b)%dim(1)
              if (sm%t(i,j,k) > 0._R8) then
                cp%source(:,i,j,k) = sm%s(:,i,j,k)/sm%t(i,j,k)
                cp%heat(i,j,k)     = sm%h(i,j,k)/sm%t(i,j,k)
                cp%birth(:,i,j,k)  = sm%b(:,i,j,k)/sm%t(i,j,k)
              endif
            enddo; enddo; enddo
          endif
          sm%s = 0._R8; sm%h = 0._R8; sm%t = 0._R8; sm%b = 0._R8
        end associate
      enddo
    enddo
    srcSteps = 0

  end subroutine ICE_Block_Mean


  ! ---------------------------------------------------------------------------------------------
  ! IGLOO's smoke, from the gas's share to ICE
  ! ---------------------------------------------------------------------------------------------

  !> The smoke an IGLOO material's burning makes in a cell, as rates [per s, the cell's]: (:,1) the smoke of its
  !  combustion record, (:,2) the smoke its cap oxide handed at burnout becomes; mass, momentum (3), energy on ICE's
  !  datum with the kinetic part, number (each drop at d-smoke); born with the parcels' velocity and temperature
  !  weighted by what they burnt there (IGLOO's columns). col: the cell's column of IGLOO's mass source. Thread-safe.
  subroutine IGLOO_Smoke_Rates ( col, imat, rate, nostate )
    use IGLOO_variables,        only: nm, smokeRow
    use ICE_Config_Types_m,     only: ICE_condensed => obj_condensed
    use ICE_Global_m,           only: mat_of
    use ICE_Lib_Properties,     only: mat_rho, mat_e, mat_cp_const
    use hydraMI2_Stoichiometry, only: productsIGLOO, Cap_Rates, T_datum_ref
    implicit none
    real(R8), intent(in)  :: col(:)
    integer,  intent(in)  :: imat
    real(R8), intent(out) :: rate(6,2)
    logical,  intent(out), optional :: nostate      !> a smoke born where no parcel's columns gave its state
    real(R8), parameter :: pi = 3.14159265358979323846_R8
    real(R8) :: md(2), w, u(3), T, es, cr(4)
    integer  :: n, ms

    rate = 0._R8
    if (present(nostate)) nostate = .false.
    if (all(smokeIG(imat)%fam == 0)) return
    if (productsIGLOO%capped(imat)) then
      cr = Cap_Rates( imat, col(imat), col(nm+2*imat-1), col(nm+2*imat) )
      md = [smokeIG(imat)%nu(1)*cr(2), smokeIG(imat)%nu(2)*cr(4)]
    else
      md = [smokeIG(imat)%nu(1)*col(imat), 0._R8]
    endif
    if (all(md == 0._R8)) return
    w = col(smokeRow(imat,1))
    if (w > 0._R8) then
      u = col(smokeRow(imat,2):smokeRow(imat,4))/w
      T = col(smokeRow(imat,5))/w
    else
      u = 0._R8
      T = T_datum_ref
      if (present(nostate)) nostate = .true.
    endif
    do n = 1, 2
      if (md(n) == 0._R8 .or. smokeIG(imat)%fam(n) == 0) cycle
      ms = mat_of(smokeIG(imat)%fam(n))
      if (ICE_condensed(ms)%cs_varies) then
        es = mat_e(ICE_condensed(ms), T)
      else
        es = mat_cp_const(ICE_condensed(ms))*T
      endif
      rate(1,n)   = md(n)
      rate(2:4,n) = md(n)*u
      rate(5,n)   = md(n)*(es + 0.5_R8*dot_product(u, u))
      rate(6,n)   = md(n)/(pi/6._R8*mat_rho(ICE_condensed(ms), T)*smokeIG(imat)%dSmoke**3)
    enddo

  end subroutine IGLOO_Smoke_Rates


  !> After a MOSE iteration: the smoke rates of IGLOO's materials times the iteration's dtlocal, added cell by cell to
  !  the amounts of the block (the gas received the same rates less, over the same step). Outside any parallel region.
  subroutine IGLOO_Smoke_Accumulate ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    use IGLOO_variables,       only: nm
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    integer  :: b, i, j, k, imat, n, f
    real(R8) :: srate(6,2), dt
    logical  :: nostate

    call Smoke_Acc_Setup( MOSE_domain )
    do b = 1, size(MOSE_domain%blk)
      if (.not. allocated(IGLOO%source(b)%sourceMass)) cycle
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        dt = MOSE_domain%blk(b)%dtlocal(i,j,k)
        do imat = 1, nm
          call IGLOO_Smoke_Rates( IGLOO%source(b)%sourceMass(:,i,j,k), imat, srate, nostate )
          if (nostate) nSmokeNoState = nSmokeNoState + 1
          do n = 1, 2
            if (srate(1,n) == 0._R8) cycle
            f = smokeIG(imat)%fam(n)
            smokeAcc(b)%a(:,f,i,j,k) = smokeAcc(b)%a(:,f,i,j,k) + srate(:,n)*dt
          enddo
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine IGLOO_Smoke_Accumulate


  !> The smoke amounts of the MOSE blocks and each smoke family's datum offset, allocated at the first call, zero.
  subroutine Smoke_Acc_Setup ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Global_m,          only: ngroups
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer :: b, imat, n

    if (allocated(smokeAcc)) return
    allocate(smokeAcc(size(MOSE_domain%blk)), dhSmokeFam(ngroups))
    dhSmokeFam = 0._R8
    do imat = 1, size(smokeIG)
      do n = 1, 2
        if (smokeIG(imat)%fam(n) > 0) dhSmokeFam(smokeIG(imat)%fam(n)) = smokeIG(imat)%dh(n)
      enddo
    enddo
    do b = 1, size(MOSE_domain%blk)
      allocate(smokeAcc(b)%a(6, ngroups, MOSE_domain%blk(b)%dim(1), MOSE_domain%blk(b)%dim(2), &
                             MOSE_domain%blk(b)%dim(3)), source=0._R8)
    enddo

  end subroutine Smoke_Acc_Setup


  !> After a MOSE block: the smoke it made (IGLOO_Smoke_Accumulate), which the gas did not get, handed to ICE in the
  !  form of handRates. Amounts: added to ICE's state cell by cell (per volume). Rates: the rates of IGLOO's sources as
  !  they stand (srcOn; else 0), per volume, set as the source of each smoke family (cond_phase%ext). sm: the totals,
  !  energies on the gas's formation scale; the amounts emptied. Outside any parallel region.
  subroutine IGLOO_Smoke_Deliver ( MOSE_domain, ICE_domain, IGLOO, srcOn, sm )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    use IGLOO_variables,       only: nm
    use ICE_Global_m,          only: ngroups, mat_of, nbase, ncond
    use ICE_Config_Types_m,    only: ICE_condensed => obj_condensed
    use ICE_Lib_Model,         only: prim_2_cons, cons_2_prim, assign_prim_2_cons, assign_cons_2_prim
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: srcOn
    type(smoke_t),          intent(out)   :: sm
    integer  :: b, g, i, j, k, ns, imat, n
    real(R8) :: vol, ac(6), srate(6,2), sr(6)
    real(R8), allocatable :: c(:)

    sm%nNoState   = nSmokeNoState
    nSmokeNoState = 0
    call Smoke_Acc_Setup( MOSE_domain )
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        if (.not. any(smokeIG(:)%fam(1) == g .or. smokeIG(:)%fam(2) == g)) cycle
        ! the conversions of family g's closure (ICE points them at the family it last stepped)
        call assign_prim_2_cons(g)
        call assign_cons_2_prim(g)
        ns = nbase(g)
        if (allocated(c)) deallocate(c)
        allocate(c(ncond(g)))
        if (handRates) then
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
          if (handRates .and. srcOn) then
            if (allocated(IGLOO%source(b)%sourceMass)) then
              sr = 0._R8
              do imat = 1, nm
                call IGLOO_Smoke_Rates( IGLOO%source(b)%sourceMass(:,i,j,k), imat, srate )
                do n = 1, 2
                  if (smokeIG(imat)%fam(n) == g) sr = sr + srate(:,n)
                enddo
              enddo
              c = 0._R8
              c(1:4)  = sr(1:4)/vol
              c(ns-1) = sr(5)/vol
              c(ns)   = sr(6)/vol
              ICE_domain%blk(b)%cond_phase(g)%ext(:,i,j,k) = c
              sm%rate(1:4) = sm%rate(1:4) + sr(1:4)
              sm%rate(5)   = sm%rate(5) + sr(5) + dhSmokeFam(g)*sr(1)
              sm%rate(6)   = sm%rate(6) + sr(6)
            endif
          endif
          ! the amounts of the block: out of the gas's share, and into ICE's state when they go as amounts
          ac = smokeAcc(b)%a(:,g,i,j,k)
          if (all(ac == 0._R8)) cycle
          if (.not. handRates) then
            c = prim_2_cons(ICE_domain%blk(b)%cond_phase(g)%prim(:,i,j,k), ICE_condensed(mat_of(g)))
            c(1:4)  = c(1:4)  + ac(1:4)/vol
            c(ns-1) = c(ns-1) + ac(5)/vol
            c(ns)   = c(ns)   + ac(6)/vol
            ICE_domain%blk(b)%cond_phase(g)%prim(:,i,j,k) = cons_2_prim(c, ICE_condensed(mat_of(g)))
          endif
          sm%amount(1:4) = sm%amount(1:4) + ac(1:4)
          sm%amount(5)   = sm%amount(5) + ac(5) + dhSmokeFam(g)*ac(1)
          smokeAcc(b)%a(:,g,i,j,k) = 0._R8
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine IGLOO_Smoke_Deliver

end module hydraMI2_Mod_Coupling
