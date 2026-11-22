!> What the coupling amounts to, for the log. Nothing here acts on the solution.
!    setup         : what got built, the overset-mesh warning, the gas inventory the books start from (Report_Start)
!    every block   : the frozen-source notices; the books of what the block handed the gas and, under [MI2-Coupling]
!                    report, its [MI2-report] lines (Report_Block_Start ... Report_Block_End)
!    end of run    : the floor-guard count and the closed-domain mass and energy balances (Report_End)
!  What a solver hands the gas is summed by hydraMI2_Mod_Coupling (IGLOO_Handed, ICE_Handed, Removals_Handed),
!  from the same cell kernels the injection applies.
module hydraMI2_Report
  use iso_fortran_env,       only: R8 => real64
  use hydraMI2_Mod_Coupling, only: handed_t

  implicit none
  private

  public :: Report_Start, Report_Frozen_Sources, Report_Values, Report_Inlet
  public :: Report_Block_Start, Report_Block_End, Report_End

  !> The books: the mass of each gas species [kg] (allocated by the first block or delivery that carries a mass
  !  coupling), the energy [J] and its datum-offset part, handed to the gas over the run.
  real(R8), allocatable :: massInjected(:)
  real(R8) :: energyInjected = 0._R8, energyDatumPart = 0._R8
  logical  :: energyAudited  = .false.
  !> Blocks whose IGLOO sources were Eulerian (mass-inert), left out of the books.
  integer  :: nAuditSkipped = 0

  !> The gas inventory at the start of the run: mass per species [kg], total energy [J], mean T [K].
  real(R8), allocatable :: massGas0(:)
  real(R8) :: enGas0 = 0._R8, Tgas0 = 0._R8

  !> The current MOSE block: its domain clock at the start and, under report, the gas momentum and energy there.
  real(R8) :: tBlock0 = 0._R8, momBlock0(3) = 0._R8, enBlock0 = 0._R8

  !> The cumulative MOSE iteration of the last frozen-source notice.
  integer :: frozenPrintIter = -1

contains

  ! ---------------------------------------------------------------------------------------------
  ! Setup
  ! ---------------------------------------------------------------------------------------------

  !> After every solver's setup: what got built, a warning on an overset mesh (gas faces of kind 102, whose overlap
  !  every sum over the domain counts once per block that covers it), and the gas inventory the balances start from.
  subroutine Report_Start ( MOSE_simulation, ICE_simulation, IGLOO )
    use MOSE_Advanced_Types_m,   only: MOSE_simulation_type
    use MOSE_Global_m,           only: nsc
    use ICE_Advanced_Types_m,    only: ICE_simulation_type
    use IGLOO_module,            only: obj_IGLOO
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    type(MOSE_simulation_type), intent(in) :: MOSE_simulation
    type(ICE_simulation_type),  intent(in) :: ICE_simulation
    type(obj_IGLOO),            intent(in) :: IGLOO
    real(R8) :: mgas, cvbar
    integer  :: i, n
    logical, allocatable :: has(:)

    associate ( gas => MOSE_simulation%domain(1) )

    allocate(massGas0(nsc))
    call Gas_Species_Mass    ( gas, massGas0 )
    call Gas_Total_Energy    ( gas, enGas0 )
    call Gas_Mean_Temperature( gas, Tgas0, mgas, cvbar )

    write(*,*)
    write(*,'(A)') ' HYDRA-MI2 solvers'
    write(*,'(A,I0)') ' -  MOSE  blocks              --> ', gas%nb
    if (hydra_MI2_obj%run_ICE)   write(*,'(A,I0)') ' -  ICE   blocks              --> ', ICE_simulation%domain(1)%nb
    if (hydra_MI2_obj%run_IGLOO) write(*,'(A,I0)') ' -  IGLOO geometry blocks     --> ', size(IGLOO%geoblock)
    if (hydra_MI2_obj%run_IGLOO) write(*,'(A,I0)') ' -  IGLOO materials           --> ', size(IGLOO%material)
    write(*,*)

    if ( .not. allocated(gas%bc) ) return
    allocate(has(gas%nb)); has = .false.
    n = 0
    do i = 1, size(gas%bc)
      if ( gas%bc(i)%type /= 102 ) cycle
      n = n + 1
      has(gas%bc(i)%b) = .true.
    enddo
    if ( n > 0 ) write(*,'(A,I0,A,I0,A)') ' [WARNING] MI2: an overset mesh: ', n, ' chimera (102) boundary cells on ', &
      count(has), ' block(s). The coupler''s sums over the domain (the report and the balances) count an overlap '// &
      'once per block that covers it, uncorrected'

    end associate

  end subroutine Report_Start


  ! ---------------------------------------------------------------------------------------------
  ! Every coupling step
  ! ---------------------------------------------------------------------------------------------

  !> Announce sources frozen this cycle (X-iter(n_cycle) = 0, sources already produced): at the first coupling step
  !  of the cycle, then at [MOSE-IO] shell-diter cadence.
  subroutine Report_Frozen_Sources ( IGLOO_src, ICE_src, n_cycle, n_couple )
    use MOSE_Config_Types_m,     only: MOSE_sim => obj_sim_param, MOSE_io => obj_io
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydraMI2_Mod_Coupling,   only: source_state_t
    implicit none
    type(source_state_t), intent(in) :: IGLOO_src, ICE_src
    integer,              intent(in) :: n_cycle, n_couple
    character(len=*), parameter :: FMT = '(A,I0,A,I0,A,I0,A,I0,A,I0,A)'
    integer :: iterNow, cadence
    logical :: iglooFrozen, iceFrozen

    ! X_iter has size 1 when its solver is off, and its sources are then never ready
    iglooFrozen = .false.
    if ( IGLOO_src%ready ) iglooFrozen = hydra_MI2_obj%IGLOO_iter(n_cycle) == 0
    iceFrozen = .false.
    if ( ICE_src%ready )   iceFrozen   = hydra_MI2_obj%ICE_iter(n_cycle) == 0
    if ( .not. (iglooFrozen .or. iceFrozen) ) return

    iterNow = MOSE_sim%iter_general
    cadence = max(MOSE_io%shell_diter, 1)
    if ( n_couple /= 1 .and. frozenPrintIter >= 0 .and. iterNow/cadence <= frozenPrintIter/cadence ) return
    frozenPrintIter = iterNow

    if ( iglooFrozen ) write(*,FMT) &
      ' [MI2] IGLOO source FROZEN this cycle (IGLOO-iter = 0): the gas keeps the sources produced in cycle ', &
      IGLOO_src%cycle, ', coupling step ', IGLOO_src%step, '  [cycle ', n_cycle, ', step ', n_couple, &
      ', MOSE iter ', iterNow, ']'
    if ( iceFrozen .and. ICE_src%cycle > 0 ) write(*,FMT) &
      ' [MI2] ICE source FROZEN this cycle (ICE-iter = 0): the gas keeps the sources produced in cycle ', &
      ICE_src%cycle, ', coupling step ', ICE_src%step, '  [cycle ', n_cycle, ', step ', n_couple, &
      ', MOSE iter ', iterNow, ']'
    if ( iceFrozen .and. ICE_src%cycle == 0 ) write(*,'(A,I0,A,I0,A,I0,A)') &
      ' [MI2] ICE source FROZEN this cycle (ICE-iter = 0): the gas keeps the sources read at the restart  [cycle ', &
      n_cycle, ', step ', n_couple, ', MOSE iter ', iterNow, ']'

  end subroutine Report_Frozen_Sources


  !> One [MI2-report] line: the tag, then the values.
  subroutine Report_Values ( tag, v )
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    character(len=*), intent(in) :: tag
    real(R8),         intent(in) :: v(:)

    if ( .not. hydra_MI2_obj%report ) return
    write(*,'(A,*(ES24.16))') ' [MI2-report] '//tag, v

  end subroutine Report_Values


  !> The [MI2-report] line of what a solver handed the gas: the mass into the gas species, momentum, energy.
  subroutine Report_Handed ( tag, h )
    implicit none
    character(len=*), intent(in) :: tag
    type(handed_t),   intent(in) :: h

    call Report_Values( tag, [Mass_Of(h), h%mom, h%energy] )

  end subroutine Report_Handed


  !> The [MI2-report] line of the inflow mass flow MOSExIGLOO handed IGLOO, and the face records carrying it.
  subroutine Report_Inlet ( mdot, n )
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    real(R8), intent(in) :: mdot
    integer,  intent(in) :: n

    if ( .not. hydra_MI2_obj%report ) return
    write(*,'(A,ES24.16,I10)') ' [MI2-report] igloo-inlet-mdot', mdot, n

  end subroutine Report_Inlet


  !> Before a MOSE block: its domain clock, and under report the gas momentum and energy.
  subroutine Report_Block_Start ( MOSE_domain )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain

    tBlock0 = MOSE_domain%time
    if ( .not. hydra_MI2_obj%report ) return
    call Gas_Momentum    ( MOSE_domain, momBlock0 )
    call Gas_Total_Energy( MOSE_domain, enBlock0 )

  end subroutine Report_Block_Start


  !> After a MOSE block: what it handed the gas, into the books and, under report, onto the log. ICE's stored
  !  removals (rm: the rates in force, or the amounts delivered) and IGLOO's smoke given to ICE (sm); then the
  !  sources in force, before the solves of this coupling step replace them, over the block's clock; and what ICE's
  !  removals still hold. ICE_simulation is read only when ICE's sources are ready.
  subroutine Report_Block_End ( MOSE_domain, ICE_simulation, IGLOO, IGLOO_src, ICE_src, n_cycle, n_couple, rm, sm )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use ICE_Advanced_Types_m,    only: ICE_simulation_type
    use ICE_Global_m,            only: ngroups
    use IGLOO_module,            only: obj_IGLOO
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydraMI2_Stoichiometry,  only: productsICE
    use hydraMI2_Mod_Coupling,   only: source_state_t, removals_t, smoke_t, handRates, smokeIGOn, &
                                       IGLOO_Handed, ICE_Handed, Removals_Handed
    implicit none
    type(MOSE_domain_type),    intent(in) :: MOSE_domain
    type(ICE_simulation_type), intent(in) :: ICE_simulation
    type(obj_IGLOO),           intent(in) :: IGLOO
    type(source_state_t),      intent(in) :: IGLOO_src, ICE_src
    integer,                   intent(in) :: n_cycle, n_couple
    type(removals_t),          intent(in) :: rm
    type(smoke_t),             intent(in) :: sm
    type(handed_t) :: hIGLOO, hICE
    real(R8) :: dt, mom1(3), en1, tau(2)
    integer  :: b, g, n(3)
    logical  :: iglooSrc

    ! ICE's removals and IGLOO's smoke, between the solvers' steps
    if ( rm%asRate .and. allocated(rm%h%mass) ) call Report_Handed( 'ice-removed-rate', rm%h )
    if ( rm%shared ) then
      ! the removals carry mass only when an ICE material converts (rm%h%mass then allocated)
      if ( allocated(rm%h%mass) ) call Book( rm%h, 1._R8 )
      call Report_Handed( 'ice-removed-delivered', rm%h )
    endif
    if ( smokeIGOn ) then
      if ( handRates ) then
        call Report_Values( 'igloo-smoke-withheld', sm%amount )
        call Report_Values( 'igloo-smoke-rate', sm%rate )
      else
        call Report_Values( 'igloo-smoke-delivered', sm%amount )
      endif
      if ( sm%nNoState > 0 ) write(*,'(A,I0,A)') ' [WARNING] MI2: in ', sm%nNoState, ' cells of this block''s '// &
        'iterations a smoke was born with no parcel''s state in IGLOO''s columns: at rest, at T_datum_ref'
    endif

    ! the sources in force over the block (IGLOO's Eulerian branch is mass-inert: left out of the books)
    dt       = MOSE_domain%time - tBlock0
    iglooSrc = IGLOO_src%ready .and. .not. IGLOO_src%eulerian
    if ( iglooSrc ) then
      hIGLOO = IGLOO_Handed( MOSE_domain, IGLOO )
      call Book( hIGLOO, dt )
    elseif ( IGLOO_src%ready ) then
      nAuditSkipped = nAuditSkipped + 1
    endif
    if ( ICE_src%ready ) then
      hICE = ICE_Handed( MOSE_domain, ICE_simulation%domain(1) )
      if ( productsICE%active ) call Book( hICE, dt )
    endif

    if ( .not. hydra_MI2_obj%report ) return
    call Gas_Momentum    ( MOSE_domain, mom1 )
    call Gas_Total_Energy( MOSE_domain, en1 )
    write(*,'(A,I0,A,I0,2(A,ES24.16))') ' [MI2-report] cycle ', n_cycle, ' step ', n_couple, &
                                         ' t0 ', tBlock0, ' t1 ', MOSE_domain%time
    call Report_Values( 'gas-momentum', [momBlock0, mom1] )
    call Report_Values( 'gas-energy',   [enBlock0, en1] )
    if ( ICE_src%ready ) then
      ! the range of ICE's relaxation time over the interior cells and the families
      tau = [huge(1._R8), -huge(1._R8)]
      do b = 1, size(MOSE_domain%blk)
        n = MOSE_domain%blk(b)%dim(1:3)
        do g = 1, ngroups
          associate ( t => ICE_simulation%domain(1)%blk(b)%cond_phase(g)%tau(1:n(1),1:n(2),1:n(3)) )
          tau = [min(tau(1), minval(t)), max(tau(2), maxval(t))]
          end associate
        enddo
      enddo
      call Report_Values( 'ice-to-gas', [hICE%mom, hICE%energy, tau] )
      if ( productsICE%active ) then
        call Report_Values( 'ice-to-gas-mass', [Mass_Of(hICE)] )
        call Report_Handed( 'ice-removed-stored', Removals_Handed( MOSE_domain, ICE_simulation%domain(1), 1._R8 ) )
      endif
    endif
    if ( iglooSrc ) call Report_Values( 'igloo-to-gas', [hIGLOO%mom, hIGLOO%energy, Mass_Of(hIGLOO)] )

  end subroutine Report_Block_End


  !> Add what h hands the gas, over dt (1 for an amount), to the books. A mass coupling opens the mass books.
  subroutine Book ( h, dt )
    implicit none
    type(handed_t), intent(in) :: h
    real(R8),       intent(in) :: dt

    if ( allocated(h%mass) .and. .not. allocated(massInjected) ) allocate(massInjected(size(h%mass)), source=0._R8)
    if ( .not. dt > 0._R8 ) return
    if ( allocated(h%mass) ) massInjected = massInjected + h%mass*dt
    energyInjected  = energyInjected  + h%energy*dt
    energyDatumPart = energyDatumPart + h%datum*dt
    energyAudited   = .true.

  end subroutine Book


  !> The mass h hands the gas species, all of them.
  pure real(R8) function Mass_Of ( h )
    implicit none
    type(handed_t), intent(in) :: h

    Mass_Of = 0._R8
    if ( allocated(h%mass) ) Mass_Of = sum(h%mass)

  end function Mass_Of


  ! ---------------------------------------------------------------------------------------------
  ! End of the run: the floor-guard count and the books closed against the gas inventory. The balances are
  ! meaningful for a time-accurate run on a domain without boundary flux.
  ! ---------------------------------------------------------------------------------------------

  subroutine Report_End ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use hydraMI2_Mod_Coupling, only: nMassLimited
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain

    if ( nMassLimited > 0 ) then
      write(*,*)
      write(*,'(A,I0,A)') ' [MI2] WARNING: the consumed-species floor guard clamped ', nMassLimited, &
        ' (cell,species) extractions: a gas species ran out where a reaction consumed it, so the reaction was '// &
        'MASS-LIMITED, not the stoichiometry asked for'
    endif

    write(*,*)
    if ( allocated(massInjected) ) call Report_Mass_Balance  ( MOSE_domain )
    if ( energyAudited )           call Report_Energy_Balance( MOSE_domain )

  end subroutine Report_End


  !> Injected mass per species and the closed-domain balance against the change in gas mass.
  subroutine Report_Mass_Balance ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use hydra_Config_Types_m,  only: hydra_obj_sim_param
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    real(R8) :: massGas1(size(massGas0)), dgas, den, floor_
    integer  :: is

    write(*,*)
    write(*,'(A)') ' [MI2] mass injected into the gas over the run [kg]'
    if ( nAuditSkipped > 0 ) write(*,'(A,I0,A)') '   [WARNING] ', nAuditSkipped, &
      ' coupling step(s) used the EULERIAN branch, which is mass-inert, and are excluded.'
    if ( .not. hydra_obj_sim_param%time_accurate ) &
      write(*,'(A)') '   (time-accurate is OFF -- these are pseudo-time integrals, not physical mass)'
    do is = 1, size(massInjected)
      if ( abs(massInjected(is)) > 0._R8 ) write(*,'(A,I0,A,ES14.6)') '   species ', is, ' : ', massInjected(is)
    enddo
    write(*,'(A,ES14.6)') '   TOTAL              : ', sum(massInjected)

    ! injected vs gas gained, both measured at the same instant
    call Gas_Species_Mass( MOSE_domain, massGas1 )
    write(*,*)
    write(*,'(A)') ' [MI2] closed-domain mass balance  (valid only if NO boundary flux)'
    write(*,'(A)') '   species        injected           gas gained         rel.err'
    ! ignore species whose change is round-off on the total inventory
    floor_ = 1.e-12_R8 * max(sum(abs(massGas0)), tiny(1._R8))
    do is = 1, size(massGas1)
      dgas = massGas1(is) - massGas0(is)
      den  = max(abs(massInjected(is)), abs(dgas))
      if ( den <= floor_ ) then
        write(*,'(A,I0,A)') '   ', is, '   (untouched; change is round-off on the inventory)'
      else
        write(*,'(A,I0,A,ES16.6,A,ES16.6,A,ES11.3)') '   ', is, '   ', massInjected(is), &
              '   ', dgas, '   ', abs(dgas-massInjected(is))/den
      endif
    enddo
    write(*,*)

  end subroutine Report_Mass_Balance


  !> Energy bookkeeping (datum-blind) and the enthalpy-datum diagnostic.
  subroutine Report_Energy_Balance ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use hydra_Config_Types_m,  only: hydra_obj_sim_param
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    real(R8) :: enGas1, Tgas1, mGas1, cvGas1, den

    call Gas_Total_Energy    ( MOSE_domain, enGas1 )
    call Gas_Mean_Temperature( MOSE_domain, Tgas1, mGas1, cvGas1 )
    write(*,*)
    write(*,'(A)') ' [MI2] closed-domain ENERGY bookkeeping  (valid only if NO boundary flux)'
    if ( .not. hydra_obj_sim_param%time_accurate ) &
      write(*,'(A)') '   (time-accurate is OFF -- pseudo-time integral, not physical energy)'
    if ( nAuditSkipped > 0 ) write(*,'(A,I0,A)') '   [WARNING] ', nAuditSkipped, ' coupling step(s) used the '// &
      'EULERIAN branch and are NOT in these totals: the balance below is incomplete'
    write(*,'(A,ES16.8)') '   handed to gas  (sourceEn + datum) [J] : ', energyInjected
    write(*,'(A,ES16.8)') '     of which the datum correction   [J] : ', energyDatumPart
    write(*,'(A,ES16.8)') '   gas total energy gained          [J] : ', enGas1 - enGas0
    den = max(abs(energyInjected), abs(enGas1-enGas0))
    if (den > 0._R8) write(*,'(A,ES11.3)') '   rel.err (bookkeeping only, datum-blind) : ', &
                                           abs((enGas1-enGas0) - energyInjected)/den

    write(*,*)
    write(*,'(A)') ' [MI2] enthalpy-datum diagnostic'
    write(*,'(A,F12.6,A,F12.6,A,ES11.3)') '   mass-weighted gas T:  ', Tgas0, ' -> ', Tgas1, '   delta = ', Tgas1-Tgas0
    if ( abs(energyDatumPart) > 0._R8 .and. mGas1*cvGas1 > 0._R8 ) then
      write(*,'(A,F12.4,A,ES12.4)') '   mass-weighted frozen cv [J/kg/K] : ', cvGas1, '   gas mass [kg] : ', mGas1
      write(*,'(A,ES16.8)') '   energy the gas would have received WITHOUT it [J] : ', energyInjected - energyDatumPart
      write(*,'(A,F12.6)')  '   => predicted delta-T without the correction  [K] : ', &
                            (Tgas1-Tgas0) - energyDatumPart/(mGas1*cvGas1)
      write(*,'(A)') '   (at a FIXED evaporated mass: a run with datum-correction = false evaporates another; the '// &
                     'invariant is the per-kg offset printed at setup)'
    endif
    write(*,*)

  end subroutine Report_Energy_Balance


  ! ---------------------------------------------------------------------------------------------
  ! Gas inventories, over the interior cells
  ! ---------------------------------------------------------------------------------------------

  !> Mass of each gas species in the domain [kg] = SUM_cells rho_s*vol.
  subroutine Gas_Species_Mass ( MOSE_domain, mass )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: mass(nsc)
    integer :: b, i, j, k, s

    mass = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        do s = 1, nsc
          mass(s) = mass(s) + MOSE_domain%blk(b)%P(s,i,j,k) * MOSE_domain%blk(b)%vol(i,j,k)
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine Gas_Species_Mass


  !> Gas momentum [kg m/s] = SUM_cells rho*u*vol.
  subroutine Gas_Momentum ( MOSE_domain, mom )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc, nu, nw
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: mom(3)
    integer :: b, i, j, k

    mom = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        mom = mom + sum(MOSE_domain%blk(b)%P(1:nsc,i,j,k)) * MOSE_domain%blk(b)%P(nu:nw,i,j,k) &
                  * MOSE_domain%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
    enddo

  end subroutine Gas_Momentum


  !> Total gas energy [J] = SUM_cells cons(np)*vol (kinetic energy included).
  subroutine Gas_Total_Energy ( MOSE_domain, energy )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: np
    use FLINT_Lib_Thermodynamic, only: prim2cons
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: energy
    integer  :: b, i, j, k
    real(R8) :: cons(np)

    energy = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        cons   = prim2cons( MOSE_domain%blk(b)%P(1:np,i,j,k) )
        energy = energy + cons(np) * MOSE_domain%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
    enddo

  end subroutine Gas_Total_Energy


  !> Mass-weighted mean gas temperature [K], gas mass [kg] and mass-weighted frozen cv [J/kg/K].
  subroutine Gas_Mean_Temperature ( MOSE_domain, Tbar, mgas, cvbar )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: nsc, np
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_tabT, cp_tab, Ri_tab
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: Tbar, mgas, cvbar
    integer  :: b, i, j, k, s
    real(R8) :: rho, R, T, dm, cvcell

    Tbar = 0._R8; mgas = 0._R8; cvbar = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        call co_rotot_Rtot ( MOSE_domain%blk(b)%P(1:nsc,i,j,k), rho, R )
        T    = EOS ( MOSE_domain%blk(b)%P(np,i,j,k), rho, R )
        dm   = rho * MOSE_domain%blk(b)%vol(i,j,k)
        cvcell = 0._R8
        do s = 1, nsc
          cvcell = cvcell + MOSE_domain%blk(b)%P(s,i,j,k) * ( f_tabT(T,s,cp_tab) - Ri_tab(s) )
        enddo
        mgas  = mgas  + dm
        Tbar  = Tbar  + dm*T
        cvbar = cvbar + cvcell * MOSE_domain%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
    enddo
    if (mgas > 0._R8) then
      Tbar  = Tbar  / mgas
      cvbar = cvbar / mgas
    endif

  end subroutine Gas_Mean_Temperature

end module hydraMI2_Report
