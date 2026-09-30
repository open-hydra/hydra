!>
!> hydra-MI2 -- MOSE (Eulerian gas) + ICE (Eulerian condensed) + IGLOO (Lagrangian) coupler.
!> OpenMP only: the root CMakeLists refuses USE_MPI for MI2.
!>
program hydra_MI2_program
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use hydra_Parameters_m
  use hydra_Config_Types_m
  use hydraMI2_Read_Ini,       only: Read_Inifile_MI2, Read_Inifile_MI2_Runtime
  use hydraMI2_Config_Types_m, only: hydra_MI2_obj
  use hydraMI2_Mod_Import,     only: Assign_Phases
  use hydraMI2_Stoichiometry,  only: Build_Species_Stoichiometry, Report_Stoichiometry, &
                                     has_stoichiometry, Build_Datum_Offsets, T_datum_ref
  use hydraMI2_Mod_Coupling,   only: MOSExIGLOO, Gas_Field_For_IGLOO, IGLOOxMOSE, &
                                     MOSExICE, ICExMOSE, &
                                     resField, Initialize_Residual_Store, Build_Family_Material_Map, &
                                     Accumulate_Injected_Mass, massInjected, Gas_Species_Mass, &
                                     nMassLimited, Accumulate_Injected_Energy, Gas_Total_Energy, &
                                     Gas_Mean_Temperature, energyInjected, energyDatumPart, &
                                     energyAudited, Check_ICE_Mass_Channel
  use hydraMI2_Build_Stamp,    only: Print_Build_Stamp
  use Lib_ORION_data,          only: orion_data

  use MOSE_Advanced_Types_m, only: MOSE_simulation_type, MOSE_domain_type
  use MOSE_Config_Types_m,   only: MOSE_obj_sim_param => obj_sim_param, MOSE_obj_io => obj_io
  use MOSE_Procedures_m,     only: MOSE_type

  use ICE_Advanced_Types_m,  only: ICE_simulation_type
  use ICE_Config_Types_m,    only: ICE_obj_sim_param => obj_sim_param
  use ICE_Procedures_m,      only: ICE_type

  use IGLOO_module,          only: obj_IGLOO

  implicit none

  type(MOSE_type)            :: MOSE
  type(MOSE_simulation_type) :: MOSE_simulation
  type(ICE_type)             :: ICE
  type(ICE_simulation_type)  :: ICE_simulation
  type(obj_IGLOO)            :: IGLOO
  ! Local
  integer :: nthreads
  integer :: n_cycle, n_couple, N_cycles
  integer :: TODO, i_MOSE, i_ICE
  logical :: IGLOO_first_sweep = .true.
  real(R8) :: t_block0 = 0._R8
  real(R8), allocatable :: massGas0(:), massGas1(:)
  real(R8) :: enGas0 = 0._R8, enGas1 = 0._R8
  real(R8) :: Tgas0  = 0._R8, Tgas1  = 0._R8, mGas0 = 0._R8, mGas1 = 0._R8
  real(R8) :: cvGas0 = 0._R8, cvGas1 = 0._R8
  type(orion_data) :: gas_for_igloo
  type(resField), allocatable :: oldResidual(:)
  !> Set once a solver has produced sources and never cleared: frozen sources keep acting on the gas.
  logical :: IGLOO_has_sources = .false.  ! true once IGLOO has produced sources MOSE may read
  logical :: ICE_has_sources   = .false.  ! ditto for ICE: cond_phase%source is only valid after a solve
  integer :: IGLOO_src_cycle = 0, IGLOO_src_step = 0   ! (cycle, coupling step) that produced IGLOO's sources
  integer :: ICE_src_cycle   = 0, ICE_src_step   = 0   ! same for ICE
  integer :: frozenPrintIter = -1                      ! cumulative MOSE iteration of the last frozen notice
  logical :: eulCoupling       = .false.   ! true when this cycle uses the Eulerian coupling mode
  !> Coupling mode in force when the sources MOSE is about to consume were produced.
  logical :: IGLOO_src_eulerian = .false.
  integer :: nAuditSkipped      = 0        ! coupling steps whose sources were Eulerian, so not audited

  call Print_Build_Stamp()

#if defined (_OPENMP)
  !$omp parallel
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  write(*,'(A)')    ' Parallel execution'
  write(*,'(A)')    ' OpenMP:'
  write(*,'(A,I4)') ' -  Number of threads --> ', nthreads
#else
  nthreads = 1
  write(*,'(A)')    ' Serial execution'
#endif
  MOSE_obj_sim_param%nthreads = nthreads
  ICE_obj_sim_param%nthreads  = nthreads

  ! Populate hydra_obj_* AND hydra_MI2_obj from input.ini
  call Read_Inifile_MI2 ()

  ! Report the parsed coupling schedule before any solver setup
  call Report_MI2_Config()

  ! Multigrid is not supported: IGLOO has no levels
  if (hydra_obj_multigrid%MGL > 1) then
    write(*,'(A,I0,A)') '[ERROR] hydra-MI2 does not support multigrid (HYDRA-Multigrid levels = ', &
                        hydra_obj_multigrid%MGL, '). IGLOO has no multigrid levels. Use levels = 1.'
    error stop 1
  endif

  ! Route each INPUT/*phase.txt to its owning solver (ICE and IGLOO each solve their own phase)
  call Assign_Phases()

  ! ---------------------------------------------------------------------------
  ! Solver setup
  ! ---------------------------------------------------------------------------
  ! ICE must know it is two-way coupled before setup (it then fills cond_phase%source)
  ICE_obj_sim_param%twcoupled = hydra_MI2_obj%run_ICE

  call MOSE%setup( MOSE_simulation )

  ! ICE and IGLOO ingest the gas field as ORION data; each is set up only when enabled
  if (hydra_MI2_obj%run_ICE)   call ICE%setup( ICE_simulation, MOSE_simulation%IOfield(1) )
  ! IGLOO%setup includes the reset_state of the first solve; every later solve needs its own
  if (hydra_MI2_obj%run_IGLOO) then
    call Gas_Field_For_IGLOO( MOSE_simulation%domain(1), MOSE_simulation%IOfield(1), gas_for_igloo )
    call IGLOO%setup( gas_for_igloo )
    call Build_Family_Material_Map( IGLOO )
    ! material -> gas-species stoichiometry for mass coupling (absent => off)
    block
      character(len=20), allocatable :: matnames(:)
      integer :: im
      allocate(matnames(size(IGLOO%material)))
      do im = 1, size(IGLOO%material)
        matnames(im) = IGLOO%material(im)%matName
      enddo
      call Build_Species_Stoichiometry( matnames, size(IGLOO%material) )
      call Report_Stoichiometry( matnames, size(IGLOO%material) )

      ! enthalpy-datum offsets; hIGLOOref is the enthalpy IGLOO credits for the material at T_ref on
      ! the datum of its table: cp*T_ref + hOff for constant cp, hTab interpolated at T_ref otherwise
      block
        real(R8), allocatable :: hIGLOOref(:), Lvmat(:)
        integer :: it
        allocate(hIGLOOref(size(IGLOO%material)), Lvmat(size(IGLOO%material)))
        do im = 1, size(IGLOO%material)
          Lvmat(im) = IGLOO%material(im)%Lv
          if ( IGLOO%material(im)%cpVariable ) then
            associate ( ht => IGLOO%material(im)%hTab )
              it = max(lbound(ht,1), min(int(T_datum_ref), ubound(ht,1)-1))
              hIGLOOref(im) = ht(it) + (ht(it+1)-ht(it))*(T_datum_ref - real(it,R8))
            end associate
          else
            hIGLOOref(im) = IGLOO%material(im)%cp * T_datum_ref + IGLOO%material(im)%hOff
          endif
          write(*,'(A,I0,A,ES14.7,A)') ' [MI2] IGLOO material ', im, ': table datum '// &
            trim(IGLOO%material(im)%hDatum)//', hOff = ', IGLOO%material(im)%hOff, ' J/kg'
        enddo
        call Build_Datum_Offsets( matnames, size(IGLOO%material), hIGLOOref, Lvmat )
        deallocate(hIGLOOref, Lvmat)
      end block

      deallocate(matnames)
    end block
    ! under-relaxation store, allocated whether or not relaxation is on
    call Initialize_Residual_Store( MOSE_simulation%domain(1), oldResidual )
  endif

  ! baseline gas inventory for the closed-domain mass and energy gates reported at the end
  block
    use MOSE_Global_m, only: nsc
    allocate(massGas0(nsc), massGas1(nsc))
    call Gas_Species_Mass  ( MOSE_simulation%domain(1), massGas0 )
    call Gas_Total_Energy  ( MOSE_simulation%domain(1), enGas0 )
    call Gas_Mean_Temperature( MOSE_simulation%domain(1), Tgas0, mGas0, cvGas0 )
  end block

  call Report_Setup()

  ! ---------------------------------------------------------------------------
  ! Coupling cycles
  ! ---------------------------------------------------------------------------
  ! Outer loop over cycles (one entry of the MOSE-iter / ICE-iter / IGLOO-iter lists), inner loop
  ! of N_Coupling(n_cycle) coupling steps. Solve order: MOSE, IGLOO, ICE.
  N_cycles = size(hydra_obj_sim_param%MOSE_iter)
  n_cycle  = 0

  do while ( n_cycle < N_cycles )
    n_cycle = n_cycle + 1
    write(*,*)
    write(*,'(A,I0,A,I0)') ' HYDRA-MI2: CYCLE ', n_cycle, ' / ', N_cycles

    n_couple = 0
    TODO = 1
    MOSE_obj_sim_param%TODO = TODO
    if (hydra_MI2_obj%run_ICE) ICE_obj_sim_param%TODO = TODO

    do while ( TODO <= 2 )

      n_couple = n_couple + 1

      ! Hot re-read of input.ini, so a long run can be steered without restarting
      if ( hydra_obj_io%ini_diter > 0 ) then
        if ( mod(n_couple, hydra_obj_io%ini_diter) == 0 ) call Read_Inifile_MI2_Runtime()
      endif

      ! --- MOSE ---------------------------------------------------------------
      ! announce sources frozen this cycle before the block that applies them
      call Report_Frozen_Sources()

      ! domain clock at block start, for the mass and energy audit
      t_block0 = MOSE_simulation%domain(1)%time

      do i_MOSE = 1, hydra_obj_sim_param%MOSE_iter(n_cycle)
        if ( MOSE_obj_sim_param%TODO <= 2 ) then
          call MOSE%solve( MOSE_simulation, MOSE_external_routine )
          if ( MOSE_obj_sim_param%TODO <= 2 ) call MOSE%postprocess( MOSE_simulation )
        endif
      enddo

      ! audit the mass and energy this block injected (source branch only: the Eulerian branch
      ! is mass-inert)
      if ( IGLOO_has_sources .and. .not. IGLOO_src_eulerian ) then
        call Accumulate_Injected_Mass  ( IGLOO, MOSE_simulation%domain(1)%time - t_block0 )
        call Accumulate_Injected_Energy( IGLOO, MOSE_simulation%domain(1)%time - t_block0 )
      elseif ( IGLOO_has_sources ) then
        nAuditSkipped = nAuditSkipped + 1
      endif

      ! --- IGLOO --------------------------------------------------------------
      if ( hydra_MI2_obj%run_IGLOO ) then
        if ( hydra_MI2_obj%IGLOO_iter(n_cycle) /= 0 ) then
          call Set_IGLOO_Switches(n_cycle)

          ! reset_state is mandatory before every solve except the first (IGLOO%setup ran it);
          ! it drops the sources MOSE has just consumed
          if ( IGLOO_first_sweep ) then
            IGLOO_first_sweep = .false.
          else
            call IGLOO%reset_state()
          endif

          ! gas state in, then solve
          call MOSExIGLOO( MOSE_simulation%domain(1), IGLOO )

          call IGLOO%solve()

          ! IGLOO%source / %euler are valid from here until the next reset_state
          IGLOO_has_sources  = .true.
          IGLOO_src_eulerian = eulCoupling
          IGLOO_src_cycle    = n_cycle
          IGLOO_src_step     = n_couple
        endif
      endif

      ! --- ICE ----------------------------------------------------------------
      if ( hydra_MI2_obj%run_ICE ) then
        ! gas state in, then solve (ICE computes its sources from gas_phase internally)
        if ( hydra_MI2_obj%ICE_iter(n_cycle) /= 0 ) then
          call MOSExICE( MOSE_simulation%domain(1), ICE_simulation%domain(1) )
        endif

        do i_ICE = 1, hydra_MI2_obj%ICE_iter(n_cycle)
          if ( ICE_obj_sim_param%TODO <= 2 ) then
            call ICE%solve( ICE_simulation, Dummy_Function )
            if ( ICE_obj_sim_param%TODO <= 2 ) call ICE%postprocess( ICE_simulation )
            ! cond_phase%source is valid after a solve with twcoupled set
            ICE_has_sources = .true.
            ICE_src_cycle   = n_cycle
            ICE_src_step    = n_couple
          endif
        enddo

        ! warn once if ICE's mass channel, dropped by ICExMOSE, becomes non-zero
        if ( ICE_has_sources ) call Check_ICE_Mass_Channel( ICE_simulation%domain(1) )
      endif

      ! --- Exit conditions ----------------------------------------------------
      ! the coupling step count is the only thing that ends a cycle: a solver's TODO = 3 is
      ! reset to 1, not propagated
      if ( n_couple >= hydra_obj_sim_param%N_Coupling(n_cycle) ) TODO = 3

      if ( TODO <= 2 ) then
        MOSE_obj_sim_param%TODO = 1
        if (hydra_MI2_obj%run_ICE) ICE_obj_sim_param%TODO = 1
      endif

    enddo

  enddo

  ! ---------------------------------------------------------------------------
  ! Final postprocess
  ! ---------------------------------------------------------------------------
  write(*,*)
  write(*,'(A)') ' HYDRA-MI2: FINAL POSTPROCESS'
  MOSE_obj_sim_param%TODO = 2
  call MOSE%postprocess( MOSE_simulation )
  if (hydra_MI2_obj%run_ICE) then
    ICE_obj_sim_param%TODO = 2
    call ICE%postprocess( ICE_simulation )
  endif
  if (hydra_MI2_obj%run_IGLOO) call IGLOO%writeout()

  if ( hydra_MI2_obj%run_IGLOO .and. has_stoichiometry .and. nMassLimited > 0 ) then
    write(*,*)
    write(*,'(A,I0,A)') ' [MI2] WARNING: the consumed-species floor guard clamped ', nMassLimited, &
      ' (cell,species) extractions.'
    write(*,'(A)') '       A gas species ran out where a reaction wanted to consume it, so the'
    write(*,'(A)') '       reaction was MASS-LIMITED: this is not the stoichiometry you asked for.'
  endif

  write(*,*)
  ! Mass-coupling audit: injected mass per species and the closed-domain balance against the
  ! change in gas mass. Meaningful only for a time-accurate run.
  if ( allocated(massInjected) ) then
    write(*,*)
    write(*,'(A)') ' [MI2] mass injected into the gas over the run [kg]'
    if ( nAuditSkipped > 0 ) then
      write(*,'(A,I0,A)') '   [WARNING] ', nAuditSkipped, ' coupling step(s) used the EULERIAN branch, '// &
                          'which is mass-inert, and are excluded.'
    endif
    if ( .not. hydra_obj_sim_param%time_accurate ) then
      write(*,'(A)') '   (time-accurate is OFF -- these are pseudo-time integrals, not physical mass)'
    endif
    block
      integer :: is
      do is = 1, size(massInjected)
        if ( abs(massInjected(is)) > 0._R8 ) &
          write(*,'(A,I0,A,ES14.6)') '   species ', is, ' : ', massInjected(is)
      enddo
      write(*,'(A,ES14.6)') '   TOTAL              : ', sum(massInjected)
    end block

    ! closed-domain balance: injected vs gas gained, both measured at the same instant
    call Gas_Species_Mass( MOSE_simulation%domain(1), massGas1 )
    write(*,*)
    write(*,'(A)') ' [MI2] closed-domain mass balance  (valid only if NO boundary flux)'
    write(*,'(A)') '   species        injected           gas gained         rel.err'
    block
      integer  :: is
      real(R8) :: dgas, den, floor_
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
    end block
    write(*,*)
  endif

  ! Energy bookkeeping (datum-blind) and the enthalpy-datum diagnostic
  if ( energyAudited ) then
    call Gas_Total_Energy    ( MOSE_simulation%domain(1), enGas1 )
    call Gas_Mean_Temperature( MOSE_simulation%domain(1), Tgas1, mGas1, cvGas1 )
    write(*,*)
    write(*,'(A)') ' [MI2] closed-domain ENERGY bookkeeping  (valid only if NO boundary flux)'
    if ( .not. hydra_obj_sim_param%time_accurate ) then
      write(*,'(A)') '   (time-accurate is OFF -- pseudo-time integral, not physical energy)'
    endif
    if ( nAuditSkipped > 0 ) then
      write(*,'(A,I0,A)') '   [WARNING] ', nAuditSkipped, ' coupling step(s) used the EULERIAN branch '// &
                          'and are NOT in these totals.'
      write(*,'(A)')      '             The balance below is therefore incomplete -- it is only a gate '// &
                          'for a pure source-mode run.'
    endif
    write(*,'(A,ES16.8)') '   handed to gas  (sourceEn + datum) [J] : ', energyInjected
    write(*,'(A,ES16.8)') '     of which the datum correction   [J] : ', energyDatumPart
    write(*,'(A,ES16.8)') '   gas total energy gained          [J] : ', enGas1 - enGas0
    block
      real(R8) :: den
      den = max(abs(energyInjected), abs(enGas1-enGas0))
      if (den > 0._R8) &
        write(*,'(A,ES11.3)') '   rel.err (bookkeeping only, datum-blind) : ', &
                              abs((enGas1-enGas0) - energyInjected)/den
    end block

    write(*,*)
    write(*,'(A)') ' [MI2] enthalpy-datum diagnostic'
    write(*,'(A,F12.6,A,F12.6,A,ES11.3)') '   mass-weighted gas T:  ', Tgas0, ' -> ', Tgas1, &
                                          '   delta = ', Tgas1-Tgas0
    if ( abs(energyDatumPart) > 0._R8 .and. mGas1*cvGas1 > 0._R8 ) then
      write(*,'(A,F12.4,A,ES12.4)') '   mass-weighted frozen cv [J/kg/K] : ', cvGas1, &
                                    '   gas mass [kg] : ', mGas1
      write(*,'(A,ES16.8)') '   energy the gas would have received WITHOUT it [J] : ', &
                            energyInjected - energyDatumPart
      write(*,'(A,F12.6)')  '   => predicted delta-T without the correction  [K] : ', &
                            (Tgas1-Tgas0) - energyDatumPart/(mGas1*cvGas1)
      write(*,'(A)') '   That prediction holds the evaporated mass FIXED. Re-running with'
      write(*,'(A)') '   [MI2-Coupling] datum-correction = false must flip the sign of delta-T, but'
      write(*,'(A)') '   will not reproduce the value: the gas state feeds back into the evaporation'
      write(*,'(A)') '   rate, so the two runs evaporate different amounts. The falsifiable invariant'
      write(*,'(A)') '   is the per-kg offset printed at setup, not this linearised delta-T.'
      write(*,'(A)') '   Scope: this validates the COUPLER''s datum bookkeeping, not the evaporation'
      write(*,'(A)') '   model -- IGLOO''s constant Lv cancels on both sides of the comparison.'
    endif
    write(*,*)
  endif

  write(*,'(A)') ' GAME OVER'
  write(*,*)


contains

  !> Source-term injection point, called by MOSE once per Runge-Kutta stage with the current
  !  domain, inside its parallel region (orphaned !$OMP DO only). Adds IGLOOxMOSE then ICExMOSE,
  !  each gated on its solver being active and having produced sources.
  subroutine MOSE_external_routine(MOSE_domain)
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain

    if ( hydra_MI2_obj%run_IGLOO .and. IGLOO_has_sources ) then
      call IGLOOxMOSE( MOSE_domain, IGLOO, eulCoupling, hydra_MI2_obj%relaxation,          &
                       i_MOSE, hydra_obj_sim_param%MOSE_iter(n_cycle),                     &
                       hydra_MI2_obj%relCoeff(n_cycle), oldResidual )
    endif

    ! MOSE_domain is the callback's argument, not the module-level simulation
    if ( hydra_MI2_obj%run_ICE .and. ICE_has_sources ) then
      call ICExMOSE( MOSE_domain, ICE_simulation%domain(1) )
    endif

  end subroutine MOSE_external_routine


  !> ICE accepts an external routine but never calls it; this satisfies the argument list.
  subroutine Dummy_Function
  end subroutine Dummy_Function


  !> Announce sources frozen this cycle (solver enabled, X-iter(n_cycle) = 0, sources already
  !  produced): at the first coupling step of the cycle, then at [MOSE-IO] shell-diter cadence.
  subroutine Report_Frozen_Sources()
    implicit none
    integer :: iterNow, cadence
    logical :: iglooFrozen, iceFrozen, due

    ! X_iter has size 1 when its solver is off
    iglooFrozen = .false.
    if ( hydra_MI2_obj%run_IGLOO .and. IGLOO_has_sources ) iglooFrozen = ( hydra_MI2_obj%IGLOO_iter(n_cycle) == 0 )
    iceFrozen = .false.
    if ( hydra_MI2_obj%run_ICE .and. ICE_has_sources )     iceFrozen   = ( hydra_MI2_obj%ICE_iter(n_cycle)   == 0 )
    if ( .not. (iglooFrozen .or. iceFrozen) ) return

    iterNow = MOSE_obj_sim_param%iter_general
    cadence = max(MOSE_obj_io%shell_diter, 1)
    due = ( n_couple == 1 ) .or. ( frozenPrintIter < 0 ) .or. ( iterNow/cadence > frozenPrintIter/cadence )
    if ( .not. due ) return
    frozenPrintIter = iterNow

    if ( iglooFrozen ) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      ' [MI2] IGLOO source FROZEN this cycle (IGLOO-iter = 0): the gas keeps the sources produced in cycle ', &
      IGLOO_src_cycle, ', coupling step ', IGLOO_src_step, '  [cycle ', n_cycle, ', step ', n_couple, &
      ', MOSE iter ', iterNow, ']'
    if ( iceFrozen ) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      ' [MI2] ICE source FROZEN this cycle (ICE-iter = 0): the gas keeps the sources produced in cycle ', &
      ICE_src_cycle, ', coupling step ', ICE_src_step, '  [cycle ', n_cycle, ', step ', n_couple, &
      ', MOSE iter ', iterNow, ']'

  end subroutine Report_Frozen_Sources


  !> Select IGLOO's coupling mode for this cycle: IGLOO-iter = 2 -> Eulerian (eulSwitch),
  !  1 -> source (srcSwitch). On the last cycle that runs IGLOO both switches are on for the writeout.
  subroutine Set_IGLOO_Switches(cycle_now)
    implicit none
    integer, intent(in) :: cycle_now
    integer :: i
    logical :: is_last

    if ( hydra_MI2_obj%IGLOO_iter(cycle_now) == 2 ) then
      eulCoupling      = .true.
      IGLOO%eulSwitch  = .true.
      IGLOO%srcSwitch  = .false.
    else
      eulCoupling      = .false.
      IGLOO%eulSwitch  = .false.
      IGLOO%srcSwitch  = .true.
    endif

    ! Is this the last cycle that runs IGLOO at all?
    is_last = .true.
    do i = cycle_now + 1, N_cycles
      if ( hydra_MI2_obj%IGLOO_iter(i) /= 0 ) then
        is_last = .false.
        exit
      endif
    enddo
    if ( is_last ) then
      IGLOO%eulSwitch = .true.
      IGLOO%srcSwitch = .true.
    endif

    ! the Eulerian branch carries no mass: warn when a stoichiometry is configured
    if ( eulCoupling .and. has_stoichiometry ) then
      write(*,'(A)') ' [MI2] WARNING: IGLOO-iter = 2 (Eulerian coupling) carries momentum and'
      write(*,'(A)') '       energy ONLY. getSourceTerms has no mass component, so the'
      write(*,'(A)') '       [MI2-Coupling] stoichiometry is NOT applied this cycle.'
      write(*,'(A)') '       Use IGLOO-iter = 1 (source coupling) for mass transfer.'
    endif

  end subroutine Set_IGLOO_Switches


  !> Echo the coupling schedule exactly as MI2 resolved it, after clamping and defaulting.
  subroutine Report_MI2_Config()
    implicit none
    integer :: i

    write(*,*)
    write(*,'(A)') ' HYDRA-MI2 coupling schedule'
    write(*,'(A,I0)')  ' -  cycles (from MOSE-iter)   --> ', size(hydra_obj_sim_param%MOSE_iter)
    write(*,'(A,L1)')  ' -  run ICE                   --> ', hydra_MI2_obj%run_ICE
    write(*,'(A,L1)')  ' -  run IGLOO                 --> ', hydra_MI2_obj%run_IGLOO
    write(*,'(A,L1)')  ' -  relaxation                --> ', hydra_MI2_obj%relaxation

    write(*,'(A)',advance='no') ' -  MOSE-iter                 --> '
    do i = 1, size(hydra_obj_sim_param%MOSE_iter)
      write(*,'(I0,1X)',advance='no') hydra_obj_sim_param%MOSE_iter(i)
    enddo
    write(*,*)

    write(*,'(A)',advance='no') ' -  N-coupling                --> '
    do i = 1, size(hydra_obj_sim_param%N_Coupling)
      write(*,'(I0,1X)',advance='no') hydra_obj_sim_param%N_Coupling(i)
    enddo
    write(*,*)

    if (hydra_MI2_obj%run_ICE) then
      write(*,'(A)',advance='no') ' -  ICE-iter                  --> '
      do i = 1, size(hydra_MI2_obj%ICE_iter)
        write(*,'(I0,1X)',advance='no') hydra_MI2_obj%ICE_iter(i)
      enddo
      write(*,*)
    endif

    if (hydra_MI2_obj%run_IGLOO) then
      write(*,'(A)',advance='no') ' -  IGLOO-iter (0/1/2)        --> '
      do i = 1, size(hydra_MI2_obj%IGLOO_iter)
        write(*,'(I0,1X)',advance='no') hydra_MI2_obj%IGLOO_iter(i)
      enddo
      write(*,*)
      write(*,'(A)',advance='no') ' -  relaxation-coeff          --> '
      do i = 1, size(hydra_MI2_obj%relCoeff)
        write(*,'(F6.3,1X)',advance='no') hydra_MI2_obj%relCoeff(i)
      enddo
      write(*,*)
    endif
    write(*,*)

  end subroutine Report_MI2_Config


  !> What actually got built, after setup.
  subroutine Report_Setup()
    implicit none

    write(*,*)
    write(*,'(A)') ' HYDRA-MI2 solvers'
    write(*,'(A,I0)') ' -  MOSE  blocks              --> ', MOSE_simulation%domain(1)%nb
    if (hydra_MI2_obj%run_ICE) then
      write(*,'(A,I0)') ' -  ICE   blocks              --> ', ICE_simulation%domain(1)%nb
    endif
    if (hydra_MI2_obj%run_IGLOO) then
      write(*,'(A,I0)') ' -  IGLOO geometry blocks     --> ', size(IGLOO%geoblock)
      write(*,'(A,I0)') ' -  IGLOO materials           --> ', size(IGLOO%material)
    endif
    write(*,*)

  end subroutine Report_Setup

end program hydra_MI2_program
