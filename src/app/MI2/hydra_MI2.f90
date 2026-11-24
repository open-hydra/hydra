!>
!> hydra-MI2 -- MOSE (Eulerian gas) + ICE (Eulerian condensed) + IGLOO (Lagrangian) coupler.
!> OpenMP only: the root CMakeLists refuses USE_MPI for MI2.
!>
!> A run is a list of cycles, one per entry of the MOSE-iter / N-coupling / ICE-iter / IGLOO-iter
!> lists of [HYDRA-Parameters]. A cycle is N-coupling coupling steps, and a coupling step is
!>   MOSE  : MOSE-iter iterations, the condensed solvers' sources injected at every Runge-Kutta stage
!>   IGLOO : one sweep on the new gas state   (IGLOO-iter: 0 none, 1 source mode, 2 Eulerian mode)
!>   ICE   : ICE-iter iterations on the new gas state
!> A solver that does not run in a cycle leaves its last sources in force: MOSE keeps reading them.
!>
!> This file is the schedule. The input and its refusals are in hydraMI2_Read_Ini; the exchange, and what lands
!> between solver steps (ICE's removals, IGLOO's smoke), in hydraMI2_Mod_Coupling; the restart rules and ICE's
!> hand-off in hydraMI2_Restart; the books and the log's report lines in hydraMI2_Report.
!>
program hydra_MI2_program
  use iso_fortran_env, only: R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use hydra_Config_Types_m,    only: hydra_obj_sim_param
  use hydraMI2_Config_Types_m, only: hydra_MI2_obj
  use hydraMI2_Build_Stamp,    only: Print_Build_Stamp
  use hydraMI2_Read_Ini,       only: Read_Inifile_MI2, Read_Inifile_MI2_Runtime, Check_MOSE_Grid_Levels,    &
                                     Check_ICE_Grid_Levels, Check_Coupled_Inputs
  use hydraMI2_Mod_Import,     only: Assign_Phases
  use hydraMI2_Mod_Coupling,   only: source_state_t, removals_t, smoke_t, handRates, smokeIGOn, mdotHanded, nHanded, &
                                     Setup_ICE_Coupling, Setup_IGLOO_Coupling, Set_IGLOO_Switches, Gas_Field_For_IGLOO, &
                                     MOSExIGLOO, MOSExICE, IGLOOxMOSE, ICExMOSE, Check_IGLOO_Sources, Check_ICE_Sources, &
                                     Open_Removals, Close_Removals, ICE_Step_Clocks, ICE_Block_Mean,                   &
                                     IGLOO_Smoke_Accumulate, IGLOO_Smoke_Deliver
  use hydraMI2_Restart,        only: Setup_Restart, Remove_Stale_Handoff, Note_Outputs, Save_Handoff_If_Fields_Current, &
                                     Write_ICE_Exchange, Made_No_Iteration
  use hydraMI2_Report,         only: Report_Start, Report_Frozen_Sources, Report_Values, Report_Inlet,                &
                                     Report_Block_Start, Report_Block_End, Report_End

  use MOSE_Advanced_Types_m, only: MOSE_simulation_type, MOSE_domain_type
  use MOSE_Config_Types_m,   only: MOSE_obj_sim_param => obj_sim_param
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
  !> The sources each condensed solver last produced: what MOSE reads at every Runge-Kutta stage.
  type(source_state_t) :: IGLOO_src, ICE_src
  !> Where the run is: the cycle, the coupling step of the cycle, the MOSE iteration of the step.
  integer :: n_cycle, n_couple, i_MOSE
  !> The solver made the iteration the schedule gives as its last of the run (its postprocess is the end-of-run output).
  logical :: MOSE_ended = .false., ICE_ended = .false.

  ! ---------------------------------------------------------------------------
  ! Setup
  ! ---------------------------------------------------------------------------
  call Print_Build_Stamp()
  call Set_Threads()

  call Read_Inifile_MI2()      ! input.ini -> hydra_obj_*, hydra_MI2_obj; echoes the schedule
  call Assign_Phases()         ! every INPUT/*phase.txt -> the solver that owns it
  call Setup_Solvers()

  ! ---------------------------------------------------------------------------
  ! Coupling cycles
  ! ---------------------------------------------------------------------------
  do n_cycle = 1, size(hydra_obj_sim_param%MOSE_iter)
    write(*,*)
    write(*,'(A,I0,A,I0)') ' HYDRA-MI2: CYCLE ', n_cycle, ' / ', size(hydra_obj_sim_param%MOSE_iter)

    ! the coupling-step count alone ends a cycle: a solver that meets its own stop condition ends
    ! its part of the step, not the run
    n_couple = 0
    do
      n_couple = n_couple + 1
      call Read_Inifile_MI2_Runtime( n_couple )

      call Advance_MOSE()
      if ( hydra_MI2_obj%run_IGLOO ) call Advance_IGLOO()
      if ( hydra_MI2_obj%run_ICE )   call Advance_ICE()

      if ( n_couple >= hydra_obj_sim_param%N_Coupling(n_cycle) ) exit
    enddo
  enddo

  call Final_Postprocess()

  write(*,'(A)') ' GAME OVER'
  write(*,*)


contains

  !> Thread count, reported and handed to the solvers.
  subroutine Set_Threads()
    implicit none
    integer :: nthreads

#if defined (_OPENMP)
    ! the team size of a parallel region, written by one thread
    !$omp parallel
    !$omp single
    nthreads = OMP_GET_NUM_THREADS()
    !$omp end single
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

  end subroutine Set_Threads


  !> Set up MOSE, then each enabled condensed solver on MOSE's gas field with its side of the coupler; then what a
  !  restart reads back, the checks that need every solver's inputs and the baseline of the books.
  subroutine Setup_Solvers()
    use Lib_ORION_data, only: orion_data
    implicit none
    type(orion_data) :: gas_for_igloo

    ! ICE must know it is two-way coupled before setup (it then fills cond_phase%source)
    ICE_obj_sim_param%twcoupled = hydra_MI2_obj%run_ICE

    call MOSE%setup( MOSE_simulation )
    call Check_MOSE_Grid_Levels()

    if ( hydra_MI2_obj%run_ICE ) then
      call ICE%setup( ICE_simulation, MOSE_simulation%IOfield(1) )
      call Check_ICE_Grid_Levels()
    endif

    call Setup_Restart( ICE_simulation, ICE_src%ready )

    if ( hydra_MI2_obj%run_ICE ) call Setup_ICE_Coupling()

    ! IGLOO%setup includes the reset_state of the first sweep
    if ( hydra_MI2_obj%run_IGLOO ) then
      call Gas_Field_For_IGLOO( MOSE_simulation%domain(1), MOSE_simulation%IOfield(1), gas_for_igloo )
      call IGLOO%setup( gas_for_igloo )
      call Setup_IGLOO_Coupling( MOSE_simulation%domain(1), IGLOO )
    endif

    call Check_Coupled_Inputs( IGLOO )
    call Report_Start( MOSE_simulation, ICE_simulation, IGLOO )

    ! past every set-up refusal: a new ICE run drops the hand-off an earlier run left
    call Remove_Stale_Handoff()

  end subroutine Setup_Solvers


  !> MOSE's part of a coupling step: MOSE-iter iterations under the sources in force; around them, ICE's stored
  !  removals and IGLOO's smoke handed on, and the books and the report of what the block handed the gas.
  subroutine Advance_MOSE()
    implicit none
    type(removals_t) :: removals
    type(smoke_t)    :: smoke
    integer :: nIter, nRun, itermaxKept
    logical :: lastIter

    nIter = hydra_obj_sim_param%MOSE_iter(n_cycle)
    call Report_Frozen_Sources( IGLOO_src, ICE_src, n_cycle, n_couple )
    call Report_Block_Start( MOSE_simulation%domain(1) )

    ! the removals ICE stores reach the gas over this block (ICExMOSE): as the rate of the store, or 1/nIter of it per
    ! iteration when both solvers are time-accurate
    if ( hydra_MI2_obj%run_ICE ) &
      call Open_Removals( MOSE_simulation%domain(1), ICE_simulation%domain(1), ICE_src%ready, nIter, removals )

    if ( MOSE_ended .and. nIter > 0 ) call Resumed( 'MOSE', MOSE_ended )
    MOSE_obj_sim_param%TODO = 1
    nRun = 0
    do i_MOSE = 1, nIter
      if ( MOSE_obj_sim_param%TODO > 2 ) exit
      ! the run's last MOSE iteration ends on MOSE's own iteration limit (TODO = 3), as a run of MOSE alone: its
      ! postprocess is the end-of-run output
      lastIter = Last_Of_Run( hydra_obj_sim_param%MOSE_iter, i_MOSE )
      if ( lastIter ) then
        itermaxKept = MOSE_simulation%domain(1)%itermax
        MOSE_simulation%domain(1)%itermax = MOSE_obj_sim_param%iter_from_call + 1
      endif
      call MOSE%solve( MOSE_simulation, MOSE_external_routine )
      if ( lastIter ) then
        MOSE_simulation%domain(1)%itermax = itermaxKept
        MOSE_ended = .true.
      endif
      nRun = nRun + 1
      if ( smokeIGOn .and. IGLOO_src%ready .and. .not. IGLOO_src%eulerian ) &
        call IGLOO_Smoke_Accumulate( MOSE_simulation%domain(1), IGLOO )
      if ( MOSE_obj_sim_param%TODO <= 2 ) call MOSE%postprocess( MOSE_simulation )
      call Note_Outputs()
    enddo

    if ( hydra_MI2_obj%run_ICE ) call Close_Removals( MOSE_simulation%domain(1), ICE_simulation%domain(1), nRun, removals )

    ! the smoke IGLOO's burning made over this block, which the gas did not get, to ICE
    if ( smokeIGOn ) call IGLOO_Smoke_Deliver( MOSE_simulation%domain(1), ICE_simulation%domain(1), IGLOO, &
                                               IGLOO_src%ready .and. .not. IGLOO_src%eulerian, smoke )

    call Report_Block_End( MOSE_simulation%domain(1), ICE_simulation, IGLOO, IGLOO_src, ICE_src, n_cycle, n_couple, &
                           removals, smoke )

  end subroutine Advance_MOSE


  !> Source-term injection point, called by MOSE once per Runge-Kutta stage with the current domain, inside its
  !  parallel region (orphaned !$OMP DO only): the sources in force of each condensed solver go into MOSE's residual.
  !  A coarse level of MOSE's multigrid cycle receives them in the residual MOSE restricts from the finest grid.
  subroutine MOSE_external_routine(MOSE_domain)
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain

    if ( MOSE_domain%mg_level /= 1 ) return

    if ( IGLOO_src%ready ) then
      call IGLOOxMOSE( MOSE_domain, IGLOO, IGLOO_src%eulerian, i_MOSE, hydra_obj_sim_param%MOSE_iter(n_cycle), &
                       merge(hydra_MI2_obj%relCoeff(n_cycle), 1._R8, hydra_MI2_obj%relaxation) )
    endif

    ! MOSE_domain is the callback's argument, not the module-level simulation
    if ( ICE_src%ready ) call ICExMOSE( MOSE_domain, ICE_simulation%domain(1) )

  end subroutine MOSE_external_routine


  !> IGLOO's part of a coupling step: one sweep on the gas state MOSE has just reached, in the coupling mode of the
  !  cycle. Nothing when IGLOO-iter = 0: its sources stay frozen.
  subroutine Advance_IGLOO()
    implicit none
    logical :: eulerian

    if ( hydra_MI2_obj%IGLOO_iter(n_cycle) == 0 ) return

    call Set_IGLOO_Switches( IGLOO, n_cycle, eulerian )

    ! reset_state drops the sources MOSE has just consumed: mandatory before every sweep but the first, which
    ! IGLOO%setup prepared
    if ( IGLOO_src%ready ) call IGLOO%reset_state()

    call MOSExIGLOO( MOSE_simulation%domain(1), IGLOO )
    call Report_Inlet( mdotHanded, nHanded )
    call IGLOO%solve()

    ! IGLOO%source / %euler are valid from here until the next reset_state
    IGLOO_src = source_state_t( ready=.true., eulerian=eulerian, cycle=n_cycle, step=n_couple )

    ! refuse a non-finite value before MOSE applies it: the source-mode rates, or the Eulerian fields
    call Check_IGLOO_Sources( MOSE_simulation%domain(1), IGLOO, eulerian )

  end subroutine Advance_IGLOO


  !> ICE's part of a coupling step: ICE-iter iterations on the gas state MOSE has just reached (ICE computes its
  !  sources from gas_phase internally); then, whether or not ICE iterated, the checks of what it hands the gas and
  !  its hand-off when both fields on disk are of this step.
  subroutine Advance_ICE()
    implicit none
    real(R8) :: t_ice0
    integer  :: i_ICE, itermaxKept
    logical  :: lastIter

    if ( hydra_MI2_obj%ICE_iter(n_cycle) /= 0 ) call MOSExICE( MOSE_simulation%domain(1), ICE_simulation%domain(1) )

    t_ice0 = ICE_simulation%domain(1)%time
    if ( ICE_ended .and. hydra_MI2_obj%ICE_iter(n_cycle) > 0 ) call Resumed( 'ICE', ICE_ended )
    ICE_obj_sim_param%TODO = 1
    do i_ICE = 1, hydra_MI2_obj%ICE_iter(n_cycle)
      if ( ICE_obj_sim_param%TODO > 2 ) exit
      ! the run's last ICE iteration ends on ICE's own iteration limit, as MOSE's
      lastIter = Last_Of_Run( hydra_MI2_obj%ICE_iter, i_ICE )
      if ( lastIter ) then
        itermaxKept = ICE_simulation%domain(1)%itermax
        ICE_simulation%domain(1)%itermax = ICE_obj_sim_param%iter_from_call + 1
      endif
      call ICE%solve( ICE_simulation, Dummy_Function )
      if ( lastIter ) then
        ICE_simulation%domain(1)%itermax = itermaxKept
        ICE_ended = .true.
      endif
      ! the time the removals of this iteration were made over (steady), the block's mean source (time-accurate)
      call ICE_Step_Clocks( ICE_simulation%domain(1) )
      if ( ICE_obj_sim_param%TODO <= 2 ) call ICE%postprocess( ICE_simulation )
      call Note_Outputs()
      ! after a solve with twcoupled set, cond_phase%source holds the rates of ICE's step
      ICE_src = source_state_t( ready=.true., cycle=n_cycle, step=n_couple )
    enddo

    call ICE_Block_Mean( ICE_simulation%domain(1) )
    ! the ICE time the smoke's rates were integrated over
    if ( handRates .and. hydra_MI2_obj%ICE_iter(n_cycle) /= 0 ) &
      call Report_Values( 'ice-time', [t_ice0, ICE_simulation%domain(1)%time] )

    ! refuse a mass source of a material without a record, and a non-finite source, before MOSE applies them
    if ( ICE_src%ready ) call Check_ICE_Sources( MOSE_simulation%domain(1), ICE_simulation%domain(1) )

    call Save_Handoff_If_Fields_Current( ICE_simulation, ICE_src%ready )

  end subroutine Advance_ICE


  !> Each solver's end-of-run output (TODO = 3): its last iteration's rows and the solution due there, written once,
  !  then its final solution; ICE's hand-off; IGLOO's fields of its last sweep; the books closed.
  subroutine Final_Postprocess()
    implicit none

    write(*,*)
    write(*,'(A)') ' HYDRA-MI2: FINAL POSTPROCESS'
    ! a restarted solver that made no iteration: its end-of-run output is the one it restarted from, not written again
    MOSE_obj_sim_param%TODO = 3
    if ( Made_No_Iteration('MOSE') ) then
      write(*,'(A)') ' [MI2] MOSE made no iteration in this run: its end-of-run output is the one it restarted from'
    else
      call MOSE%postprocess( MOSE_simulation )
    endif
    if ( hydra_MI2_obj%run_ICE ) then
      ICE_obj_sim_param%TODO = 3
      if ( Made_No_Iteration('ICE') ) then
        write(*,'(A)') ' [MI2] ICE made no iteration in this run: its end-of-run output is the one it restarted from'
      else
        call ICE%postprocess( ICE_simulation )
      endif
      call Write_ICE_Exchange( ICE_simulation, ICE_src%ready, MOSE_obj_sim_param%iter_general )
    endif
    ! IGLOO's fields are those of its last sweep: none when no cycle advanced it
    if ( IGLOO_src%ready ) call IGLOO%writeout()

    call Report_End( MOSE_simulation%domain(1) )

  end subroutine Final_Postprocess


  !> Whether iteration i of the current coupling step is the solver's last of the run: the step ends its cycle, i is
  !  the cycle's X-iter and no later cycle advances the solver (the schedule as last read).
  logical function Last_Of_Run ( Xiter, i )
    implicit none
    integer, intent(in) :: Xiter(:), i

    Last_Of_Run = n_couple >= hydra_obj_sim_param%N_Coupling(n_cycle) .and. i == Xiter(n_cycle) .and. &
                  all( Xiter(n_cycle+1:) <= 0 )

  end function Last_Of_Run


  !> A re-read of input.ini that advances a solver after the iteration the schedule had made its last: the run goes on,
  !  that iteration's rows and due solution were not written.
  subroutine Resumed ( solver, ended )
    implicit none
    character(len=*), intent(in)    :: solver
    logical,          intent(inout) :: ended

    write(*,'(A)') ' [WARNING] MI2: input.ini, read again, advances '//solver//' after the iteration the schedule had '// &
      'made its last: the rows and the solution due at that iteration were not written'
    ended = .false.

  end subroutine Resumed


  !> ICE accepts an external routine but never calls it; this satisfies the argument list.
  subroutine Dummy_Function
  end subroutine Dummy_Function

end program hydra_MI2_program
