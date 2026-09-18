program hydra_program
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use hydra_Parameters_m
  use hydra_Read_Ini,                only: Read_Inifile, Read_Inifile_Runtime
  use hydra_Config_Types_m
  use hydraAF_Mod_Import,            only: Assign_Phases
  use hydraAF_IO_Wall
  use hydraAF_Lib_Metrics,           only: Setup_Coupled_Metrics
  use hydraAF_Lib_Ghost,             only: Fill_Coupled_Ghost_Cell
  use hydraAF_Mod_dt,                only: Fix_Global_dt
  use hydraAF_Mod_BC_Coupled_Fluxes, only: BC_Coupled_Fluxes
  use FUSS_Advanced_Types_m, only: FUSS_simulation_type
  use FUSS_Config_Types_m,   only: FUSS_obj_sim_param => obj_sim_param
  use FUSS_Procedures_m,     only: FUSS_type
  use ARES_Advanced_Types_m, only: ARES_simulation_type
  use ARES_Config_Types_m,   only: ARES_obj_sim_param => obj_sim_param
  use ARES_Procedures_m,     only: ARES_type
  use FUSS_Mod_MPI,          only: FUSS_mpi_init => mpi_init_env, FUSS_mpi_finalize => mpi_finalize_env, &
                                   mpi_is_root, mpi_size_
  use ARES_Mod_MPI,          only: ARES_mpi_init => mpi_init_env
  use hydraAF_Mod_CouplingExchange, only: build_coupling_schedule, cleanup_coupling_schedule
  use hydra_Mod_Timers,      only: timer_run_begin, timer_summary, timer_report, &
                                   timer_step_begin, timer_step_end, &
                                   timer_fuss_begin, timer_fuss_end, &
                                   timer_fluid_begin, timer_fluid_end, &
                                   timer_coup_begin, timer_coup_end, &
                                   timer_fluid_name
#ifdef USE_MPI
  use FUSS_Mod_GhostExchange, only: FUSS_cleanup_ghost => cleanup_ghost_schedule
  use ARES_Mod_GhostExchange, only: ARES_cleanup_ghost => cleanup_ghost_schedule
#endif
  use IR_precision
  implicit none
  type(ARES_type)            :: ARES
  type(ARES_simulation_type) :: ARES_simulation
  type(FUSS_type)            :: FUSS
  type(FUSS_simulation_type) :: FUSS_simulation
  ! Local
  integer :: nthreads
  character(len=clen_) :: extension
  character(len=hlen_) :: hydrafile
  integer  :: n_couple, n_wall, n, nmin, nmax
  integer  :: TODO, ex_FUSS_TODO, ex_ARES_TODO
  integer  :: i_FUSS, i_ARES, FUSS_iter, ARES_iter
  integer  :: mgl, m
  real(R8) :: hydra_time, pp_time

  ! Initialize MPI (first call does MPI_INIT_THREAD; second just populates module state)
  call FUSS_mpi_init()
  call ARES_mpi_init()

#if defined (_OPENMP)
  !$omp parallel
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  FUSS_obj_sim_param%nthreads = nthreads
  ARES_obj_sim_param%nthreads = nthreads
  if (mpi_is_root) then
    write(*,'(A)')    ' Parallel execution'
    write(*,'(A)')    ' OpenMP:'
    write(*,'(A,I4)') ' -  Number of threads --> ', nthreads
  endif
#else
  if (mpi_is_root) write(*,'(A)')    ' Serial execution'
  FUSS_obj_sim_param%nthreads = 1
  ARES_obj_sim_param%nthreads = 1
#endif

#ifdef USE_MPI
  if (mpi_is_root) then
    write(*,'(A)')      ' MPI:'
    write(*,'(A,I4)')   ' -  Number of ranks   --> ', mpi_size_
  end if
#endif

  ! Populate hydra_objs
  call Read_Inifile ()

  ! Setup Time-Accurate
  hydra_time = 0.d0
  pp_time = 0.d0
  if ( hydra_obj_sim_param%time_accurate ) then
    FUSS_obj_sim_param%HYDRA_time_accurate = .true.
    ARES_obj_sim_param%HYDRA_time_accurate = .true.
  endif

  ! Setup Multigrid
  if (hydra_obj_multigrid%MGL > 1) then
    FUSS_obj_sim_param%HYDRA_MG = .true.
    ARES_obj_sim_param%HYDRA_MG = .true.
  endif

  ! Setup Solvers
  call Assign_Phases()
  call FUSS%setup( FUSS_simulation )
  call ARES%setup( ARES_simulation )

  ! Coupling Initializations
  if (hydra_obj_io%coupling_wall) then
    allocate(IOwall(1:hydra_obj_multigrid%MGL))
    do m = 1, hydra_obj_multigrid%MGL
      if (index(hydra_obj_io%coupling_format,'ascii')>0) then
        extension = '.tec'
      else
        extension = '.szplt'
      endif

      if (mpi_is_root) call Initialize_Wall_File( FUSS_simulation%domain(m), ARES_simulation%domain(m), extension, m )
      
      ! Build the cross-solver coupling schedule for this MG level BEFORE any
      ! setup routine that may iterate over type(103) BCs. The schedule builds
      ! the partner_of_* indices even in serial.
      call build_coupling_schedule( FUSS_simulation%domain(m), ARES_simulation%domain(m), m )
      call Setup_Coupled_Metrics( FUSS_simulation%domain(m), ARES_simulation%domain(m), m )
      call Fill_Coupled_Ghost_Cell( FUSS_simulation%domain(m), ARES_simulation%domain(m), m )
    enddo
  endif

  ! Solve & Postprocess
  TODO = 1
  FUSS_obj_sim_param%TODO = TODO
  ARES_obj_sim_param%TODO = TODO
  n_couple = 0
  n_wall = 0
  mgl = hydra_obj_multigrid%MGL

  ! Set-up is over: start the coupling clock
  timer_fluid_name = 'ARES'
  call timer_run_begin()

  do while ( TODO <= 2 )

    call timer_step_begin()

    ! Set iteration number for each solver
    if ( hydra_obj_sim_param%time_accurate ) then
      FUSS_iter = 1
      ARES_iter = 1
    elseif (hydra_obj_multigrid%MGL > 1) then
      FUSS_iter = hydra_obj_sim_param%FUSS_iter(1)
      ARES_iter = hydra_obj_sim_param%ARES_iter(1)
    else
      nmin = 0
      nmax = hydra_obj_sim_param%N_Coupling(1)
      FUSS_iter = hydra_obj_sim_param%FUSS_iter(1)
      ARES_iter = hydra_obj_sim_param%ARES_iter(1)
      if (size(hydra_obj_sim_param%N_Coupling) > 1) then 
        do n = 1, size(hydra_obj_sim_param%N_Coupling)
          if (n_couple >= nmin .and. n_couple < nmax ) then
            FUSS_iter = hydra_obj_sim_param%FUSS_iter(n)
            ARES_iter = hydra_obj_sim_param%ARES_iter(n)
          endif
          nmin = nmax
          nmax = hydra_obj_sim_param%N_Coupling(n+1)
        enddo
      endif
    endif

    ! Update hydra input.ini options
    if ( mod(n_couple+1, hydra_obj_io%ini_diter) == 0d0 ) then
      call Read_Inifile_Runtime()
      if ( hydra_obj_sim_param%time_accurate ) then
        FUSS_obj_sim_param%HYDRA_time_accurate = .true.
        ARES_obj_sim_param%HYDRA_time_accurate = .true.
      else
        FUSS_obj_sim_param%HYDRA_time_accurate = .false.
        ARES_obj_sim_param%HYDRA_time_accurate = .false.
      endif
    endif

    ! Set Global dt for time-accurate
    if ( hydra_obj_sim_param%time_accurate ) then
      if ( hydra_obj_multigrid%MGL > 1 ) then
        if (mpi_is_root) write(*,*) 'Multi-Grid can not be used for time-accurate simulations for now!'
        stop
      else
        call Fix_Global_dt ( FUSS_simulation, ARES_simulation, hydra_time )
      endif
    endif

    ! Solve FUSS
    do i_FUSS = 1, FUSS_iter
      if ( FUSS_obj_sim_param%TODO <= 2 ) then
        if (hydra_obj_io%coupling_wall) then
          call timer_coup_begin()
          call BC_Coupled_Fluxes ( FUSS_simulation%domain(mgl), ARES_simulation%domain(mgl), mgl )
          call timer_coup_end()
        endif
        call timer_fuss_begin()
        call FUSS%solve( FUSS_simulation, Dummy_Function )
        if ( FUSS_obj_sim_param%TODO <= 2 ) call FUSS%postprocess( FUSS_simulation )
        call timer_fuss_end()
      endif
    enddo

    ! Solve ARES
    do i_ARES = 1, ARES_iter
      if ( ARES_obj_sim_param%TODO <= 2 ) then
        if (hydra_obj_io%coupling_wall) then
          call timer_coup_begin()
          call BC_Coupled_Fluxes ( FUSS_simulation%domain(mgl), ARES_simulation%domain(mgl), mgl )
          call timer_coup_end()
        endif
        call timer_fluid_begin()
        call ARES%solve( ARES_simulation, Dummy_Function )
        if ( ARES_obj_sim_param%TODO <= 2 ) call ARES%postprocess( ARES_simulation )
        call timer_fluid_end()
      endif
    enddo
    
    n_couple = n_couple + 1
    
    ! Intermediate Postprocess
    if (hydra_obj_multigrid%MGL > 1 .and. n_couple == sum(hydra_obj_multigrid%iter_threshold(mgl:))) then
      if (mpi_is_root) print*, 'COUPLING-STEP POSTPROCESS - MG'
      if (hydra_obj_io%coupling_wall) then
        n_wall = n_wall+1
        if (mgl == 1) then
          hydrafile = 'hydra-wall'
          if (.not.hydra_obj_io%coupling_overwrite) hydrafile = trim(hydrafile)//trim(str(.true.,n_wall))
        else
          hydrafile = 'hydra-wall-level'//trim(str(.true.,mgl))
        endif
        call Write_Wall_Solution( FUSS_simulation%domain(mgl), ARES_simulation%domain(mgl), hydrafile, mgl )
      endif

      mgl = mgl - 1

    elseif ( mod(n_couple, hydra_obj_io%coupling_diter) == 0d0 .or. hydra_time >= pp_time + hydra_obj_io%coupling_dtime ) then
      if (mpi_is_root) print*, 'COUPLING-STEP POSTPROCESS'
      FUSS_obj_sim_param%HYDRA_postprocess = .true.
      ARES_obj_sim_param%HYDRA_postprocess = .true.
      ex_FUSS_TODO = FUSS_obj_sim_param%TODO
      ex_ARES_TODO = ARES_obj_sim_param%TODO
      FUSS_obj_sim_param%TODO = 2
      ARES_obj_sim_param%TODO = 2
      call FUSS%postprocess( FUSS_simulation )
      call ARES%postprocess( ARES_simulation )
      if (hydra_obj_io%coupling_wall) then
        n_wall = n_wall+1
        if (mgl == 1) then
          hydrafile = 'hydra-wall'
          if (.not.hydra_obj_io%coupling_overwrite) hydrafile = trim(hydrafile)//trim(str(.true.,n_wall))
        else
          hydrafile = 'hydra-wall-level'//trim(str(.true.,mgl))
        endif
        call Write_Wall_Solution( FUSS_simulation%domain(mgl), ARES_simulation%domain(mgl), hydrafile, mgl )
      endif
      FUSS_obj_sim_param%TODO = ex_FUSS_TODO
      ARES_obj_sim_param%TODO = ex_ARES_TODO
      FUSS_obj_sim_param%HYDRA_postprocess = .false.
      ARES_obj_sim_param%HYDRA_postprocess = .false.
      pp_time = hydra_time
    endif

    ! Check if Exit from loop
    if ( FUSS_obj_sim_param%TODO == 3 .and. ARES_obj_sim_param%TODO == 3 ) then
      TODO = 3
    endif
    if ( n_couple == hydra_obj_sim_param%N_Coupling(size(hydra_obj_sim_param%N_Coupling)) ) then
      TODO = 3
      FUSS_obj_sim_param%TODO = 3
      ARES_obj_sim_param%TODO = 3
    endif
    if ( hydra_obj_multigrid%MGL > 1 .and. mgl == 0 ) then
      TODO = 3
      FUSS_obj_sim_param%TODO = 3
      ARES_obj_sim_param%TODO = 3
    endif
    if ( hydra_obj_sim_param%time_accurate .and. hydra_time >= hydra_obj_sim_param%time_threshold) then
      TODO = 3
      FUSS_obj_sim_param%TODO = 3
      ARES_obj_sim_param%TODO = 3
    endif

    if ( TODO <= 2 ) then
      FUSS_obj_sim_param%TODO = 1
      ARES_obj_sim_param%TODO = 1
    endif

    ! Wall-clock report (collective: every rank runs the same step counter)
    call timer_step_end()
    if ( hydra_obj_io%timer_diter > 0 ) then
      if ( mod (n_couple, hydra_obj_io%timer_diter) == 0 ) call timer_report ( n_couple )
    endif

  enddo

  ! Timing summary (collective - every rank contributes its accumulators)
  call timer_summary()

  if (mpi_is_root) print*, 'FINAL POSTPROCESS'
  call FUSS%postprocess( FUSS_simulation )
  call ARES%postprocess( ARES_simulation )
  if (hydra_obj_io%coupling_wall) then
    hydrafile='hydra-wall'
    call Write_Wall_Solution( FUSS_simulation%domain(1), ARES_simulation%domain(1), hydrafile, 1 )
  endif

  ! Free persistent MPI requests before finalizing
  call cleanup_coupling_schedule()
#ifdef USE_MPI
  call FUSS_cleanup_ghost()
  call ARES_cleanup_ghost()
#endif

  ! Finalize MPI environment
  call FUSS_mpi_finalize()


contains

  subroutine Dummy_Function

  end subroutine Dummy_Function

end program hydra_program