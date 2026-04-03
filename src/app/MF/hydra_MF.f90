program hydra_program
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use hydra_Parameters_m
  use hydra_Read_Ini,                only: Read_Inifile
  use hydra_Config_Types_m
  use hydraMF_Mod_Import,            only: Assign_Phases
  use hydraMF_IO_Wall
  use hydraMF_Lib_Metrics,           only: Setup_Coupled_Metrics
  use hydraMF_Lib_Ghost,             only: Fill_Coupled_Ghost_Cell
  use hydraMF_Mod_dt,                only: Fix_Global_dt
  use hydraMF_Mod_BC_Coupled_Fluxes, only: BC_Coupled_Fluxes
  use FUSS_Advanced_Types_m, only: FUSS_simulation_type
  use FUSS_Config_Types_m,   only: FUSS_obj_sim_param => obj_sim_param, FUSS_obj_multigrid => obj_multigrid
  use FUSS_Procedures_m,     only: FUSS_type
  use MOSE_Advanced_Types_m, only: MOSE_simulation_type
  use MOSE_Config_Types_m,   only: MOSE_obj_sim_param => obj_sim_param, MOSE_obj_multigrid => obj_multigrid
  use MOSE_Procedures_m,     only: MOSE_type
  use IR_precision
  implicit none
  type(MOSE_type)            :: MOSE
  type(MOSE_simulation_type) :: MOSE_simulation
  type(FUSS_type)            :: FUSS
  type(FUSS_simulation_type) :: FUSS_simulation
  ! Local
  integer :: nthreads
  character(len=clen_) :: extension
  character(len=hlen_) :: hydrafile
  integer  :: n_couple, n_wall, n, nmin, nmax
  integer  :: TODO, ex_FUSS_TODO, ex_MOSE_TODO
  integer  :: i_FUSS, i_MOSE, FUSS_iter, MOSE_iter
  integer  :: mgl, m
  real(R8) :: hydra_time, pp_time


#if defined (_OPENMP)
  !$omp parallel
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  FUSS_obj_sim_param%nthreads = nthreads
  MOSE_obj_sim_param%nthreads = nthreads
  write(*,'(A)')    ' Parallel execution'
  write(*,'(A)')    ' OpenMP:'
  write(*,'(A,I4)') ' -  Number of threads --> ', nthreads
#else
  !if (mpi_is_root) write(*,'(A)')    ' Serial execution'
  write(*,'(A)')    ' Serial execution'
  FUSS_obj_sim_param%nthreads = 1
  MOSE_obj_sim_param%nthreads = 1
#endif

  ! Populate hydra_objs
  call Read_Inifile ()

  ! Setup Time-Accurate
  hydra_time = 0.d0
  pp_time = 0.d0
  if ( hydra_obj_sim_param%time_accurate ) then
    FUSS_obj_sim_param%HYDRA_time_accurate = .true.
    MOSE_obj_sim_param%HYDRA_time_accurate = .true.
  endif

  ! Setup Multigrid
  if (hydra_obj_multigrid%MGL > 1) then
    FUSS_obj_sim_param%HYDRA_MG = .true.
    MOSE_obj_sim_param%HYDRA_MG = .true.
    FUSS_obj_multigrid%MGL = hydra_obj_multigrid%MGL
    MOSE_obj_multigrid%MGL = hydra_obj_multigrid%MGL
  endif

  ! Setup Solvers
  call Assign_Phases()
  call FUSS%setup( FUSS_simulation )
  call MOSE%setup( MOSE_simulation )

  if (hydra_obj_multigrid%MGL > 1) then
    do m = 1, hydra_obj_multigrid%MGL
      FUSS_simulation%domain(m)%itermax = hydra_obj_multigrid%iter_threshold(m)*hydra_obj_sim_param%FUSS_iter(1)
      MOSE_simulation%domain(m)%itermax = hydra_obj_multigrid%iter_threshold(m)*hydra_obj_sim_param%MOSE_iter(1)
    enddo
  endif

  ! Coupling Initializations
  if (hydra_obj_io%coupling_wall) then
    allocate(IOwall(1:hydra_obj_multigrid%MGL))
    do m = 1, hydra_obj_multigrid%MGL
      if (index(hydra_obj_io%coupling_format,'ascii')>0) then
        extension = '.tec'
      else
        extension = '.szplt'
      endif
      call Initialize_Wall_File( FUSS_simulation%domain(m), MOSE_simulation%domain(m), extension, m )
      call Setup_Coupled_Metrics( FUSS_simulation%domain(m), MOSE_simulation%domain(m) )
      call Fill_Coupled_Ghost_Cell( FUSS_simulation%domain(m), MOSE_simulation%domain(m) )
    enddo
  endif
  
  ! Solve & Postprocess
  TODO = 1
  FUSS_obj_sim_param%TODO = TODO
  MOSE_obj_sim_param%TODO = TODO
  n_couple = 0
  n_wall = 0
  mgl = hydra_obj_multigrid%MGL

  do while ( TODO <= 2 )

    ! Set iteration number for each solver
    if ( hydra_obj_sim_param%time_accurate ) then
      FUSS_iter = 1
      MOSE_iter = 1
    elseif (hydra_obj_multigrid%MGL > 1) then
      FUSS_iter = hydra_obj_sim_param%FUSS_iter(1)
      MOSE_iter = hydra_obj_sim_param%MOSE_iter(1)
    else
      nmin = 0
      nmax = hydra_obj_sim_param%N_Coupling(1)
      do n = 1, size(hydra_obj_sim_param%N_Coupling)
        if (n_couple >= nmin .and. n_couple < nmax ) then
          FUSS_iter = hydra_obj_sim_param%FUSS_iter(n)
          MOSE_iter = hydra_obj_sim_param%MOSE_iter(n)
        endif
        nmin = nmax
        nmax = hydra_obj_sim_param%N_Coupling(n+1)
      enddo
    endif

    ! Update hydra input.ini options
    if ( mod(n_couple+1, hydra_obj_io%ini_diter) == 0d0 ) then
      call Read_Inifile_Runtime()
      if ( hydra_obj_sim_param%time_accurate ) then
        FUSS_obj_sim_param%HYDRA_time_accurate = .true.
        MOSE_obj_sim_param%HYDRA_time_accurate = .true.
      else
        FUSS_obj_sim_param%HYDRA_time_accurate = .false.
        MOSE_obj_sim_param%HYDRA_time_accurate = .false.
      endif
    endif

    ! Set Global dt for time-accurate
    if ( hydra_obj_sim_param%time_accurate ) then
      if ( hydra_obj_multigrid%MGL > 1 ) then
        write(*,*) 'Multi-Grid can not be used for time-accurate simulations for now!'
        stop
      else
        call Fix_Global_dt ( FUSS_simulation, MOSE_simulation, hydra_time )
      endif
    endif

    ! Solve FUSS
    do i_FUSS = 1, FUSS_iter
      if ( FUSS_obj_sim_param%TODO <= 2 ) then
        if (hydra_obj_io%coupling_wall) call BC_Coupled_Fluxes ( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl) )
        call FUSS%solve( FUSS_simulation, Dummy_Function )
        if ( FUSS_obj_sim_param%TODO <= 2 ) call FUSS%postprocess( FUSS_simulation )
      endif
    enddo

    ! Solve MOSE
    do i_MOSE = 1, MOSE_iter
      if ( MOSE_obj_sim_param%TODO <= 2 ) then
        if (hydra_obj_io%coupling_wall) call BC_Coupled_Fluxes ( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl) )
        call MOSE%solve( MOSE_simulation, Dummy_Function )
        if ( MOSE_obj_sim_param%TODO <= 2 ) call MOSE%postprocess( MOSE_simulation )
      endif
    enddo
    
    n_couple = n_couple + 1
    
    ! Intermediate Postprocess
    if (hydra_obj_multigrid%MGL > 1 .and. n_couple == sum(hydra_obj_multigrid%iter_threshold(mgl:))) then
      print*, 'COUPLING-STEP POSTPROCESS - MG'
      if (hydra_obj_io%coupling_wall) then
        n_wall = n_wall+1
        if (mgl == 1) then
          hydrafile = 'hydra-wall'
          if (.not.hydra_obj_io%coupling_overwrite) hydrafile = trim(hydrafile)//trim(str(.true.,n_wall))
        else
          hydrafile = 'hydra-wall-level'//trim(str(.true.,mgl))
        endif
        call Write_Wall_Solution( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl), mgl, hydrafile )
      endif

      mgl = mgl - 1

    elseif ( mod(n_couple, hydra_obj_io%coupling_diter) == 0d0 .or. hydra_time >= pp_time + hydra_obj_io%coupling_dtime ) then
      print*, 'COUPLING-STEP POSTPROCESS'
      FUSS_obj_sim_param%HYDRA_postprocess = .true.
      MOSE_obj_sim_param%HYDRA_postprocess = .true.
      ex_FUSS_TODO = FUSS_obj_sim_param%TODO
      ex_MOSE_TODO = MOSE_obj_sim_param%TODO
      FUSS_obj_sim_param%TODO = 2
      MOSE_obj_sim_param%TODO = 2
      call FUSS%postprocess( FUSS_simulation )
      call MOSE%postprocess( MOSE_simulation )
      if (hydra_obj_io%coupling_wall) then
        n_wall = n_wall+1
        if (mgl == 1) then
          hydrafile = 'hydra-wall'
          if (.not.hydra_obj_io%coupling_overwrite) hydrafile = trim(hydrafile)//trim(str(.true.,n_wall))
        else
          hydrafile = 'hydra-wall-level'//trim(str(.true.,mgl))
        endif
        call Write_Wall_Solution( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl), mgl, hydrafile )
      endif
      FUSS_obj_sim_param%TODO = ex_FUSS_TODO
      MOSE_obj_sim_param%TODO = ex_MOSE_TODO
      FUSS_obj_sim_param%HYDRA_postprocess = .false.
      MOSE_obj_sim_param%HYDRA_postprocess = .false.
      pp_time = hydra_time
    endif

    ! Check if Exit from loop
    if ( FUSS_obj_sim_param%TODO == 3 .and. MOSE_obj_sim_param%TODO == 3 ) then
      TODO = 3
    endif
    if ( n_couple == hydra_obj_sim_param%N_Coupling(size(hydra_obj_sim_param%N_Coupling)) ) then
      TODO = 3
      FUSS_obj_sim_param%TODO = 3
      MOSE_obj_sim_param%TODO = 3
    endif
    if ( hydra_obj_multigrid%MGL > 1 .and. mgl == 0 ) then
      TODO = 3
      FUSS_obj_sim_param%TODO = 3
      MOSE_obj_sim_param%TODO = 3
    endif
    if ( hydra_obj_sim_param%time_accurate .and. hydra_time >= hydra_obj_sim_param%time_threshold) then
      TODO = 3
      FUSS_obj_sim_param%TODO = 3
      MOSE_obj_sim_param%TODO = 3
    endif

    if ( TODO <= 2 ) then
      FUSS_obj_sim_param%TODO = 1
      MOSE_obj_sim_param%TODO = 1
    endif

  enddo

  print*, 'FINAL POSTPROCESS'
  call FUSS%postprocess( FUSS_simulation )
  call MOSE%postprocess( MOSE_simulation )
  if (hydra_obj_io%coupling_wall) then
    hydrafile='hydra-wall'
    call Write_Wall_Solution( FUSS_simulation%domain(1), MOSE_simulation%domain(1), 1, hydrafile, hydra_obj_io%coupling_format )
  endif


contains

  subroutine Dummy_Function

  end subroutine Dummy_Function

end program hydra_program