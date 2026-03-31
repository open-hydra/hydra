program hydra_program
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use IR_precision
  use hydra_Types_m
  use hydraMF_Mod_Import
  use hydraMF_Lib_Metrics
  use hydraMF_Lib_Ghost
  use hydraMF_Mod_Coupling
  use hydraMF_IO_Wall
  use FUSS_Procedures_m
  use FUSS_Types_m
  use MOSE_Procedures_m
  use MOSE_Types_m
  
  implicit none
  type(HYDRA_type)           :: HYDRA
  type(MOSE_type)            :: MOSE
  type(MOSE_simulation_type) :: MOSE_simulation
  type(FUSS_type)            :: FUSS
  type(FUSS_simulation_type) :: FUSS_simulation
  ! Local
  integer            :: nthreads, n_couple, n_wall, n, nmin, nmax
  integer            :: TODO, ex_FUSS_TODO, ex_MOSE_TODO
  integer            :: i_FUSS, i_MOSE, FUSS_iter, MOSE_iter
  character(len=300) :: hydrafile
  real(R8)           :: hydra_time, pp_time
  integer            :: mgl, m


#if defined (_OPENMP)
  !$omp parallel
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  FUSS_simulation%config%nthreads = nthreads
  MOSE_simulation%config%nthreads = nthreads
  write(*,'(A)')   ' OpenMP:'
  write(*,'(A,I4)')' -  Number of threads --> ', nthreads
#endif

  ! Setup HYDRA
  call HYDRA%ini_setup()

  ! Setup Time-Accurate
  hydra_time = 0.d0
  pp_time = 0.d0
  if ( HYDRA%time_accurate ) then
    FUSS_simulation%config%HYDRA_time_accurate = .true.
    MOSE_simulation%config%HYDRA_time_accurate = .true.
  endif

  ! Setup Multigrid
  if (HYDRA%MGL > 1) then
    FUSS_simulation%config%HYDRA_MG = .true.
    FUSS_simulation%config%MGL = HYDRA%MGL
    MOSE_simulation%config%HYDRA_MG = .true.
    MOSE_simulation%config%MGL = HYDRA%MGL
  endif

  ! Setup Solvers
  call Assign_Phases()
  call FUSS%setup( FUSS_simulation )
  call MOSE%setup( MOSE_simulation )

  if (HYDRA%MGL > 1) then
    do m = 1, HYDRA%MGL
      FUSS_simulation%domain(m)%itermax = HYDRA%iter_MG(m)*HYDRA%FUSS_iter(1)
      MOSE_simulation%domain(m)%itermax = HYDRA%iter_MG(m)*HYDRA%MOSE_iter(1)
    enddo
  endif

  ! Coupling Initializations
  if (HYDRA%coupling_wall) then
    allocate(IOwall(1:HYDRA%MGL))
    do m = 1, HYDRA%MGL
      call Initialize_Wall_File( FUSS_simulation%domain(m), MOSE_simulation%domain(m), m )
      call Setup_Coupled_Metrics( FUSS_simulation%domain(m), MOSE_simulation%domain(m) )
      call Fill_Coupled_Ghost_Cell( FUSS_simulation%domain(m), MOSE_simulation%domain(m) )
    enddo
  endif
  
  ! Solve & Postprocess
  TODO = 1
  FUSS_simulation%config%TODO = TODO
  MOSE_simulation%config%TODO = TODO
  n_couple = 0
  n_wall = 0
  mgl = HYDRA%MGL

  do while ( TODO <= 2 )

    ! Set iteration number for each solver
    if ( HYDRA%time_accurate ) then
      FUSS_iter = 1
      MOSE_iter = 1
    elseif (HYDRA%MGL > 1) then
      FUSS_iter = HYDRA%FUSS_iter(1)
      MOSE_iter = HYDRA%MOSE_iter(1)
    else
      nmin = 0
      nmax = HYDRA%N_Coupling(1)
      do n = 1, size(HYDRA%N_Coupling)
        if (n_couple >= nmin .and. n_couple < nmax ) then
          FUSS_iter = HYDRA%FUSS_iter(n)
          MOSE_iter = HYDRA%MOSE_iter(n)
        endif
        nmin = nmax
        nmax = HYDRA%N_Coupling(n+1)
      enddo
    endif

    ! Update HYDRA input.ini options
    if ( mod(n_couple+1, HYDRA%diter_ini) == 0d0 ) then
      call HYDRA%ini_setup()
      if ( HYDRA%time_accurate ) then
        FUSS_simulation%config%HYDRA_time_accurate = .true.
        MOSE_simulation%config%HYDRA_time_accurate = .true.
      else
        FUSS_simulation%config%HYDRA_time_accurate = .false.
        MOSE_simulation%config%HYDRA_time_accurate = .false.
      endif
    endif

    ! Set Global dt for time-accurate
    if ( HYDRA%time_accurate ) then
      if ( HYDRA%MGL > 1 ) then
        write(*,*) 'Multi-Grid can not be used for time-accurate simulations for now!'
        stop
      else
        call Fix_Global_dt ( FUSS_simulation, MOSE_simulation, hydra_time )
      endif
    endif

    ! Solve FUSS
    do i_FUSS = 1, FUSS_iter
      if ( FUSS_simulation%config%TODO <= 2 ) then
        if (HYDRA%coupling_wall) call BC_Coupled_Fluxes ( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl) )
        call FUSS%solve( FUSS_simulation )
        if ( FUSS_simulation%config%TODO <= 2 ) call FUSS%postprocess( FUSS_simulation )
      endif
    enddo

    ! Solve MOSE
    do i_MOSE = 1, MOSE_iter
      if ( MOSE_simulation%config%TODO <= 2 ) then
        if (HYDRA%coupling_wall) call BC_Coupled_Fluxes ( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl) )
        call MOSE%solve( MOSE_simulation, Dummy_Function )
        if ( MOSE_simulation%config%TODO <= 2 ) call MOSE%postprocess( MOSE_simulation )
      endif
    enddo
    
    n_couple = n_couple + 1
    
    ! Intermediate Postprocess
    if (HYDRA%MGL > 1 .and. n_couple == sum(HYDRA%iter_MG(mgl:))) then
      print*, 'COUPLING-STEP POSTPROCESS - MG'
      if (HYDRA%coupling_wall) then
        n_wall = n_wall+1
        if (mgl == 1) then
          hydrafile = 'hydra-wall'
          if (.not.HYDRA%coupling_overwrite) hydrafile = trim(hydrafile)//trim(str(.true.,n_wall))
        else
          hydrafile = 'hydra-wall-level'//trim(str(.true.,mgl))
        endif
        call Write_Wall_Solution( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl), mgl, hydrafile, HYDRA%coupling_format )
      endif

      mgl = mgl - 1

    elseif ( mod(n_couple, HYDRA%coupling_diter) == 0d0 .or. hydra_time >= pp_time + HYDRA%coupling_dtime ) then
      print*, 'COUPLING-STEP POSTPROCESS'
      FUSS_simulation%config%HYDRA_postprocess = .true.
      MOSE_simulation%config%HYDRA_postprocess = .true.
      ex_FUSS_TODO = FUSS_simulation%config%TODO
      ex_MOSE_TODO = MOSE_simulation%config%TODO
      FUSS_simulation%config%TODO = 2
      MOSE_simulation%config%TODO = 2
      call FUSS%postprocess( FUSS_simulation )
      call MOSE%postprocess( MOSE_simulation )
      if (HYDRA%coupling_wall) then
        n_wall = n_wall+1
        if (mgl == 1) then
          hydrafile = 'hydra-wall'
          if (.not.HYDRA%coupling_overwrite) hydrafile = trim(hydrafile)//trim(str(.true.,n_wall))
        else
          hydrafile = 'hydra-wall-level'//trim(str(.true.,mgl))
        endif
        call Write_Wall_Solution( FUSS_simulation%domain(mgl), MOSE_simulation%domain(mgl), mgl, hydrafile, HYDRA%coupling_format )
      endif
      FUSS_simulation%config%TODO = ex_FUSS_TODO
      MOSE_simulation%config%TODO = ex_MOSE_TODO
      FUSS_simulation%config%HYDRA_postprocess = .false.
      MOSE_simulation%config%HYDRA_postprocess = .false.
      pp_time = hydra_time
    endif

    ! Check if Exit from loop
    if ( FUSS_simulation%config%TODO == 3 .and. MOSE_simulation%config%TODO == 3 ) then
      TODO = 3
    endif
    if ( n_couple == HYDRA%N_Coupling(size(HYDRA%N_Coupling)) ) then
      TODO = 3
      FUSS_simulation%config%TODO = 3
      MOSE_simulation%config%TODO = 3
    endif
    if ( HYDRA%MGL > 1 .and. mgl == 0 ) then
      TODO = 3
      FUSS_simulation%config%TODO = 3
      MOSE_simulation%config%TODO = 3
    endif
    if ( HYDRA%time_accurate .and. hydra_time >= HYDRA%timemax) then
      TODO = 3
      FUSS_simulation%config%TODO = 3
      MOSE_simulation%config%TODO = 3
    endif

    if ( TODO <= 2 ) then
      FUSS_simulation%config%TODO = 1
      MOSE_simulation%config%TODO = 1
    endif

  enddo

  print*, 'FINAL POSTPROCESS'
  call FUSS%postprocess( FUSS_simulation )
  call MOSE%postprocess( MOSE_simulation )
  if (HYDRA%coupling_wall) then
    hydrafile='hydra-wall'
    call Write_Wall_Solution( FUSS_simulation%domain(1), MOSE_simulation%domain(1), 1, hydrafile, HYDRA%coupling_format )
  endif


contains

  subroutine Dummy_Function

  end subroutine Dummy_Function

end program hydra_program
