program hydraMI2_program
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use hydraMI2_Mod_Coupling
  use hydra_Types_m
  use hydraMI2_Mod_Import
  use MOSE_Procedures_m
  use MOSE_Types_m
  use ICE_Procedures_m
  use ICE_Mod_Type
  use IGLOO_module,  only: obj_IGLOO

  implicit none
  type(HYDRA_type)            :: HYDRA
  type(MOSE_type)             :: MOSE
  type(MOSE_simulation_type)  :: MOSE_simulation
  type(ICE_type)              :: ICE
  type(ICE_simulation_type)   :: ICE_simulation
  type(obj_IGLOO)             :: IGLOO
  type(resField), allocatable :: oldResidual(:)
  ! real(R8),       allocatable :: restotSpan(:)
  ! real(R8) :: oldResAve, newResAve
  integer  :: nthreads, n_couple, n_cycle
  integer  :: TODO
  integer  :: i_MOSE, i_ICE, i_IGLOO, i
  logical  :: eulCoupling, outCheck           


#if defined (_OPENMP)
  !$omp parallel
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  MOSE_simulation%config%nthreads = nthreads
  ICE_simulation%config%nthreads  = nthreads
  write(*,'(A)')   ' OpenMP:'
  write(*,'(A,I4)')' -  Number of threads --> ', nthreads
#endif

  ! Setup HYDRA
  call HYDRA%ini_setup()
  if (run_ICE   .and. (size(HYDRA%MOSE_iter) /= size(HYDRA%ICE_iter  ))) error stop ( 'HYDRA: MOSE and ICE iterations must have the same length  ' )
  if (run_IGLOO .and. (size(HYDRA%MOSE_iter) /= size(HYDRA%IGLOO_iter))) error stop ( 'HYDRA: MOSE and IGLOO iterations must have the same length' )
  HYDRA%N_cycles = size(HYDRA%MOSE_iter)

  ! Setup Solvers
  call Assign_Phases()
  ICE_simulation%config%twcoupled = .true.
  call MOSE%setup( MOSE_simulation )
  if (run_ICE  ) call ICE%setup  ( ICE_simulation, MOSE_simulation%IOfield(1) )
  if (run_IGLOO) call IGLOO%setup( MOSE_simulation%IOfield(1) )

  ! Solve & Postprocess
  n_cycle  = 0
  eulCoupling = .false.
  ! if (any(HYDRA%adaptive(:))) then
  !   allocate(restotSpan(1:HYDRA%iter_span))
  ! endif
  if ((.not.allocated(oldResidual)).and.HYDRA%relaxation) then
    call initializeResBlock(MOSE_simulation, oldResidual)
  endif

  do while ( n_cycle < HYDRA%N_cycles )
    n_cycle = n_cycle + 1
    write(*,*) ''
    write(*,'(A,I4)')' HYDRA: CYCLE NUMBER ==>', n_cycle
    
    n_couple = 0
    TODO = 1
    MOSE_simulation%config%TODO = TODO
    ICE_simulation%config%TODO = TODO
    
    do while ( TODO <= 2 )

      n_couple = n_couple+1

      write(*,*) ''
      write(*,'(A,I4)')' HYDRA: Coupling number --> ', n_couple
      write(*,*) ''

      ! Update HYDRA input.ini options
      if ( mod(n_couple, HYDRA%diter_ini) == 0d0 ) then
        call HYDRA%ini_setup()
      endif

      !! MOSE
      ! Solve MOSE
      ! if (HYDRA%adaptive(n_cycle)) then
      !   if (HYDRA%MOSE_iter(n_cycle) > 0) HYDRA%MOSE_iter(n_cycle) = 1
      !   restotSpan(:) = 1._R8
      !   oldResAve = 0._R8
      !   newResAve = 0._R8
      ! endif
      i = 0
      do while (i<=HYDRA%iterMax)
        do i_MOSE = 1, HYDRA%MOSE_iter(n_cycle)
          if ( MOSE_simulation%config%TODO <= 2 ) then
            call MOSE%solve( MOSE_simulation, MOSE_external_routine )
            if ( MOSE_simulation%config%TODO <= 2 ) call MOSE%postprocess( MOSE_simulation )
          endif
        enddo
        ! if (.not.HYDRA%adaptive(n_cycle)) exit
        ! restotSpan(mod(i,HYDRA%iter_span)+1) = sum(MOSE_simulation % config % residuotot(:))
        ! if (isAdapted(restotSpan, HYDRA%iter_span, oldResAve, newResAve)) then
        !   write(*,*) 'HYDRA: residuals reached asymptotic behavior'
        !   exit
        ! endif
        ! oldResAve = newResAve 
        i = i + 1
      enddo

      !! ICE
      if (run_ICE) then
        ! Update ICE gas phase
        call MOSExICE( MOSE_simulation%domain(1), ICE_simulation%grid )
        ! Solve ICE
        do i_ICE = 1, HYDRA%ICE_iter(n_cycle)
          if ( ICE_simulation%config%TODO <= 2 ) then
            call ICE%solve( ICE_simulation )
            if ( ICE_simulation%config%TODO <= 2 ) call ICE%postprocess( ICE_simulation )
          endif
        enddo
      endif

      !! IGLOO
      if (run_IGLOO .and. HYDRA%IGLOO_iter(n_cycle)/=0) then
        ! Solve IGLOO
        if ( HYDRA%IGLOO_iter(n_cycle) == 2) then
          eulCoupling     = .true.
          IGLOO%srcSwitch = .false.
          IGLOO%eulSwitch = .true.
          write(*,*) ' >> Eulerian coupling'
        else
          eulCoupling     = .false.
          IGLOO%eulSwitch = .false.
          IGLOO%srcSwitch = .true.
        endif
        outCheck = .true.
        do i = n_cycle+1, HYDRA%N_cycles
          if (HYDRA%IGLOO_iter(i)/=0) then
            outCheck = .false.
            exit
          endif
        enddo
        if (outCheck) then; IGLOO%eulSwitch = .true.; IGLOO%srcSwitch = .true.; endif
        ! Update IGLOO gas phase
        call MOSExIGLOO( MOSE_simulation%domain(1), IGLOO )
        call IGLOO%solve()
      endif

      ! Check if Exit from loop
      if ( n_couple == HYDRA%N_Coupling(n_cycle) ) then
        TODO = 3
      endif

    enddo

  enddo

  call MOSE%postprocess( MOSE_simulation )
  if (run_ICE)   call ICE%postprocess( ICE_simulation )
  if (run_IGLOO) call IGLOO%writeout()

  write(*,*)
  write(*,'(A)') 'GAME OVER'
  write(*,*)


contains

  subroutine MOSE_external_routine(MOSE_domain)
    implicit none
    type(MOSE_domain_type) :: MOSE_domain

    if (run_ICE)   call ICExMOSE(MOSE_domain,ICE_simulation%grid)
    if (run_IGLOO) call IGLOOxMOSE(MOSE_domain, IGLOO, eulCoupling, HYDRA%relaxation, i_MOSE,      &
                                   HYDRA%MOSE_iter(n_cycle), HYDRA%relCoeff(n_cycle), oldResidual)
    
  end subroutine MOSE_external_routine


  ! logical function isAdapted(residuototSpan, span, oldAverage, newAverage)
  !   implicit none
  !   real(R8), intent(in)  :: residuototSpan(:), oldAverage
  !   integer,  intent(in)  :: span
  !   real(R8), intent(out) :: newAverage
  !   real(R8) :: maxOscill, logRes(span), resRange
    
  !   logRes = log10(residuototSpan(:))
  !   newAverage = sum(logRes) / real(span)
  !   resRange  = maxval(logRes)-minval(logRes)
  !   maxOscill = maxval(abs(logRes-newAverage))
  !   isAdapted = (abs(oldAverage-newAverage) < 1.d-3/abs(oldAverage)/span &
  !                 .and. abs(resRange-2*maxOscill) < 1.d-3)

  ! end function isAdapted

end program hydraMI2_program
