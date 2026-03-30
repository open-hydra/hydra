module hydra_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  
  private :: read_ini

  logical :: run_ICE   = .true.
  logical :: run_IGLOO = .true.
  real(8) :: relaxCoeff
  ! logical :: adaptSwitch

  type :: hydra_type
    integer,  allocatable, dimension(:) :: N_Coupling, MOSE_iter, FUSS_iter, ICE_iter, IGLOO_iter
    real(R8), allocatable, dimension(:) :: relCoeff
    ! logical,  allocatable, dimension(:) :: adaptive
    integer              :: N_cycles
    logical              :: time_accurate
    real(R8)             :: timemax
    integer              :: coupling_diter
    real(R8)             :: coupling_dtime
    logical              :: coupling_wall
    logical              :: coupling_overwrite      
    character(len=32)    :: coupling_format(2)
    integer              :: diter_ini
    logical              :: relaxation
    integer              :: iterMax
    ! integer            :: iter_span
    integer              :: MGL
    integer, allocatable :: iter_MG(:)

    contains

      procedure, pass(self) :: ini_setup => read_ini

  end type hydra_type
  
contains

  subroutine read_ini(self)
    use Finer, only: file_ini
    use strings, only: parse

    implicit none
    class(hydra_type)  :: self
    ! Local
    type(file_ini)     :: fini
    integer            :: error, m, i, nRelax!, nAdaptive
    character(len=128) :: wholestring
    character(len=20)  :: option

    ! Load input.ini
    call fini%load(filename='input.ini')

    ! Coupling Parameters
    if (allocated(self%N_Coupling)) deallocate(self%N_Coupling)
    if (fini%has_option(option_name='N-coupling')) then
      allocate(self%N_Coupling(1:fini%count_values(section_name='hydra-parameters', option_name='N-coupling')))
    else
      allocate(self%N_Coupling(1))
    endif
    if (allocated(self%FUSS_iter)) deallocate(self%FUSS_iter)
    if (fini%has_option(option_name='FUSS-iter')) then
      allocate(self%FUSS_iter(1:fini%count_values(section_name='hydra-parameters', option_name='FUSS-iter')))
    else
      allocate(self%FUSS_iter(1))
    endif
    if (allocated(self%MOSE_iter)) deallocate(self%MOSE_iter)
    if (fini%has_option(option_name='MOSE-iter')) then
      allocate(self%MOSE_iter(1:fini%count_values(section_name='hydra-parameters', option_name='MOSE-iter')))
    else
      allocate(self%MOSE_iter(1))
    endif
    if (allocated(self%ICE_iter)) deallocate(self%ICE_iter)
    if (fini%has_option(option_name='ICE-iter')) then
      allocate(self%ICE_iter(1:fini%count_values(section_name='hydra-parameters', option_name='ICE-iter')))
    else
      allocate(self%ICE_iter(1))
    endif
    if (allocated(self%IGLOO_iter)) deallocate(self%IGLOO_iter)
    if (fini%has_option(option_name='IGLOO-iter')) then
      allocate(self%IGLOO_iter(1:fini%count_values(section_name='hydra-parameters', option_name='IGLOO-iter')))
    else
      allocate(self%IGLOO_iter(1))
    endif

    call fini%get(section_name='hydra-parameters', option_name='N-coupling', val=self%N_Coupling, error=error)
    if (error/=0) self%N_Coupling = huge(1)
    call fini%get(section_name='hydra-parameters', option_name='FUSS-iter', val=self%FUSS_iter, error=error)
    if (error/=0) self%FUSS_iter = 1
    call fini%get(section_name='hydra-parameters', option_name='MOSE-iter', val=self%MOSE_iter, error=error)
    if (error/=0) self%MOSE_iter = 1

    call fini%get(section_name='hydra-parameters', option_name='ICE-iter', val=self%ICE_iter, error=error)
    if (error/=0) then
      self%ICE_iter = 0
      run_ICE = .false.
    endif

    call fini%get(section_name='hydra-parameters', option_name='IGLOO-iter', val=self%IGLOO_iter, error=error)
    if (error/=0) then
      self%IGLOO_iter = 0
      run_IGLOO = .false.
    else
      do i=1, size(self%IGLOO_iter)
        if     (self%IGLOO_iter(i) < 1) then
          self%IGLOO_iter(i) = 0
        elseif (self%IGLOO_iter(i) > 1) then
          self%IGLOO_iter(i) = 2
        endif
      enddo
      call fini%get(section_name='hydra-parameters', option_name='relaxation', val=self%relaxation, error=error)
      if (error/=0) self%relaxation = .false.
      if (allocated(self%relCoeff)) deallocate(self%relCoeff)
        allocate(self%relCoeff(1:size(self%IGLOO_iter)))
      if (self%relaxation) then
        if (fini%has_option(option_name='relaxation-coeff')) then
          nRelax = fini%count_values(section_name='hydra-parameters', option_name='relaxation-coeff')
          if (nRelax/=size(self%IGLOO_iter) .and. nRelax/=1) error stop ( 'Non-constant relaxation coefficients must be defined for each cycle!' )
          if (nRelax==1) then
            call fini%get(section_name='hydra-parameters', option_name='relaxation-coeff', val=relaxCoeff, error=error)
            if (error/=0) then
              self%relCoeff = 5.d-1
            else
              self%relCoeff = relaxCoeff
            endif
          else 
            call fini%get(section_name='hydra-parameters', option_name='relaxation-coeff', val=self%relCoeff, error=error)
            if (error/=0) self%relCoeff = 5.d-1
          endif
        else
          self%relCoeff = 5.d-1
        endif
        write(*,*) ''
        write(*,*) '-  Under-relaxation coefficients --> ', self%relCoeff
      else
        self%relCoeff = 1.d0
      endif
    endif

    call fini%get(section_name='hydra-parameters', option_name='time-accurate', val=self%time_accurate, error=error)
    if (error/=0) self%time_accurate = .false.
    if (self%time_accurate ) then
      call fini%get(section_name='hydra-parameters', option_name='time-threshold', val=self%timemax, error=error)
      if (error/=0) self%timemax = huge(1d0)
    endif

    !>>> ADAPTIVE CYCLES
    ! if (allocated(self%adaptive)) deallocate(self%adaptive)
    ! allocate(self%adaptive(1:size(self%MOSE_iter)))
    ! if (fini%has_option(option_name='adaptive-cycles')) then
    !   nAdaptive = fini%count_values(section_name='hydra-parameters', option_name='adaptive-cycles')
    !   if (nAdaptive/=size(self%MOSE_iter) .and. nAdaptive/=1) error stop ( 'MOSE iteration adaptivity must be decided for each cycle!' )
    !   if (nAdaptive==1) then
    !     call fini%get(section_name='hydra-parameters', option_name='adaptive-cycles', val=adaptSwitch, error=error)
    !     if (error/=0) then
    !       self%adaptive = .false.
    !     else
    !       self%adaptive = adaptSwitch
    !     endif
    !   else 
    !     call fini%get(section_name='hydra-parameters', option_name='adaptive-cycles', val=self%adaptive, error=error)
    !     if (error/=0) self%adaptive = .false.
    !   endif
    !   call fini%get(section_name='hydra-parameters', option_name='iter-span', val=self%iter_span, error=error)
    !   if (error/=0) self%iter_span = 100
    !   call fini%get(section_name='hydra-parameters', option_name='iter-max', val=self%iterMax, error=error)
    !   if (error/=0) self%iterMax = 100000
    !   write(*,*) '-  Adaptive-cycles number        --> ', pack([(i, i=1, size(self%adaptive))], self%adaptive)
    !   write(*,*) '' 
    ! else
    !   self%adaptive = .false.  
    ! endif

    ! Multigrid
    call fini%get(section_name='hydra-multigrid', option_name='levels', val=self%MGL, error=error)
    if (error/=0) self%MGL = 1
    if (allocated(self%iter_MG)) deallocate(self%iter_MG)
    allocate(self%iter_MG(self%MGL))

    do m = 1, self%MGL
      write (option,'(A5,I0,A5)') 'level', m, '-iter'
      call fini%get(section_name='hydra-multigrid', option_name=option, &
                    val=self%iter_MG(m), error=error)
      if ( error/=0 .and. m>1 ) self%iter_MG(m) = 0 ! for grid1 itermax from -parameters section
    enddo

    ! IO
    self%coupling_diter = huge(1)
    self%coupling_dtime = huge(1.d0)
    self%coupling_wall = .false.
    self%coupling_overwrite = .true.
    self%coupling_format(1) = 'tecplot'
    self%coupling_format(2) = 'ascii'

    call fini%get(section_name='hydra-IO', option_name='coupling-diter', val=self%coupling_diter, error=error)
    call fini%get(section_name='hydra-IO', option_name='coupling-dtime', val=self%coupling_dtime, error=error)
    call fini%get(section_name='hydra-IO', option_name='coupling-wall', val=self%coupling_wall, error=error)
    if ( self%coupling_wall ) call fini%get(section_name='hydra-IO', option_name='coupling-overwrite', val=self%coupling_overwrite, error=error)
    call fini%get(section_name='hydra-IO', option_name='coupling-format', val=wholestring, error=error)
    if (error==0) call parse(wholestring,' ', self%coupling_format)

    ! Shell
    self%diter_ini = 1000
    call fini%get(section_name='hydra-IO', option_name='ini-diter', val=self%diter_ini, error=error)

  end subroutine read_ini
  
  
end module hydra_Types_m
