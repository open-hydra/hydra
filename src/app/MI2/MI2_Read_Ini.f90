module hydraMI2_Read_Ini
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Read_Inifile_MI2, Read_Inifile_MI2_Runtime

contains

  !> MI2's input reader: the common registrations plus MI2's own, populating hydra_obj_* and
  !  hydra_MI2_obj from input.ini.
  subroutine Read_Inifile_MI2()
    use Finer,                only: file_ini
    use hydra_Read_Sim_Param, only: Register_Sim_Param
    use hydra_Read_IO,        only: Register_IO_Fields
    use hydra_Read_Numerics,  only: Register_Numerics
    use hydra_Backend_INI,    only: Load_Ini, Scan_Ini
    use hydra_Input_Registry
    use hydraMI2_Read_Params, only: Register_MI2_Params
    implicit none
    type(file_ini) :: fini
    integer :: nmgl

    call fini%load(filename='input.ini')

    ! Size the common arrays and read the multigrid depth
    call Scan_Ini(fini, nmgl)
    ! Size the MI2 arrays before registration
    call Scan_Ini_MI2(fini)

    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Numerics(nmgl)
    call Register_MI2_Params()
    call Release_Unused_Keys()

    call Load_Ini(fini)

    call Normalize_MI2_Params()

  end subroutine Read_Inifile_MI2


  !> Hot re-read of input.ini: refreshes the values of the registered keys (sizes are not
  !  re-scanned).
  subroutine Read_Inifile_MI2_Runtime()
    use Finer,             only: file_ini
    use hydra_Backend_INI, only: Load_Ini
    use hydra_Input_Registry
    implicit none
    type(file_ini) :: fini
    character(len=1024) :: out

    call fini%load(filename='input.ini')
    call Load_Ini(fini)
    out = Validate_Registry()
    if (trim(out) /= '') then
      write(*,'(A)') trim(out)
      error stop 1
    endif

    call Normalize_MI2_Params()

  end subroutine Read_Inifile_MI2_Runtime


  !> The shared keys of solvers MI2 does not run (FUSS, the multigrid levels) are not required.
  subroutine Release_Unused_Keys()
    use hydra_Input_Registry, only: reg
    implicit none
    integer :: i

    do i = 1, reg%size
      if (reg%params(i)%name == 'FUSS-iter' .or. reg%params(i)%section == 'HYDRA-Multigrid') &
        reg%params(i)%required = .false.
    enddo

  end subroutine Release_Unused_Keys


  !> Allocate the MI2 arrays at the length the file carries and record whether each key is present.
  subroutine Scan_Ini_MI2(fini)
    use Finer, only: file_ini
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydra_Backend_INI,       only: Count_Values_Checked
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: n

    associate(o => hydra_MI2_obj)

    ! --- ICE-iter ---
    o%run_ICE = fini%has_option(option_name='ICE-iter')
    if (allocated(o%ICE_iter)) deallocate(o%ICE_iter)
    if (o%run_ICE) then
      n = Count_Values_Checked(fini,'HYDRA-Parameters','ICE-iter')
      allocate(o%ICE_iter(1:max(n,1)))
    else
      allocate(o%ICE_iter(1))
    endif

    ! --- IGLOO-iter ---
    o%run_IGLOO = fini%has_option(option_name='IGLOO-iter')
    if (allocated(o%IGLOO_iter)) deallocate(o%IGLOO_iter)
    if (o%run_IGLOO) then
      n = Count_Values_Checked(fini,'HYDRA-Parameters','IGLOO-iter')
      allocate(o%IGLOO_iter(1:max(n,1)))
    else
      allocate(o%IGLOO_iter(1))
    endif

    ! --- relaxation-coeff: one value, or one per cycle ---
    if (fini%has_option(option_name='relaxation-coeff')) then
      o%n_relax_given = Count_Values_Checked(fini,'HYDRA-Parameters','relaxation-coeff')
    else
      o%n_relax_given = 0
    endif
    if (allocated(o%relCoeff_raw)) deallocate(o%relCoeff_raw)
    allocate(o%relCoeff_raw(1:max(o%n_relax_given,1)))

    end associate

  end subroutine Scan_Ini_MI2


  !> Cross-key normalisation and validation after the read: IGLOO mode clamp, relaxation
  !  coefficients per cycle, cycle-length agreement with MOSE-iter.
  subroutine Normalize_MI2_Params()
    use hydra_Config_Types_m,    only: hydra_obj_sim_param
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    integer :: i, nMOSE

    associate(o => hydra_MI2_obj)

    ! --- IGLOO mode clamp: anything <1 is OFF, anything >1 is Eulerian coupling ---
    if (o%run_IGLOO) then
      do i = 1, size(o%IGLOO_iter)
        if     (o%IGLOO_iter(i) < 1) then; o%IGLOO_iter(i) = 0
        elseif (o%IGLOO_iter(i) > 1) then; o%IGLOO_iter(i) = 2
        endif
      enddo
    else
      o%IGLOO_iter = 0
    endif

    if (.not. o%run_ICE) o%ICE_iter = 0

    ! --- relaxation coefficients, one per cycle ---
    if (allocated(o%relCoeff)) deallocate(o%relCoeff)
    allocate(o%relCoeff(1:size(o%IGLOO_iter)))

    if (o%relaxation) then
      if (o%n_relax_given <= 1) then
        ! absent or a single scalar -> broadcast
        o%relCoeff = o%relCoeff_raw(1)
      elseif (o%n_relax_given == size(o%IGLOO_iter)) then
        o%relCoeff = o%relCoeff_raw(1:size(o%IGLOO_iter))
      else
        write(*,'(A)') '[ERROR] HYDRA-Parameters: relaxation-coeff must have either 1 value or '// &
                       'one per IGLOO-iter cycle.'
        error stop 1
      endif
    else
      ! relaxation off: coefficient 1
      o%relCoeff = 1._R8
    endif

    ! --- cycle-length agreement: every active solver drives the same cycle list as MOSE ---
    nMOSE = size(hydra_obj_sim_param%MOSE_iter)
    if (o%run_ICE .and. size(o%ICE_iter) /= nMOSE) then
      write(*,'(A,I0,A,I0)') '[ERROR] HYDRA-Parameters: ICE-iter has ', size(o%ICE_iter), &
                             ' values but MOSE-iter has ', nMOSE
      error stop 1
    endif
    if (o%run_IGLOO .and. size(o%IGLOO_iter) /= nMOSE) then
      write(*,'(A,I0,A,I0)') '[ERROR] HYDRA-Parameters: IGLOO-iter has ', size(o%IGLOO_iter), &
                             ' values but MOSE-iter has ', nMOSE
      error stop 1
    endif

    end associate

  end subroutine Normalize_MI2_Params

end module hydraMI2_Read_Ini
