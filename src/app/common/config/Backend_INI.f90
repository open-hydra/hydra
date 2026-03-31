module hydra_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini

contains

  subroutine Load_Ini(fini)
    use hydra_Input_Registry, only: reg
    use Finer, only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error

    do i = 1, reg%size
      error = 1

      ! Each registry entry has exactly one associated typed pointer.
      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
      else if (associated(reg%params(i)%value%iarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%iarr, error=error)
      else if (associated(reg%params(i)%value%rarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%rarr, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.
    end do

  end subroutine


  subroutine Scan_Ini(fini, nnmgl)
    use Finer, only: file_ini
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(out) :: nmgl
    ! Local
    integer :: error
    

    ! Check for N-coupling, FUSS-iter, MOSE-iter in the input.ini and allocate arrays accordingly
    if (allocated(obj_hydra_sim_param%N_Coupling)) deallocate(obj_hydra_sim_param%N_Coupling)
    if (fini%has_option(option_name='N-coupling')) then
      allocate(obj_hydra_sim_param%N_Coupling(1:fini%count_values(section_name='HYDRA-parameters', option_name='N-coupling')))
    else
      allocate(obj_hydra_sim_param%N_Coupling(1))
    endif
    if (allocated(obj_hydra_sim_param%FUSS_iter)) deallocate(obj_hydra_sim_param%FUSS_iter)
    if (fini%has_option(option_name='FUSS-iter')) then
      allocate(obj_hydra_sim_param%FUSS_iter(1:fini%count_values(section_name='HYDRA-parameters', option_name='FUSS-iter')))
    else
      allocate(obj_hydra_sim_param%FUSS_iter(1))
    endif
    if (allocated(obj_hydra_sim_param%MOSE_iter)) deallocate(obj_hydra_sim_param%MOSE_iter)
    if (fini%has_option(option_name='MOSE-iter')) then
      allocate(obj_hydra_sim_param%MOSE_iter(1:fini%count_values(section_name='HYDRA-parameters', option_name='MOSE-iter')))
    else
      allocate(obj_hydra_sim_param%MOSE_iter(1))
    endif

    ! Read the MGL
    nmgl = 1
    call fini%get(section_name='HYDRA-multigrid', option_name='levels', val=nmgl, error=error)

  end subroutine Scan_Ini

end module hydra_Backend_INI