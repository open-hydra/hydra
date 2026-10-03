module hydra_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini, Count_Values_Checked

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


  !> Number of values of an option that must be in the named section. FiNeR's count_values leaves
  !  its result unset when the section or the option is missing; section names are case-sensitive.
  integer function Count_Values_Checked(fini, section, option) result(n)
    use Finer, only: file_ini
    implicit none
    type(file_ini),   intent(in) :: fini
    character(len=*), intent(in) :: section, option

    n = 0
    if (fini%index(section_name=section, option_name=option) > 0) &
      n = fini%count_values(section_name=section, option_name=option)
    if (n < 1) then
      write(*,'(A)') '[ERROR] input.ini: "'//option//'" is not an option of section ['//section// &
                     '], or it holds no value.'
      write(*,'(A)') '        Section names are CASE-SENSITIVE: [hydra-parameters] is NOT ['//section//'].'
      error stop 1
    endif

  end function Count_Values_Checked


  subroutine Scan_Ini(fini, nmgl)
    use hydra_Config_Types_m, only: hydra_obj_sim_param
    use Finer, only: file_ini
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(out) :: nmgl
    ! Local
    integer :: error
    

    ! Check for N-coupling, FUSS-iter, MOSE-iter in the input.ini and allocate arrays accordingly
    if (allocated(hydra_obj_sim_param%N_Coupling)) deallocate(hydra_obj_sim_param%N_Coupling)
    if (fini%has_option(option_name='N-coupling')) then
      allocate(hydra_obj_sim_param%N_Coupling(1:Count_Values_Checked(fini,'HYDRA-Parameters','N-coupling')))
    else
      allocate(hydra_obj_sim_param%N_Coupling(1))
    endif
    if (allocated(hydra_obj_sim_param%FUSS_iter)) deallocate(hydra_obj_sim_param%FUSS_iter)
    if (fini%has_option(option_name='FUSS-iter')) then
      allocate(hydra_obj_sim_param%FUSS_iter(1:Count_Values_Checked(fini,'HYDRA-Parameters','FUSS-iter')))
    else
      allocate(hydra_obj_sim_param%FUSS_iter(1))
    endif
    if (allocated(hydra_obj_sim_param%MOSE_iter)) deallocate(hydra_obj_sim_param%MOSE_iter)
    if (fini%has_option(option_name='MOSE-iter')) then
      allocate(hydra_obj_sim_param%MOSE_iter(1:Count_Values_Checked(fini,'HYDRA-Parameters','MOSE-iter')))
    else
      allocate(hydra_obj_sim_param%MOSE_iter(1))
    endif
    if (allocated(hydra_obj_sim_param%ARES_iter)) deallocate(hydra_obj_sim_param%ARES_iter)
    if (fini%has_option(option_name='ARES-iter')) then
      allocate(hydra_obj_sim_param%ARES_iter(1:Count_Values_Checked(fini,'HYDRA-Parameters','ARES-iter')))
    else
      allocate(hydra_obj_sim_param%ARES_iter(1))
    endif

    ! Read the MGL
    nmgl = 1
    call fini%get(section_name='HYDRA-Multigrid', option_name='levels', val=nmgl, error=error)

  end subroutine Scan_Ini

end module hydra_Backend_INI
