module hydraMI2_Mod_Import
  implicit none

contains

  subroutine Assign_Phases()
    use strings, only: parse
    use hydra_Misc, only: Get_File_List
    use ICE_Mod_Parameters, only: ICE_phase_prefix
    use IGLOO_variables, only: IGLOO_phase_prefix
    use MOSE_Global_m, only: MOSE_phase_prefix
    implicit none
    character(len=128)              :: filename, stringa(2), substringa(2), type
    integer                         :: i, num_files, u, ios
    character(len=128), allocatable :: file_list(:)

    ! Get the list of .txt files in the current directory
    call Get_File_List(file_list, num_files)

    do i = 1, num_files
    filename = trim(adjustl(file_list(i)))
    ! Check if the phase has a name. If not, an empty string is returned
    call parse(filename,'-',stringa)
    if (stringa(2)=='') then
        stringa(2) = stringa(1)
        stringa(1) = ''
    endif
    if (index(stringa(2),'phase')>0) then
        call parse(stringa(1),'/',substringa)
        open(newunit=u, file=filename, status="old")
        read(u,'(A)',iostat=ios) type
        if (index(type,'condensed-dispersed')>0) then
          ICE_phase_prefix = trim(substringa(2))//'-'
          IGLOO_phase_prefix = trim(substringa(2))//'-'
        elseif (index(type,'ideal-gas')>0) then
          MOSE_phase_prefix = trim(substringa(2))//'-'
        else
          write(*,*) 'Error: unknown phase type'
          stop
        endif
        close(u)
    endif
    end do
    deallocate(file_list)

  end subroutine Assign_Phases

end module hydraMI2_Mod_Import