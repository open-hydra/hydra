module hydraMF_Mod_Import
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Phases

contains

  subroutine Assign_Phases()
    use hydra_Parameters_m
    use hydra_Get_File_List_m, only: Get_File_List
    use FUSS_Global_m, only: FUSS_phase_prefix
    use MOSE_Global_m, only: MOSE_phase_prefix
    use strings, only: parse
    implicit none
    character(len=llen_)              :: filename, stringa(2), substringa(2), type
    integer                           :: i, num_files, u, ios
    character(len=llen_), allocatable :: file_list(:)

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
        if (index(type,'solid')>0) then
          FUSS_phase_prefix = trim(substringa(2))//'-'
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

end module hydraMF_Mod_Import