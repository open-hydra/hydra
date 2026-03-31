module hydra_Get_File_List_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: Get_File_List

contains

  subroutine Get_File_List(file_list, num_files)
      implicit none
      character(len=*), allocatable, intent(out) :: file_list(:)
      integer, intent(out) :: num_files
      character(len=256) :: line
      integer :: ios, unit, u

      open(newunit=u, file="filelist.txt", status="old", action="read", iostat=ios)

      if (ios /= 0) then
        write(*,*) "Error opening filelist.txt"
        stop
      end if

      num_files = 0
      do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        num_files = num_files + 1
      end do

      rewind(u)
      allocate(file_list(num_files))

      do unit = 1, num_files
        read(u, '(A)', iostat=ios) file_list(unit)
      end do

      close(u)
      
  end subroutine Get_File_List

end module hydra_Get_File_List_m