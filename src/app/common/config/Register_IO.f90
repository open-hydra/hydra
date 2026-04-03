module hydra_Read_IO
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: Register_IO_Fields


contains


  subroutine Register_IO_Fields()
    use hydra_Input_Registry
    use hydra_Config_Types_m, only: hydra_obj_io
    implicit none
    ! Local
    character(len=:), allocatable :: section

    hydra_obj_io%warning_message = 'none'
    hydra_obj_io%error_message   = 'none'
    hydra_obj_io%description     = 'none'

    section = 'HYDRA-IO'
    
    ! IO options
    call reg%add( section, 'coupling-format', hydra_obj_io%coupling_format, 'tecplot ascii', 'Solution (OUTPUT/field.*) format', 'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )
    call reg%add( section, 'coupling-diter', hydra_obj_io%coupling_diter, '1000000000', 'Coupling solution output iter frequency', '> 0', .false. )
    call reg%add( section, 'coupling-dtime', hydra_obj_io%coupling_dtime, '1e30', 'Coupling solution output time frequency', '> 0', .false. )
    call reg%add( section, 'coupling-overwrite', hydra_obj_io%coupling_overwrite, 'true', 'Overwrite coupling solution files', 'true, false', .false. )
    call reg%add( section, 'coupling-wall', hydra_obj_io%coupling_wall, 'false', 'Use coupling wall boundary conditions', 'true, false', .false. )

    ! Shell options
    call reg%add( section, 'ini-diter', hydra_obj_io%ini_diter, '1000', 'input.ini update iter frequency', '> 0', .false. )

  end subroutine Register_IO_Fields


end module hydra_Read_IO
