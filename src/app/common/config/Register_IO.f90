module hydra_Read_IO
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: Register_IO_Fields


contains


  subroutine Register_IO_Fields()
    use hydra_Input_Registry
    use hydra_Config_Types_m, only: obj_hydra_io
    implicit none
    ! Local
    character(len=:), allocatable :: section

    obj_hydra_io%warning_message = 'none'
    obj_hydra_io%error_message   = 'none'
    obj_hydra_io%description     = 'none'

    section = 'HYDRA-IO'
    
    ! IO options
    call reg%add( section, 'sol-format', obj_hydra_io%sol_format, 'tecplot ascii', 'Solution (OUTPUT/field.*) format', 'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )
    call reg%add( section, 'coupling-diter', obj_hydra_io%coupling_diter, '1000000000', 'Coupling solution output iter frequency', '> 0', .false. )
    call reg%add( section, 'coupling-dtime', obj_hydra_io%coupling_dtime, '1e30', 'Coupling solution output time frequency', '> 0', .false. )
    call reg%add( section, 'coupling-overwrite', obj_hydra_io%coupling_overwrite, 'true', 'Overwrite coupling solution files', 'true, false', .false. )
    call reg%add( section, 'coupling-wall', obj_hydra_io%coupling_wall, 'false', 'Use coupling wall boundary conditions', 'true, false', .false. )

    ! Shell options
    call reg%add( section, 'ini-diter', obj_hydra_io%ini_diter, '1000', 'input.ini update iter frequency', '> 0', .false. )

  end subroutine Register_IO_Fields


end module hydra_Read_IO