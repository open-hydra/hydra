module hydra_Read_Sim_Param
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: Register_Sim_Param


contains


  subroutine Register_Sim_Param()
    use hydra_Input_Registry
    use hydra_Config_Types_m, only: hydra_obj_sim_param
    implicit none
    ! Local
    character(len=:), allocatable :: section

    hydra_obj_sim_param%warning_message = 'none'
    hydra_obj_sim_param%error_message   = 'none'
    hydra_obj_sim_param%description     = 'none'

    section = 'HYDRA-Parameters'

    call reg%add(section,'N-coupling',hydra_obj_sim_param%N_Coupling, '1000000000', 'Number of coupling iterations', 'integer array', .true. )
    call reg%add(section,'FUSS-iter',hydra_obj_sim_param%FUSS_iter, '1', 'Number of FUSS iterations per coupling iteration', 'integer array', .true. )
    call reg%add(section,'MOSE-iter',hydra_obj_sim_param%MOSE_iter, '1', 'Number of MOSE iterations per coupling iteration', 'integer array', .true. )
    call reg%add(section,'time-accurate',hydra_obj_sim_param%time_accurate, '.false.', 'Time accurate switch', 'logical', .true. )
    call reg%add(section,'time-threshold',hydra_obj_sim_param%time_threshold,'1e30','Maximum simulation time','> 0',.false.)
    call reg%add(section,'iter-threshold',hydra_obj_sim_param%iter_threshold,'1000000000','Maximum number of iterations','> 0',.false.)

  end subroutine Register_Sim_Param


end module hydra_Read_Sim_Param