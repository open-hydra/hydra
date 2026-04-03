module hydra_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Parameters_m

  implicit none
  private

  !! ------------------------------------------------------
  !! Simulation Parameters --------------------------------
  !! ------------------------------------------------------
  type :: hydra_simulation_parameters_t
    character(len=llen_) :: warning_message
    character(len=llen_) :: error_message
    character(len=llen_) :: description
    ! USER-DEFINED INPUTS
    integer,  allocatable, dimension(:) :: N_Coupling, MOSE_iter, FUSS_iter
    logical   :: time_accurate
    real(R8)  :: time_threshold   ! Tempo max per arresto esecuzione 
    integer   :: iter_threshold   ! Numero max iterate per arresto esecuzione
  end type hydra_simulation_parameters_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Input-Output -----------------------------------------
  !! ------------------------------------------------------
  type :: hydra_io_t
    character(len=llen_) :: warning_message
    character(len=llen_) :: error_message
    character(len=llen_) :: description
    ! USER-DEFINED INPUTS
    integer              :: coupling_diter         ! couplingiteration interval to save the solution
    real(R8)             :: coupling_dtime         ! coupling time interval to save the solution
    logical              :: coupling_overwrite     ! switch to overwrite the solution
    character(len=llen_) :: coupling_format        ! solution format (native,tecplot,vtk) and formatting (ascii,raw)
    logical              :: coupling_wall
    integer              :: ini_diter              ! input.ini update
    ! Useful variables
    character(len=hlen_) :: Owallnames
    integer              :: Onwall
  end type hydra_io_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Multigrid --------------------------------------------
  !! ------------------------------------------------------
  type :: hydra_multigrid_t
    character(len=llen_) :: warning_message
    character(len=llen_) :: error_message
    character(len=llen_) :: description
    ! USER-DEFINED INPUTS
    integer              :: MGL                ! Number of multigrid levels
    integer, allocatable :: iter_threshold(:)  ! Number of iterations for each level
  end type hydra_multigrid_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(hydra_simulation_parameters_t), public :: hydra_obj_sim_param
  type(hydra_io_t), public                    :: hydra_obj_io
  type(hydra_multigrid_t), public             :: hydra_obj_multigrid
  !! ------------------------------------------------------
  !! ------------------------------------------------------
  
end module hydra_Config_Types_m