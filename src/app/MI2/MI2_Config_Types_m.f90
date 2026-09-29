module hydraMI2_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Parameters_m

  implicit none
  private

  !! ------------------------------------------------------
  !! MI2 coupling parameters ------------------------------
  !! ------------------------------------------------------
  type :: hydra_MI2_parameters_t
    ! USER-DEFINED INPUTS
    integer,  allocatable, dimension(:) :: ICE_iter    ! ICE   iterations per coupling step, per cycle
    integer,  allocatable, dimension(:) :: IGLOO_iter  ! IGLOO mode per cycle: 0=off, 1=source, 2=Eulerian
    logical                             :: relaxation  ! under-relax the IGLOO->MOSE source
    real(R8), allocatable, dimension(:) :: relCoeff    ! relaxation coefficient, one per cycle
    ! Raw relaxation-coeff read buffer: one value, or one per cycle
    real(R8), allocatable, dimension(:) :: relCoeff_raw
    ! DERIVED (set by Normalize_MI2_Params, not read from file)
    logical :: run_ICE   = .false.   ! true when ICE-iter   is present in input.ini
    logical :: run_IGLOO = .false.   ! true when IGLOO-iter is present in input.ini
    integer :: n_relax_given = 0     ! how many relaxation-coeff values the file actually carried
  end type hydra_MI2_parameters_t
  !! ------------------------------------------------------

  type(hydra_MI2_parameters_t), public :: hydra_MI2_obj

end module hydraMI2_Config_Types_m
