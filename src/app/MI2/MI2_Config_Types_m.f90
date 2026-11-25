module hydraMI2_Config_Types_m
  use iso_fortran_env, only: R8 => real64

  implicit none
  private

  !! ------------------------------------------------------
  !! MI2 coupling parameters ------------------------------
  !! ------------------------------------------------------
  type :: hydra_MI2_parameters_t
    ! [HYDRA-Parameters]
    integer,  allocatable, dimension(:) :: ICE_iter    ! ICE   iterations per coupling step, per cycle
    integer,  allocatable, dimension(:) :: IGLOO_iter  ! IGLOO mode per cycle: 0=off, 1=source, 2=Eulerian
    logical                             :: relaxation  ! under-relax the IGLOO->MOSE source
    real(R8), allocatable, dimension(:) :: relCoeff    ! relaxation coefficient, one per cycle (1 when off)
    ! Raw relaxation-coeff read buffer: one value, or one per cycle
    real(R8), allocatable, dimension(:) :: relCoeff_raw
    ! [MI2-Coupling]
    logical :: report           = .false.   ! print, per coupling step, what the coupler hands the gas
    logical :: datum_correction = .true.    ! apply the enthalpy-datum offsets (false: diagnostic only)
    ! DERIVED (set by Scan_Ini_MI2 from what input.ini carries)
    logical :: run_ICE   = .false.   ! true when ICE-iter   is present in input.ini
    logical :: run_IGLOO = .false.   ! true when IGLOO-iter is present in input.ini
    integer :: n_relax_given = 0     ! how many relaxation-coeff values the file actually carried
    ! DERIVED after MOSE's set-up (Check_MOSE_Grid_Levels)
    logical :: coarse_start = .false.   ! MOSE starts on its coarse grids: a [MOSE-Multigrid] level<m>-iter > 0, m > 1
  end type hydra_MI2_parameters_t
  !! ------------------------------------------------------

  type(hydra_MI2_parameters_t), public :: hydra_MI2_obj

end module hydraMI2_Config_Types_m
