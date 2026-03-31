module hydraMF_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Fix_Global_dt
  
contains

  subroutine Fix_Global_dt ( FUSS_simulation, MOSE_simulation, hydra_time )
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    use FUSS_Mod_dt, only: FUSS_Compute_dt => Compute_dt
    use MOSE_Mod_dt, only: MOSE_Compute_dt => Compute_dt
    use FUSS_Config_Types_m, only: FUSS_obj_time_scheme => obj_time_scheme
    use MOSE_Config_Types_m, only: MOSE_obj_time_scheme => obj_time_scheme
    implicit none
    type(FUSS_simulation_type), intent(inout) :: FUSS_simulation
    type(MOSE_simulation_type), intent(inout) :: MOSE_simulation
    real(R8), intent(inout)                   :: hydra_time
    ! Local
    real(R8) :: min_dt

    ! Compute Global dt
    FUSS_simulation % domain(1) % dtglobal = 1e5
    MOSE_simulation % domain(1) % dtglobal = 1e5
    call FUSS_Compute_dt( FUSS_simulation % domain(1), FUSS_obj_time_scheme%vnn, FUSS_obj_time_scheme%rampa_vnn_iter )
    call MOSE_Compute_dt( MOSE_simulation % domain(1), MOSE_obj_time_scheme%cfl, MOSE_obj_time_scheme%vnn, MOSE_obj_time_scheme%rampa_cfl_iter )

    ! Assign Global dt
    min_dt = min(FUSS_simulation % domain(1) % dtglobal, MOSE_simulation % domain(1) % dtglobal)
    FUSS_simulation % domain(1) % dtglobal = min_dt
    MOSE_simulation % domain(1) % dtglobal = min_dt

    hydra_time = hydra_time + min_dt

  end subroutine Fix_Global_dt

end module hydraMF_Mod_dt