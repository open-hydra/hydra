module hydraAF_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Fix_Global_dt
  
contains

  subroutine Fix_Global_dt ( FUSS_simulation, ARES_simulation, hydra_time )
    use FUSS_Advanced_Types_m
    use ARES_Advanced_Types_m
    use FUSS_Mod_dt, only: FUSS_Compute_dt => Compute_dt
    use ARES_Mod_dt, only: ARES_Compute_dt => Compute_dt
    use FUSS_Config_Types_m, only: FUSS_obj_time_scheme => obj_time_scheme
    use ARES_Config_Types_m, only: ARES_obj_time_scheme => obj_time_scheme
    implicit none
    type(FUSS_simulation_type), intent(inout) :: FUSS_simulation
    type(ARES_simulation_type), intent(inout) :: ARES_simulation
    real(R8), intent(inout)                   :: hydra_time
    ! Local
    real(R8) :: min_dt

    ! Compute Global dt
    FUSS_simulation % domain(1) % dtglobal = 1e5
    ARES_simulation % domain(1) % dtglobal = 1e5
    call FUSS_Compute_dt( FUSS_simulation % domain(1), FUSS_obj_time_scheme%vnn, FUSS_obj_time_scheme%rampa_vnn_iter )
    call ARES_Compute_dt( ARES_simulation % domain(1), ARES_obj_time_scheme%cfl, ARES_obj_time_scheme%vnn, ARES_obj_time_scheme%rampa_cfl_iter )

    ! Assign Global dt
    min_dt = min(FUSS_simulation % domain(1) % dtglobal, ARES_simulation % domain(1) % dtglobal)
    FUSS_simulation % domain(1) % dtglobal = min_dt
    ARES_simulation % domain(1) % dtglobal = min_dt

    hydra_time = hydra_time + min_dt

  end subroutine Fix_Global_dt

end module hydraAF_Mod_dt