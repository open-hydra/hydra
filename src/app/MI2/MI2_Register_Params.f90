module hydraMI2_Read_Params
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Register_MI2_Params

contains

  !> Register the MI2-only keys under [HYDRA-Parameters]; all optional, an absent X-iter
  !  disables that solver. The arrays must already be sized by Scan_Ini_MI2.
  subroutine Register_MI2_Params()
    use hydra_Input_Registry
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    character(len=:), allocatable :: section

    section = 'HYDRA-Parameters'

    call reg%add(section,'ICE-iter', hydra_MI2_obj%ICE_iter, '0', &
                 'Number of ICE iterations per coupling iteration (absent => ICE disabled)', &
                 'integer array', .false. )

    call reg%add(section,'IGLOO-iter', hydra_MI2_obj%IGLOO_iter, '0', &
                 'IGLOO mode per cycle: 0=off, 1=source coupling, 2=Eulerian coupling '// &
                 '(absent => IGLOO disabled)', 'integer array', .false. )

    call reg%add(section,'relaxation', hydra_MI2_obj%relaxation, '.false.', &
                 'Under-relax the IGLOO->MOSE source term', 'logical', .false. )

    call reg%add(section,'relaxation-coeff', hydra_MI2_obj%relCoeff_raw, '0.5', &
                 'Under-relaxation coefficient: one value, or one per coupling cycle', &
                 'real array', .false. )

  end subroutine Register_MI2_Params

end module hydraMI2_Read_Params
