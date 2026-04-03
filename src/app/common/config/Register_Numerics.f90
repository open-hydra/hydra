module hydra_Read_Numerics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: Register_Numerics


contains

  subroutine Register_Numerics (nmgl)
    implicit none
    ! Local
    integer, intent(in)  :: nmgl

    !! ------------------------------------------------------
    !! Space Scheme ------------------------------------------
    !! ------------------------------------------------------

    ! Multigrid levels --------------------------------------
    call Register_Multigrid_Levels(nmgl)

    !! ------------------------------------------------------
    !! ------------------------------------------------------

  end subroutine Register_Numerics


  subroutine Register_Multigrid_Levels(nmgl)
    use hydra_Input_Registry
    use hydra_Config_Types_m, only: hydra_obj_multigrid
    use hydra_Parameters_m
    use IR_precision
    implicit none
    integer, intent(in)  :: nmgl
    character(len=llen_) :: option
    integer :: m

    hydra_obj_multigrid%warning_message = 'none'
    hydra_obj_multigrid%error_message   = 'none'
    hydra_obj_multigrid%description     = 'none'

    ! Allocate iteration threshold array for each level
    hydra_obj_multigrid%MGL = max(nmgl, 1)
    allocate( hydra_obj_multigrid%iter_threshold(hydra_obj_multigrid%MGL) )
    hydra_obj_multigrid%iter_threshold = 0

    ! Register per-level iteration parameters (level1-iter, level2-iter, level3-iter, ...)
    do m = 1, hydra_obj_multigrid%MGL
      write(option,'(A5,I0,A5)') 'level', m, '-iter'
      call reg%add( 'HYDRA-multigrid', trim(option), hydra_obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .true. )
    enddo

  end subroutine Register_Multigrid_Levels  


end module hydra_Read_Numerics