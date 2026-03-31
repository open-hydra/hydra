module hydraMF_Mod_BC_Coupled_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: BC_Coupled_Fluxes
  
contains

  subroutine BC_Coupled_Fluxes ( FUSS_domain, MOSE_domain )
    use hydraMF_Lib_Ghost
    use hydraMF_Lib_BC_Coupled_Fluxes
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain

    ! Fill Ghost Cells
    call Fill_Coupled_Ghost_Cell ( FUSS_domain, MOSE_domain )

    ! Initialize ext_flux variables
    call Initialize_BC_Coupled_Fluxes ( FUSS_domain, MOSE_domain )

    ! Compute FUSS Flux
    call Compute_FUSS_Flux ( FUSS_domain, MOSE_domain)

    ! Compute MOSE Flux
    call Compute_MOSE_Flux ( FUSS_domain, MOSE_domain)

  end subroutine BC_Coupled_Fluxes


  subroutine Initialize_BC_Coupled_Fluxes ( FUSS_domain, MOSE_domain )
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    ! Local
    integer  :: b, f, lower, upper, i
    
    !$omp parallel
    upper = 0
    do b = 1, FUSS_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + FUSS_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i)
        do i = lower, upper
          select case ( FUSS_domain%bc(i)%type )
            case (1000)
              FUSS_domain % bc(i) % ext_flux = 0.0d0
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

    !$omp parallel
    upper = 0
    do b = 1, MOSE_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + MOSE_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i)
        do i = lower, upper
          select case ( MOSE_domain%bc(i)%type )
            case (1000)
              MOSE_domain % bc(i) % ext_flux = 0.0d0
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

  end subroutine Initialize_BC_Coupled_Fluxes

end module hydraMF_Mod_BC_Coupled_Fluxes