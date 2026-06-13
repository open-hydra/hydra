module hydraAF_Mod_BC_Coupled_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: BC_Coupled_Fluxes
  
contains

  subroutine BC_Coupled_Fluxes ( FUSS_domain, ARES_domain, level )
    use hydraAF_Lib_Ghost
    use hydraAF_Lib_BC_Coupled_Fluxes
    use hydraAF_Mod_CouplingExchange, only: exchange_coupling_Tg_P, exchange_coupling_Pg_matID
    use FUSS_Advanced_Types_m
    use ARES_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    integer, intent(in) :: level

#ifdef USE_MPI
    ! Exchange cross-solver coupling data for entries whose secondary block is
    ! on a different MPI rank.  Populates bc(i)%Tg / %Pg and the remote P /
    ! matID caches before the local fill loops run.
    call exchange_coupling_Tg_P    ( FUSS_domain, ARES_domain, level )
    call exchange_coupling_Pg_matID( FUSS_domain, ARES_domain, level )
#endif

    ! Fill Ghost Cells (local-local entries only; remote entries filled above)
    call Fill_Coupled_Ghost_Cell ( FUSS_domain, ARES_domain )

    ! Initialize ext_flux variables
    call Initialize_BC_Coupled_Fluxes ( FUSS_domain, ARES_domain )

    ! Compute FUSS Flux
    call Compute_FUSS_Flux ( FUSS_domain, ARES_domain, level )

    ! Compute ARES Flux
    call Compute_ARES_Flux ( FUSS_domain, ARES_domain, level )

  end subroutine BC_Coupled_Fluxes


  subroutine Initialize_BC_Coupled_Fluxes ( FUSS_domain, ARES_domain )
    use FUSS_Advanced_Types_m
    use ARES_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    ! Local
    integer  :: b, f, lower, upper, i
    
    !$omp parallel private(b, f, lower, upper)
    upper = 0
    do b = 1, FUSS_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + FUSS_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i)
        do i = lower, upper
          select case ( FUSS_domain%bc(i)%type )
            case (103)
              FUSS_domain % bc(i) % ext_flux = 0.0d0
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

    !$omp parallel private(b, f, lower, upper)
    upper = 0
    do b = 1, ARES_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + ARES_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i)
        do i = lower, upper
          select case ( ARES_domain%bc(i)%type )
            case (103)
              ARES_domain % bc(i) % ext_flux = 0.0d0
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

  end subroutine Initialize_BC_Coupled_Fluxes

end module hydraAF_Mod_BC_Coupled_Fluxes