module hydraMF_Mod_BC_Coupled_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: BC_Coupled_Fluxes
  
contains

  subroutine BC_Coupled_Fluxes ( FUSS_domain, MOSE_domain, level )
    use hydraMF_Lib_Ghost
    use hydraMF_Lib_BC_Coupled_Fluxes
    use hydraMF_Mod_CouplingExchange, only: exchange_coupling_Tg_P, exchange_coupling_Pg_matID
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level

#ifdef USE_MPI
    ! Exchange cross-solver coupling data for entries whose secondary block is
    ! on a different MPI rank.  Populates bc(i)%Tg / %Pg and the remote P /
    ! matID caches before the local fill loops run.
    call exchange_coupling_Tg_P    ( FUSS_domain, MOSE_domain, level )
    call exchange_coupling_Pg_matID( FUSS_domain, MOSE_domain, level )
#endif

    ! Fill Ghost Cells (local-local entries only; remote entries filled above)
    call Fill_Coupled_Ghost_Cell ( FUSS_domain, MOSE_domain, level )

    ! Initialize ext_flux variables
    call Initialize_BC_Coupled_Fluxes ( FUSS_domain, MOSE_domain, level )

    ! Compute FUSS Flux
    call Compute_FUSS_Flux ( FUSS_domain, MOSE_domain, level )

    ! Compute MOSE Flux
    call Compute_MOSE_Flux ( FUSS_domain, MOSE_domain, level )

  end subroutine BC_Coupled_Fluxes


  !> Clear the accumulator on the coupling entries this rank owns.
  !>
  !> Only owned entries need clearing: every write into ext_flux -- the local
  !> update in Compute_*_Flux and the remote one applied by the reverse-flux
  !> exchange -- targets an entry whose primary block is local, and the solvers'
  !> BC loops read ext_flux only for blocks they own (Mod_BC_Fluxes cycles on
  !> is_local_block). Clearing the whole array meant sweeping it on every rank.
  subroutine Initialize_BC_Coupled_Fluxes ( FUSS_domain, MOSE_domain, level )
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    use hydraMF_Mod_CouplingExchange, only: own103_FUSS, own103_MOSE
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level
    ! Local
    integer :: k, n_F, n_M
    integer, pointer :: ownF(:), ownM(:)

    call own103_FUSS(level, n_F, ownF)
    call own103_MOSE(level, n_M, ownM)

    if (n_F > 0) then
      !$omp parallel do schedule (static) private(k)
      do k = 1, n_F
        FUSS_domain % bc(ownF(k)) % ext_flux = 0.0d0
      enddo
      !$omp end parallel do
    end if

    if (n_M > 0) then
      !$omp parallel do schedule (static) private(k)
      do k = 1, n_M
        MOSE_domain % bc(ownM(k)) % ext_flux = 0.0d0
      enddo
      !$omp end parallel do
    end if

  end subroutine Initialize_BC_Coupled_Fluxes

end module hydraMF_Mod_BC_Coupled_Fluxes