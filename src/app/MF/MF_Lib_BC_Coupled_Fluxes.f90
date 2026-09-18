module hydraMF_Lib_BC_Coupled_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_FUSS_Flux, Compute_MOSE_Flux
  
contains

  subroutine Compute_FUSS_Flux ( FUSS_domain, MOSE_domain, level )
    use hydraMF_Lib_BC_Coupled_Wall_Temperature, only: FUSS_BC_Coupled_Wall_Temperature
    use hydraMF_Mod_CouplingExchange, only: partner_of_FUSS_bc, own103_FUSS, &
                                            exchange_coupling_reverse_flux_FUSS
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: np, nsc
#ifdef USE_MPI
    use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner, mpi_rank_
    use MOSE_Mod_MPI, only: MOSE_block_owner => block_owner
#endif
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level
    ! Local
    integer  :: k, n_own, i, j, mose_bm
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    real(R8) :: FUSS_new_flux
    real(R8) :: P_con(nsc)
    integer,  pointer :: own(:)                !< type(103 entries whose Bm this rank owns
    real(R8), pointer :: new_flux_arr(:)       !< per-entry flux for remote MOSE partner update

    ! Iterate only this rank's coupling entries. The list and the bc-indexed
    ! scratch are built once per level in build_coupling_schedule; sweeping the
    ! whole bc array here made the cost grow with the block count instead of
    ! falling with the rank count.
    call own103_FUSS(level, n_own, own, new_flux_arr)

    if (n_own > 0) then
      !$omp parallel do schedule (static) &
      !$omp private (k, i, j, mose_bm, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs) &
      !$omp private (FUSS_new_flux, P_con)
      do k = 1, n_own
        i   = own(k)
        Bm  = FUSS_domain % bc(i) % b
        Im  = FUSS_domain % bc(i) % i
        Jm  = FUSS_domain % bc(i) % j
        Km  = FUSS_domain % bc(i) % k
        Fm  = FUSS_domain % bc(i) % f
        Bs  = FUSS_domain % bc(i) % bs
        Is  = FUSS_domain % bc(i) % is
        Js  = FUSS_domain % bc(i) % js
        Ks  = FUSS_domain % bc(i) % ks
        Fs  = FUSS_domain % bc(i) % fs

        ! Fetch P(1:nsc) of the MOSE secondary cell: direct access when Bs is
        ! local, or from the remote cache populated by exchange_coupling_Tg_P.
        call get_secondary_P(MOSE_domain, i, Bs, Is, Js, Ks, level, P_con)

        call FUSS_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                    Is, Js, Ks, Fs, P_con, &
                                    FUSS_domain%bc(i)%Mg(1), FUSS_domain%bc(i)%Tg, FUSS_new_flux )
        FUSS_domain % bc(i) % ext_flux = FUSS_domain % bc(i) % ext_flux + 0.5d0 * FUSS_new_flux
        new_flux_arr(i) = FUSS_new_flux

        ! Apply reverse flux into the partner MOSE bc when partner is local.
        j = partner_of_FUSS_bc(level, i)
        if (j > 0) then
#ifdef USE_MPI
          mose_bm = MOSE_domain%bc(j)%b
          if (.not. allocated(MOSE_block_owner) .or. &
              MOSE_block_owner(mose_bm) == mpi_rank_) then
            MOSE_domain%bc(j)%ext_flux(np) = MOSE_domain%bc(j)%ext_flux(np) &
                                           - 0.5d0 * FUSS_new_flux
          end if
#else
          MOSE_domain%bc(j)%ext_flux(np) = MOSE_domain%bc(j)%ext_flux(np) &
                                         - 0.5d0 * FUSS_new_flux
#endif
        end if
      enddo
      !$omp end parallel do
    end if

    ! For FUSS-primary BCs whose MOSE partner lives on a different rank, send
    ! the accumulated flux there to apply as -0.5*flux into MOSE ext_flux(np).
    ! Not guarded by n_own: a rank owning no FUSS coupling entries can still own
    ! MOSE blocks that RECEIVE reverse flux, and skipping the exchange there
    ! would leave its partners' sends unmatched.
#ifdef USE_MPI
    if (associated(new_flux_arr)) &
      call exchange_coupling_reverse_flux_FUSS(new_flux_arr, FUSS_domain, MOSE_domain, level)
#endif

  end subroutine Compute_FUSS_Flux


  subroutine Compute_MOSE_Flux ( FUSS_domain, MOSE_domain, level )
    use hydraMF_Lib_BC_Coupled_Wall_Temperature, only: MOSE_BC_Coupled_Wall_Temperature
    use hydraMF_Mod_CouplingExchange, only: partner_of_MOSE_bc, own103_MOSE, &
                                            exchange_coupling_reverse_flux_MOSE
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: nprim, np
#ifdef USE_MPI
    use MOSE_Mod_MPI, only: MOSE_block_owner => block_owner, mpi_rank_
    use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner
#endif
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level
    ! Local
    integer  :: k, n_own, i, j, fuss_bm
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    real(R8) :: MOSE_new_flux(nprim)
    real(R8) :: matID_con
    integer,  pointer :: own(:)                !< type(103 entries whose Bm this rank owns
    real(R8), pointer :: new_flux_np_arr(:)

    call own103_MOSE(level, n_own, own, new_flux_np_arr)

    if (n_own > 0) then
      !$omp parallel do schedule (static) &
      !$omp private (k, i, j, fuss_bm, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs) &
      !$omp private (MOSE_new_flux, matID_con)
      do k = 1, n_own
        i   = own(k)
        Bm  = MOSE_domain % bc(i) % b
        Im  = MOSE_domain % bc(i) % i
        Jm  = MOSE_domain % bc(i) % j
        Km  = MOSE_domain % bc(i) % k
        Fm  = MOSE_domain % bc(i) % f
        Bs  = MOSE_domain % bc(i) % bs
        Is  = MOSE_domain % bc(i) % is
        Js  = MOSE_domain % bc(i) % js
        Ks  = MOSE_domain % bc(i) % ks
        Fs  = MOSE_domain % bc(i) % fs

        call get_secondary_matID(FUSS_domain, MOSE_domain, i, Bs, Is, Js, Ks, &
                                 level, matID_con)

        call MOSE_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, MOSE_domain % blk(Bm), &
                                    Is, Js, Ks, Fs, matID_con, &
                                    MOSE_domain%bc(i)%Mg(1), MOSE_domain%bc(i)%Pg, MOSE_new_flux, &
                                    k_rough = MOSE_domain%bc(i)%k_rough )
        MOSE_domain % bc(i) % ext_flux = MOSE_domain % bc(i) % ext_flux + 0.5d0 * MOSE_new_flux
        new_flux_np_arr(i) = MOSE_new_flux(np)

        j = partner_of_MOSE_bc(level, i)
        if (j > 0) then
#ifdef USE_MPI
          fuss_bm = FUSS_domain%bc(j)%b
          if (.not. allocated(FUSS_block_owner) .or. &
              FUSS_block_owner(fuss_bm) == mpi_rank_) then
            FUSS_domain%bc(j)%ext_flux = FUSS_domain%bc(j)%ext_flux &
                                       - 0.5d0 * MOSE_new_flux(np)
          end if
#else
          FUSS_domain%bc(j)%ext_flux = FUSS_domain%bc(j)%ext_flux &
                                     - 0.5d0 * MOSE_new_flux(np)
#endif
        end if
      enddo
      !$omp end parallel do
    end if

#ifdef USE_MPI
    if (associated(new_flux_np_arr)) &
      call exchange_coupling_reverse_flux_MOSE(new_flux_np_arr, FUSS_domain, MOSE_domain, level)
#endif

  end subroutine Compute_MOSE_Flux


  !---------------------------------------------------------------------------
  !> Return P(1:nsc, Is,Js,Ks) of MOSE block Bs: direct access when Bs is local,
  !> otherwise use the remote cache populated by exchange_coupling_Tg_P.
  subroutine get_secondary_P(MOSE_domain, i, Bs, Is, Js, Ks, level, P_con)
    use MOSE_Global_m, only: nsc
    use MOSE_Advanced_Types_m
    use hydraMF_Mod_CouplingExchange, only: get_F_remote_P
#ifdef USE_MPI
    use MOSE_Mod_MPI, only: MOSE_block_owner => block_owner, mpi_rank_
#endif
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in) :: i, Bs, Is, Js, Ks, level
    real(R8), intent(out) :: P_con(:)
    real(R8), pointer :: p_remote(:)
#ifdef USE_MPI
    if (allocated(MOSE_block_owner)) then
      if (MOSE_block_owner(Bs) /= mpi_rank_) then
        p_remote => get_F_remote_P(level, i)
        if (associated(p_remote)) then
          P_con(1:nsc) = p_remote(1:nsc)
          return
        end if
      end if
    end if
#endif
    P_con(1:nsc) = MOSE_domain%blk(Bs)%P(1:nsc, Is, Js, Ks)
  end subroutine get_secondary_P


  !---------------------------------------------------------------------------
  !> Return matID(Is,Js,Ks) of FUSS block Bs: direct access when Bs is local,
  !> otherwise use the remote cache populated by exchange_coupling_Pg_matID.
  subroutine get_secondary_matID(FUSS_domain, MOSE_domain, i, Bs, Is, Js, Ks, level, matID_con)
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    use hydraMF_Mod_CouplingExchange, only: get_M_remote_matID
#ifdef USE_MPI
    use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner, mpi_rank_
#endif
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in) :: i, Bs, Is, Js, Ks, level
    real(R8), intent(out) :: matID_con
#ifdef USE_MPI
    if (allocated(FUSS_block_owner)) then
      if (FUSS_block_owner(Bs) /= mpi_rank_) then
        matID_con = get_M_remote_matID(level, i)
        return
      end if
    end if
#endif
    matID_con = FUSS_domain%blk(Bs)%matID(Is, Js, Ks)
    if (MOSE_domain%nbound < 0) return  ! silence unused-arg warning
  end subroutine get_secondary_matID

end module hydraMF_Lib_BC_Coupled_Fluxes