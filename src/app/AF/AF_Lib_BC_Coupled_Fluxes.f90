module hydraAF_Lib_BC_Coupled_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_FUSS_Flux, Compute_ARES_Flux
  
contains

  subroutine Compute_FUSS_Flux ( FUSS_domain, ARES_domain, level )
    use hydraAF_Lib_BC_Coupled_Wall_Temperature, only: FUSS_BC_Coupled_Wall_Temperature
    use hydraAF_Mod_CouplingExchange, only: partner_of_FUSS_bc, &
                                            exchange_coupling_reverse_flux_FUSS
    use FUSS_Advanced_Types_m
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: np, nh
#ifdef USE_MPI
    use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner, mpi_rank_
    use ARES_Mod_MPI, only: ARES_block_owner => block_owner
#endif
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    integer, intent(in) :: level
    ! Local
    integer  :: b, f, lower, upper, i, j, ares_bm
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    real(R8) :: FUSS_new_flux
    real(R8) :: P_con(np)
    real(R8), allocatable :: new_flux_arr(:)  !< per-entry flux for remote ARES partner update

    allocate(new_flux_arr(FUSS_domain%nbound)); new_flux_arr = 0.0d0

    !$omp parallel private(b, f, lower, upper)
    upper = 0
    do b = 1, FUSS_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + FUSS_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i,j,ares_bm,Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs) &
        !$omp private (FUSS_new_flux, P_con)
        do i = lower, upper
          select case ( FUSS_domain%bc(i)%type )
            case (103)
              Bm  = FUSS_domain % bc(i) % b
              ! Only process entries whose primary block (Bm) is owned by this rank.
#ifdef USE_MPI
              if (allocated(FUSS_block_owner)) then
                if (FUSS_block_owner(Bm) /= mpi_rank_) cycle
              end if
#endif
              Im  = FUSS_domain % bc(i) % i
              Jm  = FUSS_domain % bc(i) % j
              Km  = FUSS_domain % bc(i) % k
              Fm  = FUSS_domain % bc(i) % f
              Bs  = FUSS_domain % bc(i) % bs
              Is  = FUSS_domain % bc(i) % is
              Js  = FUSS_domain % bc(i) % js
              Ks  = FUSS_domain % bc(i) % ks
              Fs  = FUSS_domain % bc(i) % fs

              ! Fetch P(1:np) of the ARES secondary cell: direct access when Bs is
              ! local, or from the remote cache populated by exchange_coupling_Tg_P.
              call get_secondary_P(ARES_domain, i, Bs, Is, Js, Ks, level, P_con)

              call FUSS_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                          Is, Js, Ks, Fs, P_con, &
                                          FUSS_domain%bc(i)%Mg(1), FUSS_domain%bc(i)%Tg, FUSS_new_flux )
              FUSS_domain % bc(i) % ext_flux = FUSS_domain % bc(i) % ext_flux + 0.5d0 * FUSS_new_flux
              new_flux_arr(i) = FUSS_new_flux

              ! Apply reverse flux into the partner ARES bc when partner is local.
              j = partner_of_FUSS_bc(level, i)
              if (j > 0) then
#ifdef USE_MPI
                ares_bm = ARES_domain%bc(j)%b
                if (.not. allocated(ARES_block_owner) .or. &
                    ARES_block_owner(ares_bm) == mpi_rank_) then
                  ARES_domain%bc(j)%ext_flux(nh) = ARES_domain%bc(j)%ext_flux(nh) &
                                                 - 0.5d0 * FUSS_new_flux
                end if
#else
                ARES_domain%bc(j)%ext_flux(nh) = ARES_domain%bc(j)%ext_flux(nh) &
                                               - 0.5d0 * FUSS_new_flux
#endif
              end if
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

    ! For FUSS-primary BCs whose ARES partner lives on a different rank, send
    ! the accumulated flux there to apply as -0.5*flux into ARES ext_flux(nh).
#ifdef USE_MPI
    call exchange_coupling_reverse_flux_FUSS(new_flux_arr, FUSS_domain, ARES_domain, level)
#endif

    deallocate(new_flux_arr)

  end subroutine Compute_FUSS_Flux


  subroutine Compute_ARES_Flux ( FUSS_domain, ARES_domain, level )
    use hydraAF_Lib_BC_Coupled_Wall_Temperature, only: ARES_BC_Coupled_Wall_Temperature
    use hydraAF_Mod_CouplingExchange, only: partner_of_ARES_bc, &
                                            exchange_coupling_reverse_flux_ARES
    use FUSS_Advanced_Types_m
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: nprim, nh
#ifdef USE_MPI
    use ARES_Mod_MPI, only: ARES_block_owner => block_owner, mpi_rank_
    use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner
#endif
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    integer, intent(in) :: level
    ! Local
    integer  :: b, f, lower, upper, i, j, fuss_bm
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    real(R8) :: ARES_new_flux(nprim)
    real(R8) :: matID_con
    real(R8), allocatable :: new_flux_nh_arr(:)

    allocate(new_flux_nh_arr(ARES_domain%nbound)); new_flux_nh_arr = 0.0d0

    !$omp parallel private(b, f, lower, upper)
    upper = 0
    do b = 1, ARES_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + ARES_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i,j,fuss_bm,Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs) &
        !$omp private (ARES_new_flux, matID_con)
        do i = lower, upper
          select case ( ARES_domain%bc(i)%type )
            case (103)
              Bm  = ARES_domain % bc(i) % b
#ifdef USE_MPI
              if (allocated(ARES_block_owner)) then
                if (ARES_block_owner(Bm) /= mpi_rank_) cycle
              end if
#endif
              Im  = ARES_domain % bc(i) % i
              Jm  = ARES_domain % bc(i) % j
              Km  = ARES_domain % bc(i) % k
              Fm  = ARES_domain % bc(i) % f
              Bs  = ARES_domain % bc(i) % bs
              Is  = ARES_domain % bc(i) % is
              Js  = ARES_domain % bc(i) % js
              Ks  = ARES_domain % bc(i) % ks
              Fs  = ARES_domain % bc(i) % fs

              call get_secondary_matID(FUSS_domain, ARES_domain, i, Bs, Is, Js, Ks, &
                                       level, matID_con)

              call ARES_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, ARES_domain % blk(Bm), &
                                          Is, Js, Ks, Fs, matID_con, &
                                          ARES_domain%bc(i)%Mg(1), ARES_domain%bc(i)%Pg, ARES_new_flux )
              ARES_domain % bc(i) % ext_flux = ARES_domain % bc(i) % ext_flux + 0.5d0 * ARES_new_flux
              new_flux_nh_arr(i) = ARES_new_flux(nh)

              j = partner_of_ARES_bc(level, i)
              if (j > 0) then
#ifdef USE_MPI
                fuss_bm = FUSS_domain%bc(j)%b
                if (.not. allocated(FUSS_block_owner) .or. &
                    FUSS_block_owner(fuss_bm) == mpi_rank_) then
                  FUSS_domain%bc(j)%ext_flux = FUSS_domain%bc(j)%ext_flux &
                                             - 0.5d0 * ARES_new_flux(nh)
                end if
#else
                FUSS_domain%bc(j)%ext_flux = FUSS_domain%bc(j)%ext_flux &
                                           - 0.5d0 * ARES_new_flux(nh)
#endif
              end if
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

#ifdef USE_MPI
    call exchange_coupling_reverse_flux_ARES(new_flux_nh_arr, FUSS_domain, ARES_domain, level)
#endif

    deallocate(new_flux_nh_arr)

  end subroutine Compute_ARES_Flux


  !---------------------------------------------------------------------------
  !> Return P(1:np, Is,Js,Ks) of ARES block Bs: direct access when Bs is local,
  !> otherwise use the remote cache populated by exchange_coupling_Tg_P.
  subroutine get_secondary_P(ARES_domain, i, Bs, Is, Js, Ks, level, P_con)
    use ARES_Global_m, only: np
    use ARES_Advanced_Types_m
    use hydraAF_Mod_CouplingExchange, only: get_F_remote_P
#ifdef USE_MPI
    use ARES_Mod_MPI, only: ARES_block_owner => block_owner, mpi_rank_
#endif
    implicit none
    type(ARES_domain_type), intent(in) :: ARES_domain
    integer, intent(in) :: i, Bs, Is, Js, Ks, level
    real(R8), intent(out) :: P_con(:)
    real(R8), pointer :: p_remote(:)
#ifdef USE_MPI
    if (allocated(ARES_block_owner)) then
      if (ARES_block_owner(Bs) /= mpi_rank_) then
        p_remote => get_F_remote_P(level, i)
        if (associated(p_remote)) then
          P_con(1:np) = p_remote(1:np)
          return
        end if
      end if
    end if
#endif
    P_con(1:np) = ARES_domain%blk(Bs)%P(1:np, Is, Js, Ks)
  end subroutine get_secondary_P


  !---------------------------------------------------------------------------
  !> Return matID(Is,Js,Ks) of FUSS block Bs: direct access when Bs is local,
  !> otherwise use the remote cache populated by exchange_coupling_Pg_matID.
  subroutine get_secondary_matID(FUSS_domain, ARES_domain, i, Bs, Is, Js, Ks, level, matID_con)
    use FUSS_Advanced_Types_m
    use ARES_Advanced_Types_m
    use hydraAF_Mod_CouplingExchange, only: get_M_remote_matID
#ifdef USE_MPI
    use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner, mpi_rank_
#endif
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(ARES_domain_type), intent(in) :: ARES_domain
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
    if (ARES_domain%nbound < 0) return  ! silence unused-arg warning
  end subroutine get_secondary_matID

end module hydraAF_Lib_BC_Coupled_Fluxes