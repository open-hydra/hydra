module hydraMF_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m
  use hydra_Parameters_m
  use MOSE_Advanced_Types_m
  use FUSS_Advanced_Types_m
  use hydraMF_Mod_CouplingExchange, only: is_FUSS_bc_local_local, is_MOSE_bc_local_local, &
                                          own103_FUSS, own103_MOSE, &
                                          set_partner_half_FUSS, set_partner_half_MOSE
  use hydraMF_Lib_Metrics,          only: FUSS_Wall_Half, MOSE_Wall_Half

  implicit none
  private
  public :: Fill_Coupled_Ghost_Cell

contains


  !> Fill the coupling ghost cells for entries this rank can serve locally: the
  !> temperatures of the partner cells across the face (Tg / Pg, 1..gc_) and the
  !> half (A, Teff) of the interface balance of the partner face, computed in the
  !> partner's own index frame with its own metric (hydraMF_Lib_Metrics).
  !>
  !> Driven by the per-rank type(103) lists rather than by a sweep of the whole
  !> bc array: the array is global and grows with the block count, so the old
  !> form did more work on every rank as ranks were added. The local-local test
  !> stays, because an entry whose primary block is local may still have its
  !> secondary block on another rank -- those are served by the exchange that
  !> runs before this routine.
  subroutine Fill_Coupled_Ghost_Cell ( FUSS_domain, MOSE_domain, level )
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level
    ! Local
    integer :: k, n_F, n_M, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    real(R8) :: A, Teff
    integer, pointer :: ownF(:), ownM(:)

    call own103_FUSS(level, n_F, ownF)
    call own103_MOSE(level, n_M, ownM)

    ! FUSS: only fill locally when both Bm and Bs are owned by this rank.
    ! Remote-Bs entries are handled by exchange_coupling_Tg_P (called beforehand),
    ! and remote-Bm entries are not in ownF at all.
    if (n_F > 0) then
      !$omp parallel do schedule (static) &
      !$omp private(k, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
      do k = 1, n_F
        i = ownF(k)
        if (.not. is_FUSS_bc_local_local(FUSS_domain, i)) cycle
        Bm = FUSS_domain % bc(i) % b
        Im = FUSS_domain % bc(i) % i
        Jm = FUSS_domain % bc(i) % j
        Km = FUSS_domain % bc(i) % k
        Fm = FUSS_domain % bc(i) % f
        Bs = FUSS_domain % bc(i) % bs
        Is = FUSS_domain % bc(i) % is
        Js = FUSS_domain % bc(i) % js
        Ks = FUSS_domain % bc(i) % ks
        Fs = FUSS_domain % bc(i) % fs
        call Ghost_Coupling_FUSS ( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                   Is, Js, Ks, Fs, MOSE_domain % blk(Bs), &
                                   FUSS_domain % bc(i) % Tg )
      enddo
      !$omp end parallel do

      ! Half of the balance of the MOSE partner face
      !$omp parallel do schedule (static) &
      !$omp private(k, i, Bs, Is, Js, Ks, Fs, A, Teff)
      do k = 1, n_F
        i = ownF(k)
        if (.not. is_FUSS_bc_local_local(FUSS_domain, i)) cycle
        Bs = FUSS_domain % bc(i) % bs
        Is = FUSS_domain % bc(i) % is
        Js = FUSS_domain % bc(i) % js
        Ks = FUSS_domain % bc(i) % ks
        Fs = FUSS_domain % bc(i) % fs
        call MOSE_Wall_Half ( MOSE_domain % blk(Bs), Is, Js, Ks, Fs, A, Teff )
        call set_partner_half_FUSS ( level, i, A, Teff )
      enddo
      !$omp end parallel do
    end if

    ! MOSE: same guard -- skip entries whose secondary block is remote.
    if (n_M > 0) then
      !$omp parallel do schedule (static) &
      !$omp private(k, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
      do k = 1, n_M
        i = ownM(k)
        if (.not. is_MOSE_bc_local_local(MOSE_domain, i)) cycle
        Bm = MOSE_domain % bc(i) % b
        Im = MOSE_domain % bc(i) % i
        Jm = MOSE_domain % bc(i) % j
        Km = MOSE_domain % bc(i) % k
        Fm = MOSE_domain % bc(i) % f
        Bs = MOSE_domain % bc(i) % bs
        Is = MOSE_domain % bc(i) % is
        Js = MOSE_domain % bc(i) % js
        Ks = MOSE_domain % bc(i) % ks
        Fs = MOSE_domain % bc(i) % fs
        call Ghost_Coupling_MOSE ( Im, Jm, Km, Fm, MOSE_domain % blk(Bm), &
                                   Is, Js, Ks, Fs, FUSS_domain % blk(Bs), &
                                   MOSE_domain % bc(i) % Pg )
      enddo
      !$omp end parallel do

      ! Half of the balance of the FUSS partner face
      !$omp parallel do schedule (static) &
      !$omp private(k, i, Bs, Is, Js, Ks, Fs, A, Teff)
      do k = 1, n_M
        i = ownM(k)
        if (.not. is_MOSE_bc_local_local(MOSE_domain, i)) cycle
        Bs = MOSE_domain % bc(i) % bs
        Is = MOSE_domain % bc(i) % is
        Js = MOSE_domain % bc(i) % js
        Ks = MOSE_domain % bc(i) % ks
        Fs = MOSE_domain % bc(i) % fs
        call FUSS_Wall_Half ( FUSS_domain % blk(Bs), Is, Js, Ks, Fs, A, Teff )
        call set_partner_half_MOSE ( level, i, A, Teff )
      enddo
      !$omp end parallel do
    end if

  end subroutine Fill_Coupled_Ghost_Cell


  subroutine Ghost_Coupling_FUSS ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, Tg )
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, EOS
    use MOSE_Global_m, only: nsc, np
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(FUSS_block_type), intent(in) :: blkm
    type(MOSE_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Tg(6)
    ! Local
    integer  :: g, Is_, Js_, Ks_
    real(R8) :: rho, Rgas, Ts

    do g = 1, gc_

      Is_ = Is + guide_(Fs,1)*(g-1)
      Js_ = Js + guide_(Fs,2)*(g-1)
      Ks_ = Ks + guide_(Fs,3)*(g-1)
      
      call co_rotot_Rtot ( blks % P(1:nsc,Is_,Js_,Ks_), rho, Rgas )
      Ts = EOS ( p=blks%P(np,Is_,Js_,Ks_), rho=rho, R=Rgas )
      Tg(g) = Ts

    enddo

  end subroutine Ghost_Coupling_FUSS
  

  subroutine Ghost_Coupling_MOSE ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, Pg )
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, EOS
    use MOSE_Global_m, only: nsc, np
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(MOSE_block_type), intent(in) :: blkm
    type(FUSS_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Pg(1,6)
    ! Local
    integer  :: g, Is_, Js_, Ks_

    do g = 1, gc_

      Is_ = Is + guide_(Fs,1)*(g-1)
      Js_ = Js + guide_(Fs,2)*(g-1)
      Ks_ = Ks + guide_(Fs,3)*(g-1)
      
      Pg(1,g) = blks % T (Is_,Js_,Ks_)

    enddo

  end subroutine Ghost_Coupling_MOSE


end module hydraMF_Lib_Ghost