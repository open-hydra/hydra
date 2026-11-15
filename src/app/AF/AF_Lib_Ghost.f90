module hydraAF_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m
  use hydra_Parameters_m
  use ARES_Advanced_Types_m
  use FUSS_Advanced_Types_m
  use hydraAF_Mod_CouplingExchange, only: is_FUSS_bc_local_local, is_ARES_bc_local_local, &
                                          own103_FUSS, own103_ARES, &
                                          set_partner_half_FUSS, set_partner_half_ARES
  use hydraAF_Lib_Metrics,          only: FUSS_Wall_Half, ARES_Wall_Half

  implicit none
  private
  public :: Fill_Coupled_Ghost_Cell

contains


  !> Fill the coupling ghost cells for entries this rank can serve locally: the
  !> temperatures of the partner cells across the face (Tg / Pg, 1..gc_) and the
  !> half (A, Teff) of the interface balance of the partner face, computed in the
  !> partner's own index frame with its own metric (hydraAF_Lib_Metrics).
  !>
  !> Driven by the per-rank type(103) lists rather than by a sweep of the whole
  !> bc array: the array is global and grows with the block count, so the old
  !> form did more work on every rank as ranks were added. The local-local test
  !> stays, because an entry whose primary block is local may still have its
  !> secondary block on another rank -- those are served by the exchange that
  !> runs before this routine.
  subroutine Fill_Coupled_Ghost_Cell ( FUSS_domain, ARES_domain, level )
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    integer, intent(in) :: level
    ! Local
    integer :: k, n_F, n_A, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    real(R8) :: A, Teff
    integer, pointer :: ownF(:), ownA(:)

    call own103_FUSS(level, n_F, ownF)
    call own103_ARES(level, n_A, ownA)

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
                                   Is, Js, Ks, Fs, ARES_domain % blk(Bs), &
                                   FUSS_domain % bc(i) % Tg )
      enddo
      !$omp end parallel do

      ! Half of the balance of the ARES partner face
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
        call ARES_Wall_Half ( ARES_domain % blk(Bs), Is, Js, Ks, Fs, A, Teff )
        call set_partner_half_FUSS ( level, i, A, Teff )
      enddo
      !$omp end parallel do
    end if

    ! ARES: same guard -- skip entries whose secondary block is remote.
    if (n_A > 0) then
      !$omp parallel do schedule (static) &
      !$omp private(k, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
      do k = 1, n_A
        i = ownA(k)
        if (.not. is_ARES_bc_local_local(ARES_domain, i)) cycle
        Bm = ARES_domain % bc(i) % b
        Im = ARES_domain % bc(i) % i
        Jm = ARES_domain % bc(i) % j
        Km = ARES_domain % bc(i) % k
        Fm = ARES_domain % bc(i) % f
        Bs = ARES_domain % bc(i) % bs
        Is = ARES_domain % bc(i) % is
        Js = ARES_domain % bc(i) % js
        Ks = ARES_domain % bc(i) % ks
        Fs = ARES_domain % bc(i) % fs
        call Ghost_Coupling_ARES ( Im, Jm, Km, Fm, ARES_domain % blk(Bm), &
                                   Is, Js, Ks, Fs, FUSS_domain % blk(Bs), &
                                   ARES_domain % bc(i) % Pg )
      enddo
      !$omp end parallel do

      ! Half of the balance of the FUSS partner face
      !$omp parallel do schedule (static) &
      !$omp private(k, i, Bs, Is, Js, Ks, Fs, A, Teff)
      do k = 1, n_A
        i = ownA(k)
        if (.not. is_ARES_bc_local_local(ARES_domain, i)) cycle
        Bs = ARES_domain % bc(i) % bs
        Is = ARES_domain % bc(i) % is
        Js = ARES_domain % bc(i) % js
        Ks = ARES_domain % bc(i) % ks
        Fs = ARES_domain % bc(i) % fs
        call FUSS_Wall_Half ( FUSS_domain % blk(Bs), Is, Js, Ks, Fs, A, Teff )
        call set_partner_half_ARES ( level, i, A, Teff )
      enddo
      !$omp end parallel do
    end if

  end subroutine Fill_Coupled_Ghost_Cell


  subroutine Ghost_Coupling_FUSS ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, Tg )
    use FLINT_Lib_Thermodynamic
    use ARES_Global_m, only: np, nh
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(FUSS_block_type), intent(in) :: blkm
    type(ARES_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Tg(6)
    ! Local
    integer  :: g, Is_, Js_, Ks_
    real(R8) :: Ts

    do g = 1, gc_

      Is_ = Is + guide_(Fs,1)*(g-1)
      Js_ = Js + guide_(Fs,2)*(g-1)
      Ks_ = Ks + guide_(Fs,3)*(g-1)
      
      Ts = ph2vars(blks%P(np,Is_,Js_,Ks_), blks%P(nh,Is_,Js_,Ks_), T_tab )
      Tg(g) = Ts

    enddo

  end subroutine Ghost_Coupling_FUSS
  

  subroutine Ghost_Coupling_ARES ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, Pg )
    use FLINT_Lib_Thermodynamic
    use ARES_Global_m, only: np, nh
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(ARES_block_type), intent(in) :: blkm
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

  end subroutine Ghost_Coupling_ARES


end module hydraAF_Lib_Ghost