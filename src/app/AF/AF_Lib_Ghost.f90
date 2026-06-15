module hydraAF_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m
  use hydra_Parameters_m
  use ARES_Advanced_Types_m
  use FUSS_Advanced_Types_m
  use hydraAF_Mod_CouplingExchange, only: is_FUSS_bc_local_local, is_ARES_bc_local_local

  implicit none
  private
  public :: Fill_Coupled_Ghost_Cell

contains


  subroutine Fill_Coupled_Ghost_Cell ( FUSS_domain, ARES_domain )
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    ! Local
    integer :: i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s

    ! FUSS: only fill locally when both Bm and Bs are owned by this rank.
    ! Remote-Bs entries are handled by exchange_coupling_Tg_P (called beforehand),
    ! and remote-Bm entries belong to another rank's FUSS partition.
    !$omp parallel
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
    do i = 1, FUSS_domain % nbound
      select case ( FUSS_domain % bc(i) % type)
        case(103)
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
      end select
    enddo
    !$omp end parallel

    !$omp parallel
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, FUSS_domain % nbound
      select case ( FUSS_domain % bc(i) % type)
        case(103)
          if (.not. is_FUSS_bc_local_local(FUSS_domain, i)) cycle
          Bs = FUSS_domain % bc(i) % bs
          Is = FUSS_domain % bc(i) % is 
          Js = FUSS_domain % bc(i) % js
          Ks = FUSS_domain % bc(i) % ks
          Fs = FUSS_domain % bc(i) % fs
          d11s = FUSS_domain % bc(i) % d11
          d12s = FUSS_domain % bc(i) % d12
          d21s = FUSS_domain % bc(i) % d21
          d22s = FUSS_domain % bc(i) % d22
          call Fill_BC_Ghost_Coupling_FUSS ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, &
                                             ARES_domain % blk(Bs), FUSS_domain % bc(i) % Tg )
      end select
    enddo
    !$omp end parallel

    ! ARES: same guard — skip entries where Bm or Bs is remote.
    !$omp parallel
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
    do i = 1, ARES_domain % nbound
      select case ( ARES_domain % bc(i) % type)
        case(103)
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
      end select
    enddo
    !$omp end parallel

    !$omp parallel
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, ARES_domain % nbound
      select case ( ARES_domain % bc(i) % type)
        case(103)
          if (.not. is_ARES_bc_local_local(ARES_domain, i)) cycle
          Bs = ARES_domain % bc(i) % bs
          Is = ARES_domain % bc(i) % is 
          Js = ARES_domain % bc(i) % js
          Ks = ARES_domain % bc(i) % ks
          Fs = ARES_domain % bc(i) % fs
          d11s = ARES_domain % bc(i) % d11
          d12s = ARES_domain % bc(i) % d12
          d21s = ARES_domain % bc(i) % d21
          d22s = ARES_domain % bc(i) % d22
          call Fill_BC_Ghost_Coupling_ARES ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, &
                                             FUSS_domain % blk(Bs), ARES_domain % bc(i) % Pg )
      end select
    enddo
    !$omp end parallel

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
  

  subroutine Fill_BC_Ghost_Coupling_FUSS ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, blks, Tg )
    use FLINT_Lib_Thermodynamic
    use ARES_Global_m, only: np, nh
    implicit none
    integer, intent(in)               :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(ARES_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Tg(6)
    ! Local
    integer  :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
    real(R8) :: Ts1, Ts2, Ts3, Ts4

    select case(Fs)
      case(1:2)
        i1 = Is
        j1 = Js - d11s
        k1 = Ks - d12s
        i2 = Is
        j2 = Js + d11s
        k2 = Ks + d12s
        i3 = Is
        j3 = Js - d21s
        k3 = Ks - d22s
        i4 = Is
        j4 = Js + d21s
        k4 = Ks + d22s
      case(3:4)
        i1 = Is - d11s
        j1 = Js 
        k1 = Ks - d12s
        i2 = Is + d11s
        j2 = Js 
        k2 = Ks + d12s
        i3 = Is - d21s
        j3 = Js
        k3 = Ks - d22s
        i4 = Is + d21s
        j4 = Js 
        k4 = Ks + d22s
      case(5:6)
        i1 = Is - d11s
        j1 = Js - d12s
        k1 = Ks 
        i2 = Is + d11s
        j2 = Js + d12s
        k2 = Ks
        i3 = Is - d21s
        j3 = Js - d22s
        k3 = Ks
        i4 = Is + d21s
        j4 = Js + d22s
        k4 = Ks
    end select

    Ts1 = ph2vars( blks%P(np,i1,j1,k1), blks%P(nh,i1,j1,k1), T_tab )
    Ts2 = ph2vars( blks%P(np,i2,j2,k2), blks%P(nh,i2,j2,k2), T_tab )
    Ts3 = ph2vars( blks%P(np,i3,j3,k3), blks%P(nh,i3,j3,k3), T_tab )
    Ts4 = ph2vars( blks%P(np,i4,j4,k4), blks%P(nh,i4,j4,k4), T_tab )

    Tg (3) = Ts1
    Tg (4) = Ts2
    Tg (5) = Ts3
    Tg (6) = Ts4

  end subroutine Fill_BC_Ghost_Coupling_FUSS


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


  subroutine Fill_BC_Ghost_Coupling_ARES ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, blks, Pg )
    implicit none
    integer, intent(in)               :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(FUSS_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Pg(1,6)
    ! Local
    integer  :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
    
    select case(Fs)
      case(1:2)
        i1 = Is
        j1 = Js - d11s
        k1 = Ks - d12s
        i2 = Is
        j2 = Js + d11s
        k2 = Ks + d12s
        i3 = Is
        j3 = Js - d21s
        k3 = Ks - d22s
        i4 = Is
        j4 = Js + d21s
        k4 = Ks + d22s
      case(3:4)
        i1 = Is - d11s
        j1 = Js 
        k1 = Ks - d12s
        i2 = Is + d11s
        j2 = Js 
        k2 = Ks + d12s
        i3 = Is - d21s
        j3 = Js
        k3 = Ks - d22s
        i4 = Is + d21s
        j4 = Js 
        k4 = Ks + d22s
      case(5:6)
        i1 = Is - d11s
        j1 = Js - d12s
        k1 = Ks 
        i2 = Is + d11s
        j2 = Js + d12s
        k2 = Ks
        i3 = Is - d21s
        j3 = Js - d22s
        k3 = Ks
        i4 = Is + d21s
        j4 = Js + d22s
        k4 = Ks
    end select

    Pg (1,3) = blks % T (i1,j1,k1)
    Pg (1,4) = blks % T (i2,j2,k2)
    Pg (1,5) = blks % T (i3,j3,k3)
    Pg (1,6) = blks % T (i4,j4,k4)

  end subroutine Fill_BC_Ghost_Coupling_ARES
    
end module hydraAF_Lib_Ghost