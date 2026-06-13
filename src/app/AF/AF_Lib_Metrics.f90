module hydraAF_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m
  use ARES_Base_Types_m
  use ARES_Advanced_Types_m
  use FUSS_Base_Types_m
  use FUSS_Advanced_Types_m
  use hydraAF_Mod_CouplingExchange, only: is_FUSS_bc_local_local, is_ARES_bc_local_local

  implicit none
  private
  public :: Setup_Coupled_Metrics

contains

  subroutine Setup_Coupled_Metrics ( FUSS_domain, ARES_domain )
    ! NOTE: For MPI runs where a type(103) BC has its primary block (Bm) and
    ! secondary block (Bs) owned by different ranks, the Mg/dlg/volg for that
    ! entry are NOT populated here (the remote rank's node coords are needed).
    ! A follow-up metrics-node exchange in hydraAF_Mod_CouplingExchange is
    ! required for full cross-rank support.
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    ! Local
    integer :: i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s

    !$omp parallel
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, FUSS_domain % nbound
      select case ( FUSS_domain % bc(i) % type )
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
          d11s = FUSS_domain % bc(i) % d11
          d12s = FUSS_domain % bc(i) % d12
          d21s = FUSS_domain % bc(i) % d21
          d22s = FUSS_domain % bc(i) % d22
          call BC_Coupled_Metrics_FUSS ( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                         Is, Js, Ks, Fs, ARES_domain % blk(Bs), &
                                         d11s, d12s, d21s, d22s, FUSS_domain % bc(i) % Mg, &
                                         FUSS_domain % bc(i) % dlg, FUSS_domain % bc(i) % volg )
      end select
    enddo
    !$omp end parallel

    !$omp parallel
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, ARES_domain % nbound
      select case ( ARES_domain % bc(i) % type )
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
          d11s = ARES_domain % bc(i) % d11
          d12s = ARES_domain % bc(i) % d12
          d21s = ARES_domain % bc(i) % d21
          d22s = ARES_domain % bc(i) % d22
          call BC_Coupled_Metrics_ARES ( Im, Jm, Km, Fm, ARES_domain % blk(Bm), &
                                         Is, Js, Ks, Fs, FUSS_domain % blk(Bs), &
                                         d11s, d12s, d21s, d22s, ARES_domain % bc(i) % Mg, &
                                         ARES_domain % bc(i) % dlg, ARES_domain % bc(i) % volg )
      end select
    enddo
    !$omp end parallel

  end subroutine Setup_Coupled_Metrics


  subroutine Compute_Metric_Tensor_FUSS( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol )

    implicit none
    type(FUSS_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(FUSS_tensor_3D_type), intent(out) :: M
    type(FUSS_vector_3D_type), intent(out) :: dl
    real(R8), intent(out)             :: vol
    ! Local
    integer  :: h
    real(R8) :: det, A(3,3), cofactor(3,3)
    real(R8) :: vx(8), vy(8), vz(8)

    ! Compute average cell dimension. A is M^-1, inverse of metric tensor.
    ! M1 = [ xcs, ycs, zcs ; xet, yet, zet ; xzi, yzi, zzi ].
    
    A(1,:) = N5 % c - N1 % c + &
             N6 % c - N2 % c + &
             N7 % c - N3 % c + &
             N8 % c - N4 % c
  
    A(2,:) = N3 % c - N1 % c + &
             N4 % c - N2 % c + &
             N7 % c - N5 % c + &
             N8 % c - N6 % c
             
    A(3,:) = N2 % c - N1 % c + &
             N4 % c - N3 % c + &
             N6 % c - N5 % c + &
             N8 % c - N7 % c

    A = 0.25d0 * A

    ! Determinant of M^-1.
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    if ( abs(det) == 0d0 ) then
      write(*,'(A90)') ' WARNING - Metric tensor det=0. Should not happen, but going on with M==I'
      stop

      M % c = 0d0
      do h = 1, 3
        M % c(h,h) = 1d0
      enddo

    else
      
      cofactor(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))
      cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      cofactor(1,3) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))
      cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      cofactor(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))
      cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      cofactor(3,1) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))
      cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      cofactor(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))

      M % c = cofactor / det

    endif

    ! Cell length in i,j,k
    do h = 1, 3
      dl % c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
    enddo

    ! Cell volume computation
    vx(1) = N1%c(1)
    vy(1) = N1%c(2)
    vz(1) = N1%c(3)

    vx(2) = N5%c(1)
    vy(2) = N5%c(2)
    vz(2) = N5%c(3)

    vx(3) = N3%c(1)
    vy(3) = N3%c(2)
    vz(3) = N3%c(3)

    vx(4) = N7%c(1)
    vy(4) = N7%c(2)
    vz(4) = N7%c(3)
  
    vx(5) = N2%c(1)
    vy(5) = N2%c(2)
    vz(5) = N2%c(3)

    vx(6) = N6%c(1)
    vy(6) = N6%c(2)
    vz(6) = N6%c(3)

    vx(7) = N4%c(1)
    vy(7) = N4%c(2)
    vz(7) = N4%c(3)
  
    vx(8) = N8%c(1)
    vy(8) = N8%c(2)
    vz(8) = N8%c(3)
  
    vol = tvol(vx,vy,vz,1,2,3,5) + tvol(vx,vy,vz,2,4,3,8) &
        + tvol(vx,vy,vz,5,8,6,2) + tvol(vx,vy,vz,5,7,8,3) &
        + tvol(vx,vy,vz,5,8,2,3)
  
    if( vol <= 0d0 ) then
      write(*,*) 'Negative volume'
      stop
    endif
  

    contains
  
      pure function tvol(vx, vy, vz, i1, i2, i3, i4) result(volume)
        implicit none
        real(R8), intent(in) :: vx(8), vy(8), vz(8)
        integer, intent(in)  :: i1, i2, i3, i4
        real(R8)             :: volume
  
        volume = abs(((vx(i2)-vx(i1))* &
          ((vy(i3)-vy(i1))*(vz(i4)-vz(i1))-(vy(i4)-vy(i1))*(vz(i3)-vz(i1)))+ &
                                  (vy(i2)-vy(i1))* &
          ((vx(i4)-vx(i1))*(vz(i3)-vz(i1))-(vx(i3)-vx(i1))*(vz(i4)-vz(i1)))+ &
                                  (vz(i2)-vz(i1))* &
          ((vx(i3)-vx(i1))*(vy(i4)-vy(i1))-(vx(i4)-vx(i1))*(vy(i3)-vy(i1)))) &
          /6.d0)
  
      end function tvol
  
  end subroutine Compute_Metric_Tensor_FUSS


  subroutine Compute_Metric_Tensor_ARES( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol )

    implicit none
    type(ARES_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(ARES_tensor_3D_type), intent(out) :: M
    type(ARES_vector_3D_type), intent(out) :: dl
    real(R8), intent(out)             :: vol
    ! Local
    integer  :: h
    real(R8) :: det, A(3,3), cofactor(3,3)
    real(R8) :: vx(8), vy(8), vz(8)

    ! Compute average cell dimension. A is M^-1, inverse of metric tensor.
    ! M1 = [ xcs, ycs, zcs ; xet, yet, zet ; xzi, yzi, zzi ].
    
    A(1,:) = N5 % c - N1 % c + &
             N6 % c - N2 % c + &
             N7 % c - N3 % c + &
             N8 % c - N4 % c
  
    A(2,:) = N3 % c - N1 % c + &
             N4 % c - N2 % c + &
             N7 % c - N5 % c + &
             N8 % c - N6 % c
             
    A(3,:) = N2 % c - N1 % c + &
             N4 % c - N3 % c + &
             N6 % c - N5 % c + &
             N8 % c - N7 % c

    A = 0.25d0 * A

    ! Determinant of M^-1.
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    if ( abs(det) == 0d0 ) then
      write(*,'(A90)') ' WARNING - Metric tensor det=0. Should not happen, but going on with M==I'
      stop

      M % c = 0d0
      do h = 1, 3
        M % c(h,h) = 1d0
      enddo

    else
      
      cofactor(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))
      cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      cofactor(1,3) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))
      cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      cofactor(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))
      cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      cofactor(3,1) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))
      cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      cofactor(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))

      M % c = cofactor / det

    endif

    ! Cell length in i,j,k
    do h = 1, 3
      dl % c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
    enddo

    ! Cell volume computation
    vx(1) = N1%c(1)
    vy(1) = N1%c(2)
    vz(1) = N1%c(3)

    vx(2) = N5%c(1)
    vy(2) = N5%c(2)
    vz(2) = N5%c(3)

    vx(3) = N3%c(1)
    vy(3) = N3%c(2)
    vz(3) = N3%c(3)

    vx(4) = N7%c(1)
    vy(4) = N7%c(2)
    vz(4) = N7%c(3)
  
    vx(5) = N2%c(1)
    vy(5) = N2%c(2)
    vz(5) = N2%c(3)

    vx(6) = N6%c(1)
    vy(6) = N6%c(2)
    vz(6) = N6%c(3)

    vx(7) = N4%c(1)
    vy(7) = N4%c(2)
    vz(7) = N4%c(3)
  
    vx(8) = N8%c(1)
    vy(8) = N8%c(2)
    vz(8) = N8%c(3)
  
    vol = tvol(vx,vy,vz,1,2,3,5) + tvol(vx,vy,vz,2,4,3,8) &
        + tvol(vx,vy,vz,5,8,6,2) + tvol(vx,vy,vz,5,7,8,3) &
        + tvol(vx,vy,vz,5,8,2,3)
  
    if( vol <= 0d0 ) then
      write(*,*) 'Negative volume'
      stop
    endif
  

    contains
  
      pure function tvol(vx, vy, vz, i1, i2, i3, i4) result(volume)
        implicit none
        real(R8), intent(in) :: vx(8), vy(8), vz(8)
        integer, intent(in)  :: i1, i2, i3, i4
        real(R8)             :: volume
  
        volume = abs(((vx(i2)-vx(i1))* &
          ((vy(i3)-vy(i1))*(vz(i4)-vz(i1))-(vy(i4)-vy(i1))*(vz(i3)-vz(i1)))+ &
                                  (vy(i2)-vy(i1))* &
          ((vx(i4)-vx(i1))*(vz(i3)-vz(i1))-(vx(i3)-vx(i1))*(vz(i4)-vz(i1)))+ &
                                  (vz(i2)-vz(i1))* &
          ((vx(i3)-vx(i1))*(vy(i4)-vy(i1))-(vx(i4)-vx(i1))*(vy(i3)-vy(i1)))) &
          /6.d0)
  
      end function tvol
  
  end subroutine Compute_Metric_Tensor_ARES


  subroutine BC_Coupled_Metrics_FUSS ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg )

    implicit none
    integer, intent(in)                    :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(FUSS_block_type), intent(inout)   :: blkm
    type(ARES_block_type), intent(in)      :: blks
    type(FUSS_tensor_3D_type), intent(out) :: Mg(2)
    type(FUSS_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    ! Local
    integer      :: g, Is1, Js1, Ks1
    integer      :: II, JJ, KK, III, JJJ, KKK      
    integer      :: guidem(6,3), guides(6,3), guidem2(6,3), guides2(6,3)
    type(FUSS_vector_3D_type), dimension(0:gc_):: N1, N2, N3, N4, N5, N6, N7, N8
    
    ! ---------------------------------------------------------------------------------------------
    ! Preliminary definitions
    guidem(1,1)=-2 ; guidem(1,2)= 0 ; guidem(1,3)= 0 ; guidem(2,1)=1 ; guidem(2,2)=0 ; guidem(2,3)=0
    guidem(3,1)= 0 ; guidem(3,2)=-2 ; guidem(3,3)= 0 ; guidem(4,1)=0 ; guidem(4,2)=1 ; guidem(4,3)=0
    guidem(5,1)= 0 ; guidem(5,2)=0  ; guidem(5,3)=-2 ; guidem(6,1)=0 ; guidem(6,2)=0 ; guidem(6,3)=1

    guidem2(1,1)=-3 ; guidem2(1,2)= 0 ; guidem2(1,3)= 0 ; guidem2(2,1)=2 ; guidem2(2,2)=0 ; guidem2(2,3)=0
    guidem2(3,1)= 0 ; guidem2(3,2)=-3 ; guidem2(3,3)= 0 ; guidem2(4,1)=0 ; guidem2(4,2)=2 ; guidem2(4,3)=0
    guidem2(5,1)= 0 ; guidem2(5,2)= 0 ; guidem2(5,3)=-3 ; guidem2(6,1)=0 ; guidem2(6,2)=0 ; guidem2(6,3)=2

    guides(1,1)= 0 ; guides(1,2)= 0 ; guides(1,3)= 0 ; guides(2,1)=-1 ; guides(2,2)= 0 ; guides(2,3)= 0
    guides(3,1)= 0 ; guides(3,2)= 0 ; guides(3,3)= 0 ; guides(4,1)= 0 ; guides(4,2)=-1 ; guides(4,3)= 0
    guides(5,1)= 0 ; guides(5,2)= 0 ; guides(5,3)= 0 ; guides(6,1)= 0 ; guides(6,2)= 0 ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0 ; guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0 ; guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1 ; guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2
    ! ---------------------------------------------------------------------------------------------

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
      case(2)
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(3)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
      case(4)
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(5)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
      case(6)
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
    end select

    do g = 1, gc_

      ! connected cell node index
      Is1 = Is + (gc_-g)*guides(Fs,1)+(g-1)*guides2(Fs,1)
      Js1 = Js + (gc_-g)*guides(Fs,2)+(g-1)*guides2(Fs,2)
      Ks1 = Ks + (gc_-g)*guides(Fs,3)+(g-1)*guides2(Fs,3)

      ! compute the indexes in the connected block corresponding to the ghost nodes to update in the local block 
      select case(Fs)
        case(1:2) ! i-faces
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4) ! j-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6) ! k-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      if (Fm == 1) then
          N4(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 2) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 3) then
          N6(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 4) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 5) then
          N7(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 6) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      endif

      ! update the remaining 3 indexes (only necessary in the corners, but for good measure)
      ! node i1-1,i2 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s
          KKK=KK-d12s
        case(3:4)
          III=II-d11s
          JJJ=JJ
          KKK=KK-d12s
        case(5:6)
          III=II-d11s
          JJJ=JJ-d12s
          KKK=KK
      end select
      if (Fm == 1) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d21s
          KKK=KK-d22s
        case(3:4)
          III=II-d21s
          JJJ=JJ
          KKK=KK-d22s
        case(5:6)
          III=II-d21s
          JJJ=JJ-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1-1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s-d21s
          KKK=KK-d12s-d22s
        case(3:4)
          III=II-d11s-d21s
          JJJ=JJ
          KKK=KK-d12s-d22s
        case(5:6)
          III=II-d11s-d21s
          JJJ=JJ-d12s-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! compute metric variables
      select case(Fm)
        case(1)
          call Compute_Metric_Tensor_FUSS ( N1(g), N2(g), N3(g), N4(g), N1(g-1), N2(g-1), N3(g-1), N4(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(2)
          call Compute_Metric_Tensor_FUSS ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        case(3)
          call Compute_Metric_Tensor_FUSS ( N1(g), N2(g), N1(g-1), N2(g-1), N5(g), N6(g), N5(g-1), N6(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(4)
          call Compute_Metric_Tensor_FUSS ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        case(5)
          call Compute_Metric_Tensor_FUSS ( N1(g), N1(g-1), N3(g), N3(g-1), N5(g), N5(g-1), N7(g), N7(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(6)
          call Compute_Metric_Tensor_FUSS ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        end select

    enddo ! End loop over ghost cell           

  end subroutine BC_Coupled_Metrics_FUSS


  subroutine BC_Coupled_Metrics_ARES ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg )

    implicit none
    integer, intent(in)                    :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(ARES_block_type), intent(inout)   :: blkm
    type(FUSS_block_type), intent(in)      :: blks
    type(ARES_tensor_3D_type), intent(out) :: Mg(2)
    type(ARES_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    ! Local
    integer      :: g, Is1, Js1, Ks1
    integer      :: II, JJ, KK, III, JJJ, KKK      
    integer      :: guidem(6,3), guides(6,3), guidem2(6,3), guides2(6,3)
    type(ARES_vector_3D_type), dimension(0:gc_):: N1, N2, N3, N4, N5, N6, N7, N8
    
    ! ---------------------------------------------------------------------------------------------
    ! Preliminary definitions
    guidem(1,1)=-2 ; guidem(1,2)= 0 ; guidem(1,3)= 0 ; guidem(2,1)=1 ; guidem(2,2)=0 ; guidem(2,3)=0
    guidem(3,1)= 0 ; guidem(3,2)=-2 ; guidem(3,3)= 0 ; guidem(4,1)=0 ; guidem(4,2)=1 ; guidem(4,3)=0
    guidem(5,1)= 0 ; guidem(5,2)=0  ; guidem(5,3)=-2 ; guidem(6,1)=0 ; guidem(6,2)=0 ; guidem(6,3)=1

    guidem2(1,1)=-3 ; guidem2(1,2)= 0 ; guidem2(1,3)= 0 ; guidem2(2,1)=2 ; guidem2(2,2)=0 ; guidem2(2,3)=0
    guidem2(3,1)= 0 ; guidem2(3,2)=-3 ; guidem2(3,3)= 0 ; guidem2(4,1)=0 ; guidem2(4,2)=2 ; guidem2(4,3)=0
    guidem2(5,1)= 0 ; guidem2(5,2)= 0 ; guidem2(5,3)=-3 ; guidem2(6,1)=0 ; guidem2(6,2)=0 ; guidem2(6,3)=2

    guides(1,1)= 0 ; guides(1,2)= 0 ; guides(1,3)= 0 ; guides(2,1)=-1 ; guides(2,2)= 0 ; guides(2,3)= 0
    guides(3,1)= 0 ; guides(3,2)= 0 ; guides(3,3)= 0 ; guides(4,1)= 0 ; guides(4,2)=-1 ; guides(4,3)= 0
    guides(5,1)= 0 ; guides(5,2)= 0 ; guides(5,3)= 0 ; guides(6,1)= 0 ; guides(6,2)= 0 ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0 ; guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0 ; guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1 ; guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2
    ! ---------------------------------------------------------------------------------------------

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
      case(2)
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(3)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
      case(4)
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(5)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
      case(6)
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
    end select

    do g = 1, gc_

      ! connected cell node index
      Is1 = Is + (gc_-g)*guides(Fs,1)+(g-1)*guides2(Fs,1)
      Js1 = Js + (gc_-g)*guides(Fs,2)+(g-1)*guides2(Fs,2)
      Ks1 = Ks + (gc_-g)*guides(Fs,3)+(g-1)*guides2(Fs,3)

      ! compute the indexes in the connected block corresponding to the ghost nodes to update in the local block 
      select case(Fs)
        case(1:2) ! i-faces
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4) ! j-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6) ! k-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      if (Fm == 1) then
          N4(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 2) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 3) then
          N6(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 4) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 5) then
          N7(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 6) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      endif

      ! update the remaining 3 indexes (only necessary in the corners, but for good measure)
      ! node i1-1,i2 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s
          KKK=KK-d12s
        case(3:4)
          III=II-d11s
          JJJ=JJ
          KKK=KK-d12s
        case(5:6)
          III=II-d11s
          JJJ=JJ-d12s
          KKK=KK
      end select
      if (Fm == 1) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d21s
          KKK=KK-d22s
        case(3:4)
          III=II-d21s
          JJJ=JJ
          KKK=KK-d22s
        case(5:6)
          III=II-d21s
          JJJ=JJ-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1-1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s-d21s
          KKK=KK-d12s-d22s
        case(3:4)
          III=II-d11s-d21s
          JJJ=JJ
          KKK=KK-d12s-d22s
        case(5:6)
          III=II-d11s-d21s
          JJJ=JJ-d12s-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! compute metric variables
      select case(Fm)
        case(1)
          call Compute_Metric_Tensor_ARES ( N1(g), N2(g), N3(g), N4(g), N1(g-1), N2(g-1), N3(g-1), N4(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(2)
          call Compute_Metric_Tensor_ARES ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        case(3)
          call Compute_Metric_Tensor_ARES ( N1(g), N2(g), N1(g-1), N2(g-1), N5(g), N6(g), N5(g-1), N6(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(4)
          call Compute_Metric_Tensor_ARES ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        case(5)
          call Compute_Metric_Tensor_ARES ( N1(g), N1(g-1), N3(g), N3(g-1), N5(g), N5(g-1), N7(g), N7(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(6)
          call Compute_Metric_Tensor_ARES ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        end select

    enddo ! End loop over ghost cell           

  end subroutine BC_Coupled_Metrics_ARES


end module hydraAF_Lib_Metrics