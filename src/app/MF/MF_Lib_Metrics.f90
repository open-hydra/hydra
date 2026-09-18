module hydraMF_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m
  use MOSE_Base_Types_m
  use MOSE_Advanced_Types_m
  use FUSS_Base_Types_m
  use FUSS_Advanced_Types_m
  use hydraMF_Mod_CouplingExchange, only: is_FUSS_bc_local_local, is_MOSE_bc_local_local, &
                                          coupled_secondary_node_indices,               &
                                          exchange_coupling_metric_nodes,               &
                                          free_coupling_metric_nodes,                   &
                                          n_F_meta_entries, n_M_meta_entries,           &
                                          get_F_meta_entry, get_M_meta_entry

  implicit none
  private
  public :: Setup_Coupled_Metrics

contains

  !> Build Mg/dlg/volg for every type(103) coupling BC.
  !>
  !> The two loops below can only serve entries whose secondary block is owned
  !> by this rank, because they read that block's node coordinates directly.
  !> Entries whose two blocks were partitioned onto different ranks are then
  !> completed from coordinates fetched over MPI; skipping them leaves Mg = 0,
  !> which makes the wall-temperature routines drop the secondary side out of
  !> the interface temperature without any error being raised.
  subroutine Setup_Coupled_Metrics ( FUSS_domain, MOSE_domain, level )
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in) :: level
    ! Local
    integer :: i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    integer :: k
    real(R8) :: snode(3,4,gc_)

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
                                         Is, Js, Ks, Fs, MOSE_domain % blk(Bs), &
                                         d11s, d12s, d21s, d22s, FUSS_domain % bc(i) % Mg, &
                                         FUSS_domain % bc(i) % dlg, FUSS_domain % bc(i) % volg )
      end select
    enddo
    !$omp end parallel

    !$omp parallel
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, MOSE_domain % nbound
      select case ( MOSE_domain % bc(i) % type )
        case(103)
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
          d11s = MOSE_domain % bc(i) % d11
          d12s = MOSE_domain % bc(i) % d12
          d21s = MOSE_domain % bc(i) % d21
          d22s = MOSE_domain % bc(i) % d22
          call BC_Coupled_Metrics_MOSE ( Im, Jm, Km, Fm, MOSE_domain % blk(Bm), &
                                         Is, Js, Ks, Fs, FUSS_domain % blk(Bs), &
                                         d11s, d12s, d21s, d22s, MOSE_domain % bc(i) % Mg, &
                                         MOSE_domain % bc(i) % dlg, MOSE_domain % bc(i) % volg )
      end select
    enddo
    !$omp end parallel

    ! Cross-rank entries: fetch the secondary block's node coordinates from the
    ! rank that owns it, then build the same metric tensor from them.
    call exchange_coupling_metric_nodes ( FUSS_domain, MOSE_domain, level )

    do k = 1, n_F_meta_entries(level)
      call get_F_meta_entry ( level, k, i, snode )
      if (i <= 0) cycle
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
                                     Is, Js, Ks, Fs, MOSE_domain % blk(Bs), &
                                     d11s, d12s, d21s, d22s, FUSS_domain % bc(i) % Mg, &
                                     FUSS_domain % bc(i) % dlg, FUSS_domain % bc(i) % volg, &
                                     snode )
    enddo

    do k = 1, n_M_meta_entries(level)
      call get_M_meta_entry ( level, k, i, snode )
      if (i <= 0) cycle
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
      d11s = MOSE_domain % bc(i) % d11
      d12s = MOSE_domain % bc(i) % d12
      d21s = MOSE_domain % bc(i) % d21
      d22s = MOSE_domain % bc(i) % d22
      call BC_Coupled_Metrics_MOSE ( Im, Jm, Km, Fm, MOSE_domain % blk(Bm), &
                                     Is, Js, Ks, Fs, FUSS_domain % blk(Bs), &
                                     d11s, d12s, d21s, d22s, MOSE_domain % bc(i) % Mg, &
                                     MOSE_domain % bc(i) % dlg, MOSE_domain % bc(i) % volg, &
                                     snode )
    enddo

    call free_coupling_metric_nodes ( level )

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


  subroutine Compute_Metric_Tensor_MOSE( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol )

    implicit none
    type(MOSE_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(MOSE_tensor_3D_type), intent(out) :: M
    type(MOSE_vector_3D_type), intent(out) :: dl
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
  
  end subroutine Compute_Metric_Tensor_MOSE


  !> When snode is present the secondary block's node coordinates are taken
  !> from it instead of from blks, which lets a rank build Mg for an entry
  !> whose secondary block belongs to another rank. blks is then untouched.
  subroutine BC_Coupled_Metrics_FUSS ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg, snode )

    implicit none
    integer, intent(in)                    :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(FUSS_block_type), intent(inout)   :: blkm
    type(MOSE_block_type), intent(in)      :: blks
    type(FUSS_tensor_3D_type), intent(out) :: Mg(2)
    type(FUSS_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    real(R8), intent(in), optional         :: snode(:,:,:)  !< (3, node, ghost level)
    ! Local
    integer      :: g, n, idx(3,4)
    real(R8)     :: nc(3,4)
    type(FUSS_vector_3D_type), dimension(0:gc_):: N1, N2, N3, N4, N5, N6, N7, N8

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

      ! Node indices in the connected (secondary) block, then the four node
      ! coordinates themselves: read straight from that block when it is local,
      ! or taken from what its owning rank sent when it is not.
      call coupled_secondary_node_indices ( g, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, idx )
      if ( present(snode) ) then
        nc(:,1:4) = snode(:,1:4,g)
      else
        do n = 1, 4
          nc(:,n) = blks % node(idx(1,n),idx(2,n),idx(3,n)) % c
        enddo
      endif

      ! connected node
      if (Fm == 1) then
          N4(g)%c = nc(:,1)
      elseif (Fm == 2) then
          N8(g)%c = nc(:,1)
      elseif (Fm == 3) then
          N6(g)%c = nc(:,1)
      elseif (Fm == 4) then
          N8(g)%c = nc(:,1)
      elseif (Fm == 5) then
          N7(g)%c = nc(:,1)
      elseif (Fm == 6) then
          N8(g)%c = nc(:,1)
      endif

      ! the remaining 3 nodes (only necessary in the corners, but for good measure)
      ! node i1-1,i2
      if (Fm == 1) then
        N2(g)%c = nc(:,2)
      elseif (Fm == 2) then
        N6(g)%c = nc(:,2)
      elseif (Fm == 3) then
        N2(g)%c = nc(:,2)
      elseif (Fm == 4) then
        N4(g)%c = nc(:,2)
      elseif (Fm == 5) then
        N3(g)%c = nc(:,2)
      elseif (Fm == 6) then
        N4(g)%c = nc(:,2)
      endif

      ! node i1,i2-1
      if (Fm == 1) then
        N3(g)%c = nc(:,3)
      elseif (Fm == 2) then
        N7(g)%c = nc(:,3)
      elseif (Fm == 3) then
        N5(g)%c = nc(:,3)
      elseif (Fm == 4) then
        N7(g)%c = nc(:,3)
      elseif (Fm == 5) then
        N5(g)%c = nc(:,3)
      elseif (Fm == 6) then
        N6(g)%c = nc(:,3)
      endif

      ! node i1-1,i2-1
      if (Fm == 1) then
        N1(g)%c = nc(:,4)
      elseif (Fm == 2) then
        N5(g)%c = nc(:,4)
      elseif (Fm == 3) then
        N1(g)%c = nc(:,4)
      elseif (Fm == 4) then
        N3(g)%c = nc(:,4)
      elseif (Fm == 5) then
        N1(g)%c = nc(:,4)
      elseif (Fm == 6) then
        N2(g)%c = nc(:,4)
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


  !> When snode is present the secondary block's node coordinates are taken
  !> from it instead of from blks, which lets a rank build Mg for an entry
  !> whose secondary block belongs to another rank. blks is then untouched.
  subroutine BC_Coupled_Metrics_MOSE ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg, snode )

    implicit none
    integer, intent(in)                    :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(MOSE_block_type), intent(inout)   :: blkm
    type(FUSS_block_type), intent(in)      :: blks
    type(MOSE_tensor_3D_type), intent(out) :: Mg(2)
    type(MOSE_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    real(R8), intent(in), optional         :: snode(:,:,:)  !< (3, node, ghost level)
    ! Local
    integer      :: g, n, idx(3,4)
    real(R8)     :: nc(3,4)
    type(MOSE_vector_3D_type), dimension(0:gc_):: N1, N2, N3, N4, N5, N6, N7, N8

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

      ! Node indices in the connected (secondary) block, then the four node
      ! coordinates themselves: read straight from that block when it is local,
      ! or taken from what its owning rank sent when it is not.
      call coupled_secondary_node_indices ( g, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, idx )
      if ( present(snode) ) then
        nc(:,1:4) = snode(:,1:4,g)
      else
        do n = 1, 4
          nc(:,n) = blks % node(idx(1,n),idx(2,n),idx(3,n)) % c
        enddo
      endif

      ! connected node
      if (Fm == 1) then
          N4(g)%c = nc(:,1)
      elseif (Fm == 2) then
          N8(g)%c = nc(:,1)
      elseif (Fm == 3) then
          N6(g)%c = nc(:,1)
      elseif (Fm == 4) then
          N8(g)%c = nc(:,1)
      elseif (Fm == 5) then
          N7(g)%c = nc(:,1)
      elseif (Fm == 6) then
          N8(g)%c = nc(:,1)
      endif

      ! the remaining 3 nodes (only necessary in the corners, but for good measure)
      ! node i1-1,i2
      if (Fm == 1) then
        N2(g)%c = nc(:,2)
      elseif (Fm == 2) then
        N6(g)%c = nc(:,2)
      elseif (Fm == 3) then
        N2(g)%c = nc(:,2)
      elseif (Fm == 4) then
        N4(g)%c = nc(:,2)
      elseif (Fm == 5) then
        N3(g)%c = nc(:,2)
      elseif (Fm == 6) then
        N4(g)%c = nc(:,2)
      endif

      ! node i1,i2-1
      if (Fm == 1) then
        N3(g)%c = nc(:,3)
      elseif (Fm == 2) then
        N7(g)%c = nc(:,3)
      elseif (Fm == 3) then
        N5(g)%c = nc(:,3)
      elseif (Fm == 4) then
        N7(g)%c = nc(:,3)
      elseif (Fm == 5) then
        N5(g)%c = nc(:,3)
      elseif (Fm == 6) then
        N6(g)%c = nc(:,3)
      endif

      ! node i1-1,i2-1
      if (Fm == 1) then
        N1(g)%c = nc(:,4)
      elseif (Fm == 2) then
        N5(g)%c = nc(:,4)
      elseif (Fm == 3) then
        N1(g)%c = nc(:,4)
      elseif (Fm == 4) then
        N3(g)%c = nc(:,4)
      elseif (Fm == 5) then
        N1(g)%c = nc(:,4)
      elseif (Fm == 6) then
        N2(g)%c = nc(:,4)
      endif

      ! compute metric variables
      select case(Fm)
        case(1)
          call Compute_Metric_Tensor_MOSE ( N1(g), N2(g), N3(g), N4(g), N1(g-1), N2(g-1), N3(g-1), N4(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(2)
          call Compute_Metric_Tensor_MOSE ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        case(3)
          call Compute_Metric_Tensor_MOSE ( N1(g), N2(g), N1(g-1), N2(g-1), N5(g), N6(g), N5(g-1), N6(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(4)
          call Compute_Metric_Tensor_MOSE ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        case(5)
          call Compute_Metric_Tensor_MOSE ( N1(g), N1(g-1), N3(g), N3(g-1), N5(g), N5(g-1), N7(g), N7(g-1), &
                                            Mg(g), dlg(g), volg(g))
        case(6)
          call Compute_Metric_Tensor_MOSE ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                            Mg(g), dlg(g), volg(g))
        end select

    enddo ! End loop over ghost cell           

  end subroutine BC_Coupled_Metrics_MOSE


end module hydraMF_Lib_Metrics