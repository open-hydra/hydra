module hydraMQ2D_Lib_Extrapolation
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m
  use Q2D_Types_m
  use MOSKA_types_m
  use hydraMQ2D_Mod_Coupling, only: Face_Index_1D, Face_Index, Compute_Modfm, &
                                    Q2D_to_MOSKA_map, thread_Pg

  implicit none
  private
  public :: Compute_Extrapolation_Coupled

contains

  subroutine Compute_Extrapolation_Coupled ( Q2D_domain, MOSKA_domain )

    implicit none
    type(Q2D_domain_type), intent(inout)      :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain

    ! Initialize MOSKA Pg values
    call Initialize_BC_Coupled_Pg ( Q2D_domain, MOSKA_domain )

    ! Fill Ghost Cells
    call Fill_Coupled_Ghost_Cell ( Q2D_domain, MOSKA_domain )

  end subroutine Compute_Extrapolation_Coupled


  subroutine Initialize_BC_Coupled_Pg ( Q2D_domain, MOSKA_domain ) 

    implicit none
    type(Q2D_domain_type), intent(inout)     :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain
    ! Local
    integer  :: b, f, lower, upper, i
    
    upper = 0
    do b = 1, Q2D_domain%nb
      do f = 1, 4
        lower = upper + 1
        upper = upper + Q2D_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i)
        do i = lower, upper
          if ( Q2D_domain%bc(i)%type == 999 ) then
            Q2D_domain % bc(i) % Pg = 0.0d0
          endif
        enddo
      enddo
    enddo

    upper = 0
    do b = 1, MOSKA_domain%nb
      do f = 1, 2
        lower = upper + 1
        upper = upper + MOSKA_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i)
        do i = lower, upper
          if ( MOSKA_domain%bc(i)%type == 999 ) then
            MOSKA_domain % bc(i) % Pg = 0.0d0
          endif
        enddo
      enddo
    enddo

  end subroutine Initialize_BC_Coupled_Pg


  subroutine Fill_Coupled_Ghost_Cell ( Q2D_domain, MOSKA_domain )
    use MOSKA_Global_m, only: nprim_M => nprim
#if defined (_OPENMP)
    use omp_lib
#endif

    implicit none
    type(Q2D_domain_type), intent(inout)   :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout) :: MOSKA_domain
    ! Local
    integer               :: i, j, g, Bm, Im, Jm, Fm, Bs, Fs
    integer               :: i_moska, tid, nthreads
    integer               :: Fm_, Im_
    integer               :: modfm_, modfm1_, modfm2_, modfm3_
    real(R8)              :: dPg_local(nprim_M, gc_)
    real(R8)              :: z_hydra_local

#if defined (_OPENMP)
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    ! Allocate thread-local Pg arrays (module variable with SAVE attribute)
    !$omp single
    if (.not. allocated(thread_Pg)) then
      allocate(thread_Pg(nprim_M, gc_, MOSKA_domain%nbound, nthreads))
    endif
    thread_Pg = 0.0d0
    !$omp end single

    !$omp do schedule (dynamic) private(i, j, g, Bm, Im, Jm, Fm, Bs, Fs, i_moska, tid) &
    !$omp private(dPg_local, z_hydra_local)
    do i = 1, Q2D_domain % nbound
      if ( Q2D_domain % bc(i) % type == 999 ) then
#if defined (_OPENMP)
        tid = omp_get_thread_num() + 1
#else
        tid = 1
#endif
        Bm = Q2D_domain % bc(i) % b
        Im = Q2D_domain % bc(i) % i
        Jm = Q2D_domain % bc(i) % j
        Fm = Q2D_domain % bc(i) % f
        Bs = Q2D_domain % bc(i) % bs
        Fs = Q2D_domain % bc(i) % fs

        ! Fill Q2D ghost cells (no change needed)
        call Ghost_Coupling_Q2D ( Im, Jm, Fm, Q2D_domain % blk(Bm), &
                                  Fs, MOSKA_domain % blk(Bs) % nsign_global, MOSKA_domain % blk(Bs), &
                                  Q2D_domain % bc(i) % Pg )

        ! *** OPTIMIZED: O(1) lookup instead of O(N) search ***
        i_moska = Q2D_to_MOSKA_map(i)%moska_bc_index

        if (i_moska > 0) then
          z_hydra_local = Q2D_to_MOSKA_map(i)%z_hydra

          call Ghost_Coupling_MOSKA ( Fs, MOSKA_domain % blk(Bs), MOSKA_domain % blk(Bs) % nsign_global, &
                                      Im, Jm, Fm, Q2D_domain % blk(Bm), &
                                      z_hydra_local, dPg_local )

          ! *** OPTIMIZED: Thread-local accumulation (no critical section!) ***
          do g = 1, gc_
            do j = 1, nprim_M
              thread_Pg(j, g, i_moska, tid) = thread_Pg(j, g, i_moska, tid) + dPg_local(j, g)
            enddo
          enddo
        endif
      endif
    enddo

    ! Reduce thread-local arrays to final MOSKA Pg and update ghost cells
    !$omp do schedule(static) private(i, tid, j, g, Bs, Fm_, Im_, modfm_, modfm1_, modfm2_, modfm3_)
    do i = 1, MOSKA_domain%nbound
      if (MOSKA_domain%bc(i)%type == 999) then
        Bs = MOSKA_domain%bc(i)%b
        Fm_ = MOSKA_domain%bc(i)%f
        call Compute_Modfm ( Fm_, modfm_, modfm1_, modfm2_, modfm3_ )
        Im_ = modfm_ + modfm1_ * MOSKA_domain%blk(Bs)%dim

        do tid = 1, nthreads
          do g = 1, gc_
            do j = 1, nprim_M
              MOSKA_domain%bc(i)%Pg(j,g) = MOSKA_domain%bc(i)%Pg(j,g) + thread_Pg(j, g, i, tid)
            enddo
          enddo
        enddo

        ! Update ghost cells in block
        do g = 1, gc_
          do j = 1, nprim_M
            MOSKA_domain%blk(Bs)%P(j, Im_ - guide_(Fm_,1)*g) = MOSKA_domain%bc(i)%Pg(j,g)
          enddo
        enddo
      endif
    enddo

  end subroutine Fill_Coupled_Ghost_Cell


  subroutine Ghost_Coupling_Q2D ( Im, Jm, Fm, blkm, Fs, sign_n, blks, Pg ) 
    use MOSKA_Global_m, only: nsc_M => nsc, nu_M => nu, np_M => np, nprim_M => nprim
    use Q2D_Global_m, only: nsc_Q => nsc, nu_Q => nu, nv_Q => nv, np_Q => np, nprim_Q => nprim, nt_Q => nt, nrans_Q => nrans

    implicit none
    integer, intent(in)                :: Im, Jm, Fm, Fs, sign_n
    type(Q2D_block_type), intent(inout)   :: blkm
    type(MOSKA_block_type), intent(in) :: blks
    real(R8), intent(inout)            :: Pg(nprim_Q,4)
    ! Local
    integer  :: g, Ig, Jg, Is, Is_
    integer  :: modfs, modfs1, modfs2, modfs3 

    call Compute_Modfm ( Fs, modfs, modfs1, modfs2, modfs3 )
    Is = modfs + modfs1 * blks%dim

    do g = 1, gc_

      Ig = Im - guide_(Fm,1)*(g)
      Jg = Jm - guide_(Fm,2)*(g)
  
      Is_ = Is + guide_(Fs,1)*(g-1)
           
      Pg(1:nsc_Q,g) = blks%P(1:nsc_M,Is_) ! copio specie da MOSKA
      Pg(nu_Q,g)    = abs(guide_(Fm,1)) * blks%P(nu_M,Is_) * sign_n   ! copio velocita' da MOSKA, per ogni iniettore puo essere una componente diversa 
      Pg(nv_Q,g)    = abs(guide_(Fm,2)) * blks%P(nu_M,Is_) * sign_n   ! v = 0 da MOSKA --> RIVEDI LA TRATTAZIONE DEI FLUSSI
      Pg(np_Q,g)    = blks%P(np_M,Is_)    ! copio pressione da MOSKA
        if (nrans_Q>0) then
            Pg(nt_Q:nprim_Q,g) = blkm%P(nt_Q:nprim_Q,Im,Jm) ! copio variabili turbolente dalla cella Im,Jm di Q2D
        end if

      blkm % P(:,Ig,Jg) = Pg(:,g) ! assegno alle ghost del blocco, mi servono per ricostruzione
    enddo

  end subroutine Ghost_Coupling_Q2D
  

  subroutine Ghost_Coupling_MOSKA ( Fm, blkm, sign_n, Is, Js, Fs, blks, z_hydra, Pg) 
    use MOSKA_Global_m, only: nsc_M => nsc, nu_M => nu, np_M => np, nprim_M => nprim ! sign_n lo calcoli in MOSKA, negativo se la normale è negativa
    use Q2D_Global_m, only: nsc_Q => nsc, nu_Q => nu, nv_Q => nv, np_Q => np, nprim_Q => nprim, nt_Q => nt, nrans_Q => nrans

    implicit none
    integer, intent(in)                :: Fm, Is, Js, Fs, sign_n
    type(MOSKA_block_type), intent(in) :: blkm
    type(Q2D_block_type), intent(in)   :: blks
    real(R8), intent(in)               :: z_hydra
    real(R8), intent(inout)            :: Pg(nprim_M,gc_)
    ! Local
    real(R8) :: Area_Q2D, Area_MOSKA

    integer  :: g, Im, Ig, Is_, Js_, face_i, face_j, dirs
    integer  :: modfm, modfm1, modfm2, modfm3


    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )

    Im = modfm + modfm1 * blkm%dim
    call Face_Index_1D( Fm, Im, face_i)
    Area_MOSKA = blkm%dir%f(face_i)%A

    call Face_Index( Fs, dirs, Is, Js, face_i, face_j)
    Area_Q2D = blks%dir(dirs)%f(face_i,face_j)%int_l * z_hydra ! Devo avere z_hydra, magari in global

    do g = 1, gc_

      Ig = Im - guide_(Fm,1)*(g)

      Is_ = Is + guide_(Fs,1)*(g-1)
      Js_ = Js + guide_(Fs,2)*(g-1)
      
      ! Fondamentale inizializzare a 0 Pg di MOSKA, altrimenti sommo a valori spuri
      Pg(1:nsc_M,g) = blks % P (1:nsc_Q,Is_,Js_) * (Area_Q2D/Area_MOSKA) ! sommo specie da Q2D pesate per l'area
      Pg(nu_M,g)    = (abs(guide_(Fs,1)) * blks % P (nu_Q,Is_,Js_) * sign_n + &
                      abs(guide_(Fs,2)) * blks % P (nv_Q,Is_,Js_) * sign_n) * (Area_Q2D/Area_MOSKA)
      Pg(np_M,g)    = blks % P (np_Q,Is_,Js_) * (Area_Q2D/Area_MOSKA)

      !blkm % P(:,Ig) = Pg(:,g) ! assegno alle ghost del blocco, mi servono per ricostruzione

    enddo

  end subroutine Ghost_Coupling_MOSKA


  ! subroutine Fill_BC_Ghost_Coupling_MOSE ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, blks, Pg )

  !   implicit none
  !   integer, intent(in)               :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
  !   type(FUSS_block_type), intent(in) :: blks
  !   real(R8), intent(inout)           :: Pg(1,6)
  !   ! Local
  !   integer  :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
    
  !   select case(Fs)
  !     case(1:2)
  !       i1 = Is
  !       j1 = Js - d11s
  !       k1 = Ks - d12s
  !       i2 = Is
  !       j2 = Js + d11s
  !       k2 = Ks + d12s
  !       i3 = Is
  !       j3 = Js - d21s
  !       k3 = Ks - d22s
  !       i4 = Is
  !       j4 = Js + d21s
  !       k4 = Ks + d22s
  !     case(3:4)
  !       i1 = Is - d11s
  !       j1 = Js 
  !       k1 = Ks - d12s
  !       i2 = Is + d11s
  !       j2 = Js 
  !       k2 = Ks + d12s
  !       i3 = Is - d21s
  !       j3 = Js
  !       k3 = Ks - d22s
  !       i4 = Is + d21s
  !       j4 = Js 
  !       k4 = Ks + d22s
  !     case(5:6)
  !       i1 = Is - d11s
  !       j1 = Js - d12s
  !       k1 = Ks 
  !       i2 = Is + d11s
  !       j2 = Js + d12s
  !       k2 = Ks
  !       i3 = Is - d21s
  !       j3 = Js - d22s
  !       k3 = Ks
  !       i4 = Is + d21s
  !       j4 = Js + d22s
  !       k4 = Ks
  !   end select

  !   Pg (1,3) = blks % T (i1,j1,k1)
  !   Pg (1,4) = blks % T (i2,j2,k2)
  !   Pg (1,5) = blks % T (i3,j3,k3)
  !   Pg (1,6) = blks % T (i4,j4,k4)

  ! end subroutine Fill_BC_Ghost_Coupling_MOSE    

end module hydraMQ2D_Lib_Extrapolation
