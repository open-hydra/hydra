module hydraMQ2D_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m
  use Q2D_Types_m
  use MOSKA_types_m
  use hydraMQ2D_Mod_Coupling, only: Face_Index_1D, Face_Index, Compute_Modfm, &
                                    Q2D_to_MOSKA_map, thread_dlg

  implicit none
  private
  public :: Compute_Metrics_Coupled

contains

  subroutine Compute_Metrics_Coupled ( Q2D_domain, MOSKA_domain )

    implicit none
    type(Q2D_domain_type), intent(inout)      :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain

    !$omp parallel
    ! Initialize dlg values
    call Initialize_BC_Coupled_dlg ( Q2D_domain, MOSKA_domain )

    ! Fill ghost metrics
    call Fill_Metrics_Ghost_Cell ( Q2D_domain, MOSKA_domain )
    !$omp end parallel

  end subroutine Compute_Metrics_Coupled

  subroutine Initialize_BC_Coupled_dlg ( Q2D_domain, MOSKA_domain ) 

    implicit none
    type(Q2D_domain_type), intent(inout)        :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain
    ! Local
    integer  :: b, f, lower, upper, i, g
    
    upper = 0
    do b = 1, Q2D_domain%nb
      do f = 1, 4
        lower = upper + 1
        upper = upper + Q2D_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i, g)
        do i = lower, upper
          if ( Q2D_domain%bc(i)%type == 999 ) then
            do g = 1, gc_
            Q2D_domain % bc(i) % dlg(g) % c = 0.0d0
            enddo
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
            MOSKA_domain % bc(i) % dlg = 0.0d0
          endif
        enddo
      enddo
    enddo

  end subroutine Initialize_BC_Coupled_dlg

  subroutine Fill_Metrics_Ghost_Cell ( Q2D_domain, MOSKA_domain )
#if defined (_OPENMP)
    use omp_lib
#endif

    implicit none
    type(Q2D_domain_type), intent(inout) :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout) :: MOSKA_domain
    ! Local
    integer :: i, g, Bm, Im, Jm, Fm, Bs, Fs
    integer :: i_moska, tid, nthreads
    integer :: Fm_, Im_
    integer :: modfm_, modfm1_, modfm2_, modfm3_
    real(R8) :: dlg_local(gc_)
    real(R8) :: z_hydra_local

#if defined (_OPENMP)
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    ! Allocate thread-local dlg arrays (module variable with SAVE attribute)
    !$omp single
    if (.not. allocated(thread_dlg)) then
      allocate(thread_dlg(gc_, MOSKA_domain%nbound, nthreads))
    endif
    thread_dlg = 0.0d0
    !$omp end single

    !$omp do schedule (dynamic) private(i, g, Bm, Im, Jm, Fm, Bs, Fs, i_moska, tid) &
    !$omp private(dlg_local, z_hydra_local)
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

        ! Fill Q2D ghost cell metrics (no change needed)
        call Metrics_Coupling_Q2D ( Im, Jm, Fm, Q2D_domain % blk(Bm), &
                                    Fs, MOSKA_domain % blk(Bs), &
                                    Q2D_domain % bc(i) % dlg )

        ! *** OPTIMIZED: O(1) lookup instead of O(N) search ***
        i_moska = Q2D_to_MOSKA_map(i)%moska_bc_index

        if (i_moska > 0) then
          z_hydra_local = Q2D_to_MOSKA_map(i)%z_hydra

          call Metrics_Coupling_MOSKA ( Fs, MOSKA_domain % blk(Bs), &
                                        Im, Jm, Fm, Q2D_domain % blk(Bm), &
                                        z_hydra_local, dlg_local )

          ! *** OPTIMIZED: Thread-local accumulation (no critical section!) ***
          do g = 1, gc_
            thread_dlg(g, i_moska, tid) = thread_dlg(g, i_moska, tid) + dlg_local(g)
          enddo
        endif
      endif
    enddo

    ! Reduce thread-local arrays to final MOSKA dlg and update ghost cells
    !$omp do schedule(static) private(i, tid, g, Bs, Fm_, Im_, modfm_, modfm1_, modfm2_, modfm3_)
    do i = 1, MOSKA_domain%nbound
      if (MOSKA_domain%bc(i)%type == 999) then
        Bs = MOSKA_domain%bc(i)%b
        Fm_ = MOSKA_domain%bc(i)%f
        call Compute_Modfm ( Fm_, modfm_, modfm1_, modfm2_, modfm3_ )
        Im_ = modfm_ + modfm1_ * MOSKA_domain%blk(Bs)%dim

        do tid = 1, nthreads
          do g = 1, gc_
            MOSKA_domain%bc(i)%dlg(g) = MOSKA_domain%bc(i)%dlg(g) + thread_dlg(g, i, tid)
          enddo
        enddo

        ! Update ghost cell metrics in block
        do g = 1, gc_
          MOSKA_domain%blk(Bs)%dl(Im_ - guide_(Fm_,1)*g) = MOSKA_domain%bc(i)%dlg(g)
        enddo
      endif
    enddo

  end subroutine Fill_Metrics_Ghost_Cell


  subroutine Metrics_Coupling_Q2D ( Im, Jm, Fm, blkm, Fs, blks, dlg )

    implicit none
    integer, intent(in)                     :: Im, Jm, Fm, Fs
    type(Q2D_block_type), intent(inout)        :: blkm
    type(MOSKA_block_type), intent(in)      :: blks
    type(Q2D_vector_2D_type), intent(out)   :: dlg(2)
    ! Local
    integer  :: g, Ig, Jg, Is, Is_, face_i, face_j, dirm
    integer :: modfs, modfs1, modfs2, modfs3  

    call Compute_Modfm ( Fs, modfs, modfs1, modfs2, modfs3 )
    Is = modfs + modfs1 * blks%dim

    do g = 1, gc_

      call Face_Index( Fm, dirm, Im, Jm, face_i, face_j)
      Ig = Im - guide_(Fm,1)*(g)
      Jg = Jm - guide_(Fm,2)*(g)
  
      Is_ = Is + guide_(Fs,1)*(g-1) 
      
      blkm % dl(Ig,Jg) % c(dirm) = blks % dl(Is_) ! Non mi interessa dell'altra direzione
      dlg(g) % c(dirm) = blkm % dl(Ig,Jg) % c(dirm)
    enddo

  end subroutine Metrics_Coupling_Q2D


  subroutine Metrics_Coupling_MOSKA ( Fm, blkm, Is, Js, Fs, blks, z_hydra, dlg ) 
    implicit none
    integer, intent(in)                :: Fm, Is, Js, Fs
    type(MOSKA_block_type), intent(in) :: blkm
    type(Q2D_block_type), intent(in)   :: blks
    real(R8), intent(in)               :: z_hydra
    real(R8), intent(inout)            :: dlg(gc_)
    ! Local
    real(R8) :: Area_Q2D, Area_MOSKA
    integer  :: g, Im, Ig, Is_, Js_, face_i, face_j, dirs
    integer  :: modfm, modfm1, modfm2, modfm3


    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )

    Im = modfm + modfm1 * blkm%dim
    call Face_Index_1D( Fm, Im, face_i)
    Area_MOSKA = blkm%dir%f(face_i)%A

    call Face_Index( Fs, dirs, Is, Js, face_i, face_j)
    Area_Q2D = blks%dir(dirs)%f(face_i,face_j)%int_l * z_hydra

    do g = 1, gc_

      Ig = Im - guide_(Fm,1)*(g)

      Is_ = Is + guide_(Fs,1)*(g-1)
      Js_ = Js + guide_(Fs,2)*(g-1)
      
      dlg(g) = blks % dl(Is_,Js_) % c(dirs) * (Area_Q2D/Area_MOSKA)
      !blkm % dl(Ig) = dlg(g) 
    enddo

  end subroutine Metrics_Coupling_MOSKA

 

end module hydraMQ2D_Lib_Metrics
