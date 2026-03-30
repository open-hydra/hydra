module hydraMQ2D_Mod_Coupling
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m, only: guide_
  use Q2D_Types_m, only: Q2D_domain_type, Q2D_config_type => config_type, Q2D_block_type
  use MOSKA_Types_m, only: MOSKA_1Ddomain_type, MOSKA_config_type => config_type

  implicit none
  private
  public :: Set_Global_dt, BC_Coupled_Fluxes, Face_Index_1D, Face_Index, Compute_Modfm
  public :: Build_Coupling_Lookup_Table, Free_Coupling_Lookup_Table
  public :: Q2D_to_MOSKA_map, thread_Pg, thread_dlg

  ! Lookup table for O(1) Q2D->MOSKA boundary mapping
  type :: bc_lookup_type
    integer :: moska_bc_index    ! Index of corresponding MOSKA boundary (-1 if not coupled)
    integer :: moska_face        ! MOSKA face index
    real(R8) :: z_hydra          ! Area scaling factor
  end type bc_lookup_type

  type(bc_lookup_type), allocatable, save :: Q2D_to_MOSKA_map(:)

  ! Thread-local flux accumulation array for parallel coupling
  real(R8), allocatable, save :: thread_flux(:,:,:)  ! (nprim, nbound_MOSKA, nthreads)

  ! Thread-local Pg accumulation array for parallel extrapolation
  real(R8), allocatable, save :: thread_Pg(:,:,:,:)  ! (nprim_M, gc_, nbound_MOSKA, nthreads)

  ! Thread-local dlg accumulation array for parallel metrics
  real(R8), allocatable, save :: thread_dlg(:,:,:)   ! (gc_, nbound_MOSKA, nthreads)

contains

  subroutine Set_Global_dt ( Q2D_domain, MOSKA_domain, Q2D_config, MOSKA_config ) 
    use Q2D_Mod_dt, only: Q2D_Compute_dt => Compute_dt, Q2D_Set_dt => Set_Global_dt
    use MOSKA_Mod_dt, only: MOSKA_Compute_dt => Compute_dt1D, MOSKA_Set_dt => Set_Global_dt1D

    implicit none
    type(Q2D_domain_type), intent(inout)       :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout)   :: MOSKA_domain
    type(Q2D_config_type), intent(inout)       :: Q2D_config
    type(MOSKA_config_type), intent(inout)     :: MOSKA_config
    ! Local
    real(R8) :: min_dt

    ! Compute Global dt
    Q2D_domain % dtglobal = 1e5
    MOSKA_domain % dtglobal = 1e5
    call Q2D_Compute_dt( Q2D_domain, Q2D_config )
    call MOSKA_Compute_dt( MOSKA_domain, MOSKA_config )

    ! Assign Global dt
    min_dt = min(Q2D_domain % dtglobal, MOSKA_domain % dtglobal)
    Q2D_domain % dtglobal = min_dt
    MOSKA_domain % dtglobal = min_dt

    call Q2D_Set_dt ( Q2D_domain ) 
    call MOSKA_Set_dt ( MOSKA_domain )

  end subroutine Set_Global_dt


  subroutine BC_Coupled_Fluxes ( Q2D_domain, MOSKA_domain )

    implicit none
    type(Q2D_domain_type), intent(inout)        :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain

    ! Initialize ext_flux variables
    call Initialize_BC_Coupled_Fluxes ( Q2D_domain, MOSKA_domain )

    ! Compute BC fluxes
    call Compute_Coupled_Fluxes ( Q2D_domain, MOSKA_domain )

  end subroutine BC_Coupled_Fluxes



  subroutine Initialize_BC_Coupled_Fluxes ( Q2D_domain, MOSKA_domain ) !FATTO

    implicit none
    type(Q2D_domain_type), intent(inout)        :: Q2D_domain
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
            Q2D_domain % bc(i) % ext_flux = 0.0d0
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
            MOSKA_domain % bc(i) % ext_flux = 0.0d0
          endif
        enddo
      enddo
    enddo

  end subroutine Initialize_BC_Coupled_Fluxes


  subroutine Compute_Coupled_Fluxes ( Q2D_domain, MOSKA_domain )
    use MOSKA_Global_m, only: nprim_M => nprim
    use Q2D_Global_m, only: nprim_Q => nprim
#if defined (_OPENMP)
    use omp_lib
#endif

    implicit none
    type(Q2D_domain_type),   intent(inout)    :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain
    ! Local
    integer  :: b, f, lower, upper, i, j, i_moska, tid, nthreads
    integer  :: Bm, Im, Jm, Fm, Bs, Fs
    real(R8) :: Q2D_new_flux(nprim_Q)
    real(R8) :: MOSKA_new_flux(nprim_M)
    real(R8) :: z_hydra_local

#if defined (_OPENMP)
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    ! Allocate thread-local flux arrays (module variable with SAVE attribute)
    ! Use single to ensure only one thread allocates
    !$omp single
    if (.not. allocated(thread_flux)) then
      allocate(thread_flux(nprim_M, MOSKA_domain%nbound, nthreads))
    endif
    thread_flux = 0.0d0
    !$omp end single

    upper = 0
    do b = 1, Q2D_domain%nb
      do f = 1, 4
        lower = upper + 1
        upper = upper + Q2D_domain%n_bf(b,f)
        !$omp do schedule ( guided, 8 ) private (i, j, Bm, Im, Jm, Fm, Bs, Fs, i_moska, tid) &
        !$omp private (Q2D_new_flux, MOSKA_new_flux, z_hydra_local)
        do i = lower, upper
          if (Q2D_domain%bc(i)%type == 999) then
#if defined (_OPENMP)
              tid = omp_get_thread_num() + 1
#else
              tid = 1
#endif
              Bm  = Q2D_domain % bc(i) % b
              Im  = Q2D_domain % bc(i) % i
              Jm  = Q2D_domain % bc(i) % j
              Fm  = Q2D_domain % bc(i) % f
              Bs  = Q2D_domain % bc(i) % bs
              Fs  = Q2D_domain % bc(i) % fs

              ! Compute Q2D flux (no change)
              call Q2D_Moska_Connection_Eul( Im, Jm, Fm, Q2D_domain % blk(Bm), Q2D_new_flux )
              Q2D_domain % bc(i) % ext_flux = Q2D_domain % bc(i) % ext_flux + Q2D_new_flux

              ! *** OPTIMIZED: O(1) lookup instead of O(N) search ***
              i_moska = Q2D_to_MOSKA_map(i)%moska_bc_index

              if (i_moska > 0) then
                ! Found matching MOSKA boundary via lookup table
                ! Cache lookup table values for better memory access
                z_hydra_local = Q2D_to_MOSKA_map(i)%z_hydra
                call Moska_Q2D_Connection_Eul ( Im, Jm, Fm, Q2D_domain % blk(Bm), &
                                               Fs, &
                                               MOSKA_domain % blk(Bs) % nsign_global, MOSKA_new_flux)

                ! *** OPTIMIZED: Thread-local accumulation (no atomics!) ***
                thread_flux(:, i_moska, tid) = thread_flux(:, i_moska, tid) + &
                                               MOSKA_new_flux(:) * z_hydra_local
              endif

          endif
        enddo
      enddo
    enddo
    ! Reduce thread-local arrays to final MOSKA flux
    ! Loop order optimized for cache locality: innermost loop over j (contiguous in thread_flux)
    !$omp do schedule(static) private(i, tid, j)
    do i = 1, MOSKA_domain%nbound
      if (MOSKA_domain%bc(i)%type == 999) then
        do tid = 1, nthreads
          do j = 1, nprim_M
            MOSKA_domain % bc(i) % ext_flux(j) = MOSKA_domain % bc(i) % ext_flux(j) + thread_flux(j, i, tid)
          enddo
        enddo
      endif
    enddo

    ! Note: thread_flux is not deallocated here as it's a module variable
    ! It will be deallocated in Free_Coupling_Lookup_Table

  end subroutine Compute_Coupled_Fluxes


  subroutine Q2D_Moska_Connection_Eul ( Im, Jm, Fm, blkm, Q2D_flux)
    use Q2D_Global_m
    use Q2D_Lib_Limiters
    use Q2D_Lib_Riemann
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Fm
    type(Q2D_block_type), intent(inout)        :: blkm 
    real(R8), intent(inout)                    :: Q2D_flux(nprim)
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dirm, Face_i, Face_j, s
    integer :: I1, J1, I2, J2
    integer :: I3, J3, I4, J4
    real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, diff21, diff32, diff43, Slope
    real(R8), dimension(nprim) :: Prim_L, Prim_R, Flux
    real(R8) :: Normal(2), int_l, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
    real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_e, Fmass
    real(R8) :: su, sel_L, sel_R


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dirm, Im, Jm, Face_i, Face_j )

    ! Stencil construction: 4 points for 2nd order reconstruction from both sides at boundary face
    I1 = modfm*( Im - guide_(Fm,1) * gc) + modfm1*( Im + guide_(Fm,1) )
    J1 = modFm*( Jm - guide_(Fm,2) * gc) + modFm1*( Jm + guide_(Fm,2) )

    I2 = modFm1*Im + modFm*( Im - guide_(Fm,1))
    J2 = modFm1*Jm + modFm*( Jm - guide_(Fm,2))

    I3 = modFm1*( Im - guide_(Fm,1)) + modFm*Im
    J3 = modFm1*( Jm - guide_(Fm,2)) + modFm*Jm

    I4 = modFm1*( Im - guide_(Fm,1) * gc) + modFm*(Im + guide_(Fm,1))
    J4 = modFm1*( Jm - guide_(Fm,2) * gc) + modFm*(Jm + guide_(Fm,2))

    ! Primitive variables for the stencil
    Prim1 = blkm % P(:,I1,J1)
    Prim2 = blkm % P(:,I2,J2)
    Prim3 = blkm % P(:,I3,J3)
    Prim4 = blkm % P(:,I4,J4)
    
    Normal = blkm % dir(Dirm) % f(Face_i,Face_j) % n
    int_l = blkm % dir(Dirm) % f(Face_i,Face_j) % int_l
    Dl1 = blkm % Dl(I1,J1) % c(Dirm) * 0.5d0
    Dl2 = blkm % Dl(I2,J2) % c(Dirm) * 0.5d0
    Dl3 = blkm % Dl(I3,J3) % c(Dirm) * 0.5d0
    Dl4 = blkm % Dl(I4,J4) % c(Dirm) * 0.5d0

    Dist23 = Dl2 + Dl3
    Dist12 = Dl1 + Dl2
    Dist34 = Dl3 + Dl4

    ! Reconstruction in point (2)
    diff32 = ( Prim3 - Prim2 ) / Dist23
    diff21 = ( Prim2 - Prim1 ) / Dist12

    ! Limiter for point (2) slope => state (1) for the riemann solver
    do s = 1, nprim
      Slope(s) = rlimiter ( diff32(s), diff21(s) )
    end do
    
    ! Extrapolation from point (2) to the interface => state (1)
    Prim_L = Prim2 + Slope * Dl2

    ! Auxiliary variables for state (1)
    call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rtot_L )
    a_L = f_ss ( Prim_L(1:nsc), Prim_L(np), rho_L, Rtot_L )
    
    ! Reconstruction in point (3)
    diff43 = ( Prim4 - Prim3 ) / Dist34

    ! Limiter for point (3) slope => state (4) for the riemann solver
    do s = 1, nprim
      Slope(s) = rlimiter ( diff43(s), diff32(s) )
    end do

    ! Extrapolation from point (3) to the interface => state (4)
    Prim_R = Prim3 - Slope * Dl3

    ! Auxiliary variables for state (4)
    call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rtot_R)
    a_R = f_ss ( Prim_R(1:nsc), Prim_R(np), rho_R, Rtot_R )
    
    call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(np), a_L, rho_L, &
                  Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(np), a_R, rho_R, &
                  normal(1),normal(2), F_r, F_u, F_v, F_E)

    su = sign ( 0.5d0, F_r )
    Sel_L = 0.5d0 + su
    Sel_R = su - 0.5d0

    ! Fluxes
    Fmass = F_r * int_l
    Flux(1:nsc) = Fmass*( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
    Flux(nu) = abs(guide_(Fm,1)) * F_u * int_l
    Flux(nv) = abs(guide_(Fm,2)) * F_v * int_l 
    Flux(np) = F_E * int_l

    ! Per turbolenza, per ora lasciamo cosi, trasporto e ho messo nella ghost il valore estrapolato
    if (RANS) &
    Flux(nt:nprim) = Fmass*(sel_L*Prim_L(nt:nprim)/rho_L - sel_R*Prim_R(nt:nprim)/rho_R)

    Q2D_flux = modfm2 * Flux
    !MOSKA_flux = -modfm2 * Flux * z_hydra ! N.B. non generale, funziona solo nel caso di connessione dispari/pari, pari/dispari

  end subroutine Q2D_Moska_Connection_Eul


  subroutine Moska_Q2D_Connection_Eul ( Im, Jm, Fm, blkm, Fs, sign_n, MOSKA_flux)
    use Q2D_Global_m
    use Q2D_Lib_Limiters
    use Q2D_Lib_Riemann
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in)                        :: Im, Jm, Fm, Fs, sign_n
    type(Q2D_block_type), intent(inout)        :: blkm 
    real(R8), intent(inout)                    :: MOSKA_flux(np-1)
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dirm, Face_i, Face_j, s
    integer :: modfs, modfs1, modfs2, modfs3
    integer :: I1, J1, I2, J2
    integer :: I3, J3, I4, J4
    real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, diff21, diff32, diff43, Slope
    real(R8), dimension(nprim) :: Prim_L, Prim_R
    real(R8), dimension(np-1) :: Flux
    real(R8) :: Normal(2), int_l, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
    real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_e, Fmass
    real(R8) :: su, sel_L, sel_R
    integer  :: n_sign_loc, Fs_, Im_, Jm_ ! sign_n che usi in estrapolazione

    ! Qua lavori sempre con blocco Q2D, ma devi stare molto attento a prendere le normali
    ! come se stessi lavorando in MOSKA; devi quindi vedere quando Fm e Fs sono entrambi pari
    ! o dispari, in tal caso le normali sono opposte rispetto a quelle di Q2D e di conseguenza si inverte
    ! la costruzione del problema di Riemann

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Compute_Modfm ( Fs, modfs, modfs1, modfs2, modfs3 )
    call Face_Index ( Fm, Dirm, Im, Jm, Face_i, Face_j )

    n_sign_loc = -modfs2 * modfm2 ! se sono entrambi pari o dispari n_sign = -1, altrimenti 1
    Normal = n_sign_loc * blkm % dir(Dirm) % f(Face_i,Face_j) % n
    int_l = blkm % dir(Dirm) % f(Face_i,Face_j) % int_l

    ! Ridefinisco Im e Jm mettendomi nel dominio MOSKA
    Im_ = Im - guide_(Fm,1)
    Jm_ = Jm - guide_(Fm,2)
    Fs_ = Fs + 2 * abs(guide_(Fm,2)) ! Se Fm è 3 o 4, incremento Fs di 2: Fs = 1 --> 3, Fs = 2 --> 4

    ! Stencil construction: 4 points for 2nd order reconstruction from both sides at boundary face
    I1 = ((1+n_sign_loc)/2) * (modfs*( Im_ - guide_(Fs_,1) * gc) + modfs1*( Im_ + guide_(Fs_,1) )) + &
          ((1-n_sign_loc)/2) * (modfs*( Im_ + guide_(Fs_,1) * gc) + modfs1*( Im_ - guide_(Fs_,1) ))
    J1 = ((1+n_sign_loc)/2) * (modfs*( Jm_ - guide_(Fs_,2) * gc) + modfs1*( Jm_ + guide_(Fs_,2) )) + &
          ((1-n_sign_loc)/2) * (modfs*( Jm_ + guide_(Fs_,2) * gc) + modfs1*( Jm_ - guide_(Fs_,2) ))

    I2 = ((1+n_sign_loc)/2) * (modfs1*Im_ + modfs*( Im_ - guide_(Fs_,1))) + &
          ((1-n_sign_loc)/2) * (modfs1*Im_ + modfs*( Im_ + guide_(Fs_,1)))
    J2 = ((1+n_sign_loc)/2) * (modfs1*Jm_ + modfs*( Jm_ - guide_(Fs_,2))) + &
          ((1-n_sign_loc)/2) * (modfs1*Jm_ + modfs*( Jm_ + guide_(Fs_,2)))

    I3 = ((1+n_sign_loc)/2) * (modfs1*( Im_ - guide_(Fs_,1)) + modfs*Im_) + &
          ((1-n_sign_loc)/2) * (modfs1*( Im_ + guide_(Fs_,1)) + modfs*Im_)
    J3 = ((1+n_sign_loc)/2) * (modfs1*( Jm_ - guide_(Fs_,2)) + modfs*Jm_) + &
          ((1-n_sign_loc)/2) * (modfs1*( Jm_ + guide_(Fs_,2)) + modfs*Jm_)

    I4 = ((1+n_sign_loc)/2) * (modfs1*( Im_ - guide_(Fs_,1) * gc) + modfs*(Im_ + guide_(Fs_,1))) + &
          ((1-n_sign_loc)/2) * (modfs1*( Im_ + guide_(Fs_,1) * gc) + modfs*(Im_ - guide_(Fs_,1)))
    J4 = ((1+n_sign_loc)/2) * (modfs1*( Jm_ - guide_(Fs_,2) * gc) + modfs*(Jm_ + guide_(Fs_,2))) + &
          ((1-n_sign_loc)/2) * (modfs1*( Jm_ + guide_(Fs_,2) * gc) + modfs*(Jm_ - guide_(Fs_,2)))      

    ! Primitive variables for the stencil
    Prim1 = blkm % P(:,I1,J1)
    Prim2 = blkm % P(:,I2,J2)
    Prim3 = blkm % P(:,I3,J3)
    Prim4 = blkm % P(:,I4,J4)
    

    Dl1 = blkm % Dl(I1,J1) % c(Dirm) * 0.5d0
    Dl2 = blkm % Dl(I2,J2) % c(Dirm) * 0.5d0
    Dl3 = blkm % Dl(I3,J3) % c(Dirm) * 0.5d0
    Dl4 = blkm % Dl(I4,J4) % c(Dirm) * 0.5d0

    Dist23 = Dl2 + Dl3
    Dist12 = Dl1 + Dl2
    Dist34 = Dl3 + Dl4

    ! Reconstruction in point (2)
    diff32 = ( Prim3 - Prim2 ) / Dist23
    diff21 = ( Prim2 - Prim1 ) / Dist12

    ! Limiter for point (2) slope => state (1) for the riemann solver
    do s = 1, nprim
      Slope(s) = rlimiter ( diff32(s), diff21(s) )
    end do
    
    ! Extrapolation from point (2) to the interface => state (1)
    Prim_L = Prim2 + Slope * Dl2

    ! Auxiliary variables for state (1)
    call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rtot_L )
    a_L = f_ss ( Prim_L(1:nsc), Prim_L(np), rho_L, Rtot_L )
    
    ! Reconstruction in point (3)
    diff43 = ( Prim4 - Prim3 ) / Dist34

    ! Limiter for point (3) slope => state (4) for the riemann solver
    do s = 1, nprim
      Slope(s) = rlimiter ( diff43(s), diff32(s) )
    end do

    ! Extrapolation from point (3) to the interface => state (4)
    Prim_R = Prim3 - Slope * Dl3

    ! Auxiliary variables for state (4)
    call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rtot_R)
    a_R = f_ss ( Prim_R(1:nsc), Prim_R(np), rho_R, Rtot_R )
    
    call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(np), a_L, rho_L, &
                  Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(np), a_R, rho_R, &
                  normal(1),normal(2), F_r, F_u, F_v, F_E)

    su = sign ( 0.5d0, F_r )
    Sel_L = 0.5d0 + su
    Sel_R = su - 0.5d0

    ! Fluxes
    Fmass = F_r * int_l
    Flux(1:nsc) = Fmass*( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
    Flux(nu) = sign_n * (abs(guide_(Fm,1))*F_u + abs(guide_(Fm,2))*F_v) * int_l
    Flux(np-1) = F_E * int_l

    MOSKA_flux = modfs2 * Flux 

  end subroutine Moska_Q2D_Connection_Eul

  ! subroutine Moska_Q2D_Connection_Eul ( Im, Jm, Fm, blkm, Is, Fs, blks, Moska_flux)
  !     use Q2D_Lib_Limiters
  !     use MOSE_Lib_Thermodynamic_Laws
  !     use MOSE_Lib_Riemann

  !     implicit none
  !     integer, intent(in) :: Im, Jm, Fm, Is, Fs
  !     type(Q2D_block_type), intent(inout)   :: blkm ! Q2D block
  !     type(MOSKA_block_type), intent(inout) :: blks ! Connected MOSKA block
  !     real(R8), intent(inout)                    :: Moska_flux(nprim)
  !     ! Local
  !     integer :: modfm, modfm1, modfm2, modfm3, Dirm, Face_i, Face_j, s
  !     integer :: I1, J1, I2, J2
  !     integer :: I3, J3, I4, J4
  !     real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, diff21, diff32, diff43, Slope
  !     real(R8), dimension(nprim) :: Prim_L, Prim_R, Flux
  !     real(R8) :: Normal(2), Area, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
  !     real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_e, Fmass
  !     real(R8) :: su, sel_L, sel_R


  !     call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
  !     call Compute_Modfm ( Fs, modfs, modfs1, modfs2, modfs3 )

  !     call Face_Index ( Fm, Dirm, Im, Jm, Face_i, Face_j )

  !     ! Stencil construction: 4 points for 2nd order reconstruction from both sides at boundary face
  !     I1 = modfm*( Im - guide_(Fm,1) * gc) + modfm1*( Im + guide_(Fm,1) )
  !     J1 = modFm*( Jm - guide_(Fm,2) * gc) + modFm1*( Jm + guide_(Fm,2) )

  !     I2 = modFm1*Im + modFm*( Im - guide_(Fm,1))
  !     J2 = modFm1*Jm + modFm*( Jm - guide_(Fm,2))

  !     I3 = modFm1*( Im - guide_(Fm,1)) + modFm*Im
  !     J3 = modFm1*( Jm - guide_(Fm,2)) + modFm*Jm

  !     I4 = modFm1*( Im - guide_(Fm,1) * gc) + modFm*(Im + guide_(Fm,1))
  !     J4 = modFm1*( Jm - guide_(Fm,2) * gc) + modFm*(Jm + guide_(Fm,2))

  !     ! Primitive variables for the stencil
  !     Prim1 = blkm % P(:,I1,J1)
  !     Prim2 = blkm % P(:,I2,J2)
  !     Prim3 = blkm % P(:,I3,J3)
  !     Prim4 = blkm % P(:,I4,J4)
      
  !     Normal = blkm % dir(Dirm) % f(Face_i,Face_j) % n
  !     Area = blkm % dir(Dirm) % f(Face_i,Face_j) % int_l
  !     Dl1 = blkm % Dl(I1,J1) % c(Dirm) * 0.5d0
  !     Dl2 = blkm % Dl(I2,J2) % c(Dirm) * 0.5d0
  !     Dl3 = blkm % Dl(I3,J3) % c(Dirm) * 0.5d0
  !     Dl4 = blkm % Dl(I4,J4) % c(Dirm) * 0.5d0

  !     Dist23 = Dl2 + Dl3
  !     Dist12 = Dl1 + Dl2
  !     Dist34 = Dl3 + Dl4

  !     ! Reconstruction in point (2)
  !     diff32 = ( Prim3 - Prim2 ) / Dist23
  !     diff21 = ( Prim2 - Prim1 ) / Dist12

  !     ! Limiter for point (2) slope => state (1) for the riemann solver
  !     do s = 1, nprim
  !       Slope(s) = rlimiter ( diff32(s), diff21(s) )
  !     end do
      
  !     ! Extrapolation from point (2) to the interface => state (1)
  !     Prim_L = Prim2 + Slope * Dl2

  !     ! Auxiliary variables for state (1)
  !     call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rtot_L )
  !     a_L = co_a ( Prim_L(1:nsc), Prim_L(np), rho_L, Rtot_L )
      
  !     ! Reconstruction in point (3)
  !     diff43 = ( Prim4 - Prim3 ) / Dist34

  !     ! Limiter for point (3) slope => state (4) for the riemann solver
  !     do s = 1, nprim
  !       Slope(s) = rlimiter ( diff43(s), diff32(s) )
  !     end do

  !     ! Extrapolation from point (3) to the interface => state (4)
  !     Prim_R = Prim3 - Slope * Dl3

  !     ! Auxiliary variables for state (4)
  !     call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rtot_R)
  !     a_R = co_a ( Prim_R(1:nsc), Prim_R(np), rho_R, Rtot_R )
      
  !     call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), a_L, rho_L, &
  !                 Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), a_R, rho_R, &
  !                 normal(1),normal(2),normal(3), F_r, F_u, F_v, F_w, F_E)

  !     su = sign ( 0.5d0, F_r )
  !     Sel_L = 0.5d0 + su
  !     Sel_R = su - 0.5d0

  !     ! Fluxes
  !     Fmass = F_r * int_l
  !     Flux(1:nsc) = Fmass*( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
  !     Flux(nu) = F_u * int_l
  !     Flux(np) = F_E * int_l
  !     Flux = - modfs2 * Flux ! Attento che deve essere modfs2, orientamento del blocco bs di MOSKA

  !     Moska_flux = Flux + Moska_flux

  ! end subroutine Q2D_Moska_Connection_Eul

  subroutine Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )

      implicit none
      integer, intent(in) :: Fm
      integer, intent(out) :: modfm, modfm1, modfm2, modfm3
  
      modFm  = mod(Fm,2)
      modfm1 = 1-modFm  
      modfm2 = 1-2*modFm
      modfm3 = 2*( 2*modFm-1 )
  
  end subroutine
  
    
  subroutine Face_Index ( face, dir, cell_i, cell_j, face_i, face_j )
  
      implicit none
      integer, intent(in) :: face, cell_i, cell_j
      integer, intent(out) :: dir, face_i, face_j
  
      select case ( face )
        case(1:2)
          dir = 1
          face_i = cell_i - mod ( face, 2 )
          face_j = cell_j
        case(3:4)
          dir = 2
          face_i = cell_i
          face_j = cell_j - mod ( face, 2 )
      end select
  
  end subroutine Face_Index

  subroutine Face_Index_1D ( face, cell_i, face_i)

      implicit none
      integer, intent(in) :: face, cell_i
      integer, intent(out) :: face_i

      select case ( face )
        case(1:2)
          face_i = cell_i - mod ( face, 2 )
      end select

  end subroutine Face_Index_1D


  !> Build lookup table for Q2D->MOSKA boundary mapping (call once during initialization)
  subroutine Build_Coupling_Lookup_Table ( Q2D_domain, MOSKA_domain )

    implicit none
    type(Q2D_domain_type), intent(in)       :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(in)   :: MOSKA_domain
    ! Local
    integer :: i, i_, b_, f_, lower_, upper_
    integer :: Bs, Fs, Bm_, Fm_
    logical :: found

    ! Allocate lookup table
    if (allocated(Q2D_to_MOSKA_map)) deallocate(Q2D_to_MOSKA_map)
    allocate(Q2D_to_MOSKA_map(Q2D_domain%nbound))

    ! Initialize all entries to "not coupled"
    do i = 1, Q2D_domain%nbound
      Q2D_to_MOSKA_map(i)%moska_bc_index = -1
      Q2D_to_MOSKA_map(i)%moska_face = -1
      Q2D_to_MOSKA_map(i)%z_hydra = 0.0d0
    enddo

    ! Build mapping for coupled boundaries (type 999)
    do i = 1, Q2D_domain%nbound
      if (Q2D_domain%bc(i)%type == 999) then
        Bs = Q2D_domain%bc(i)%bs
        Fs = Q2D_domain%bc(i)%fs

        ! Search for matching MOSKA boundary (do this ONCE during initialization)
        found = .false.
        upper_ = 0
        do b_ = 1, MOSKA_domain%nb
          do f_ = 1, 2
            lower_ = upper_ + 1
            upper_ = upper_ + MOSKA_domain%n_bf(b_,f_)
            do i_ = lower_, upper_
              if (MOSKA_domain%bc(i_)%type == 999) then
                Bm_ = MOSKA_domain%bc(i_)%b
                Fm_ = MOSKA_domain%bc(i_)%f
                if (Bs == Bm_ .and. Fs == Fm_) then
                  ! Found matching MOSKA boundary!
                  Q2D_to_MOSKA_map(i)%moska_bc_index = i_
                  Q2D_to_MOSKA_map(i)%moska_face = Fm_
                  Q2D_to_MOSKA_map(i)%z_hydra = MOSKA_domain%bc(i_)%z_hydra
                  found = .true.
                  exit
                endif
              endif
            enddo
            if (found) exit
          enddo
          if (found) exit
        enddo

        if (.not. found) then
          write(*,*) 'WARNING: Q2D boundary', i, 'with type 999 has no matching MOSKA boundary!'
          write(*,*) '         Bs =', Bs, ', Fs =', Fs
        endif
      endif
    enddo

  end subroutine Build_Coupling_Lookup_Table


  !> Free lookup table and thread-local arrays
  subroutine Free_Coupling_Lookup_Table()
    if (allocated(Q2D_to_MOSKA_map)) deallocate(Q2D_to_MOSKA_map)
    if (allocated(thread_flux)) deallocate(thread_flux)
    if (allocated(thread_Pg)) deallocate(thread_Pg)
    if (allocated(thread_dlg)) deallocate(thread_dlg)
  end subroutine Free_Coupling_Lookup_Table


end module hydraMQ2D_Mod_Coupling
