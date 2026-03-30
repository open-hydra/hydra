module hydraMM_Mod_Coupling
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Global_m, only: guide_
  use MOSKA_Types_m,  only:  MOSKA_config_type => config_type , MOSKA_1Ddomain_type, MOSKA_block_type
  use MOSE_Types_m ,  only: MOSE_config_type => config_type, MOSE_domain_type, MOSE_block_type
  use MOSE_Global_m

  implicit none
  
  contains

    subroutine Fix_Global_dt ( MOSKA_domain, MOSE_domain, MOSKAconfig, MOSEconfig )
        use MOSKA_Mod_dt, only: MOSKA_Compute_dt => Compute_dt1D, MOSKA_Set_dt => Set_Global_dt1D
        use MOSE_Mod_dt, only: MOSE_Compute_dt => Compute_dt, MOSE_Set_dt => Set_Global_dt

        implicit none
        type(MOSKA_1Ddomain_type), intent(inout)   :: MOSKA_domain
        type(MOSE_domain_type), intent(inout)    :: MOSE_domain
        type(MOSE_config_type), intent(inout)      :: MOSEconfig
        type(MOSKA_config_type), intent(inout)      :: MOSKAconfig
        ! Local
        real(R8) :: min_dt

        ! Compute Global dt
        MOSKA_domain % dtglobal = 1e5
        MOSE_domain % dtglobal = 1e5
        
        call MOSKA_Compute_dt( MOSKA_domain, MOSKAconfig)       
         
        call MOSE_Compute_dt( MOSE_domain, MOSEconfig)
        ! Assign Global dt
        min_dt = min(MOSKA_domain % dtglobal, MOSE_domain % dtglobal)

        MOSKA_domain % dtglobal = min_dt
        MOSE_domain % dtglobal = min_dt
        call MOSKA_Set_dt (MOSKA_domain)
        call MOSE_Set_dt (MOSE_domain)
       

    end subroutine Fix_Global_dt

    subroutine Inj_Area(MOSKA_domain, MOSE_domain)
            implicit none
            type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain
            type(MOSE_domain_type), intent(inout)    :: MOSE_domain
            integer :: b,f,upper,lower,Bm,Im,Jm,Km,Fm,i
            integer :: InjID, face_i,face_j, face_k, dir
            real(8) :: Area
       
            upper = 0
   
            do b = 1, MOSE_domain%nb
                do f = 1, 6
                    lower = upper + 1
                    upper = upper + MOSE_domain%n_bf(b,f)
                   !$omp parallel
                   !$omp do schedule(dynamic) &
                   !$omp private(i,Bm,Im,Jm,Km,Fm,InjID,face_i,face_j,face_k,dir,Area) 
                   
                    do i = lower, upper
                        select case ( MOSE_domain%bc(i)%type )
                            case (7)
                                Bm  = MOSE_domain % bc(i) % b
                                Im  = MOSE_domain % bc(i) % i
                                Jm  = MOSE_domain % bc(i) % j 
                                Km  = MOSE_domain % bc(i) % k 
                                Fm  = MOSE_domain % bc(i) % f
                                InjID  = MOSE_domain % bc(i) % bs
                               
                                
                            call Face_Index( Fm, dir, Im, Jm, Km, face_i, face_j, face_k )
                            !Compute Area of all the cells facing the 1D injector (The mesh has to be built accordingly, this is done to double-check)
                            !$omp critical   
                                Area = MOSE_domain%blk(Bm)%dir(dir)%f(face_i,face_j,face_k)%A
                                           
                                MOSKA_domain%Blk(InjID)%tmpA = MOSKA_domain%Blk(InjID)%tmpA + Area
                            !$omp end critical
                                

                        end select
                    enddo
                   
                    !$omp end do 
                     !$omp end parallel
        
                enddo
            enddo
     
    end subroutine Inj_Area

    subroutine Fill_MOSExMOSKA ( MOSKA_domain, MOSE_domain )
        implicit none
        type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain
        type(MOSE_domain_type), intent(inout)    :: MOSE_domain
        ! Local
        integer  :: b, f, lower, upper, i
        integer  :: Bm, Im, Jm, Km, Fm, Bs
        integer  ::  InjID

        upper = 0
        lower = 0
        do b = 1, MOSE_domain%nb
            do f = 1, 6
                
                lower = sum(MOSE_domain%n_bf(1:b-1, 1:6))+sum(MOSE_domain%n_bf(b,1:f-1)) + 1
                upper = sum(MOSE_domain%n_bf(b,1:f))+sum(MOSE_domain%n_bf(1:b-1, 1:6))
                
                !$omp do  schedule(dynamic) &
                !$omp private(i,Bm,Im,Jm,Km,Fm,InjID) 
              

                    do i = lower, upper
                        select case ( MOSE_domain%bc(i)%type )
                        case (7)
                            Bm  = MOSE_domain % bc(i) % b
                            Im  = MOSE_domain % bc(i) % i
                            Jm  = MOSE_domain % bc(i) % j 
                            Km  = MOSE_domain % bc(i) % k 
                            Fm  = MOSE_domain % bc(i) % f
                            InjID  = MOSE_domain % bc(i) % bs

                            !Filling multi-D domain ghost cell with 1D state and viceversa
                            call Hybrid_Coupling_Fill( Im, Jm, Km, Fm,  MOSE_domain % blk(Bm),&
                             MOSKA_domain % blk(InjID) )
                           
                
                
                        end select
                    enddo
                !$omp end  do     
            enddo
        enddo
    
     
    end subroutine Fill_MOSExMOSKA

    subroutine Flux_MOSExMOSKA ( MOSKA_domain, MOSE_domain )
        implicit none
        type(MOSKA_1Ddomain_type), intent(inout)    :: MOSKA_domain
        type(MOSE_domain_type), intent(inout)    :: MOSE_domain
        ! Local
        integer  :: b, f, lower, upper, i
        integer  :: Bm, Im, Jm, Km, Fm, Bs
        integer  ::  InjID


        upper = 0
        do b = 1, MOSE_domain%nb
            do f = 1, 6
                
                lower = sum(MOSE_domain%n_bf(1:b-1, 1:6))+sum(MOSE_domain%n_bf(b,1:f-1)) + 1
                upper = sum(MOSE_domain%n_bf(b,1:f))+sum(MOSE_domain%n_bf(1:b-1, 1:6))

         !$omp do schedule(dynamic) private(i,Bm,Im,Jm,Km,Fm,InjID,lower,upper)
            do i = lower, upper
                select case ( MOSE_domain%bc(i)%type )
            
                case (7)
                    Bm  = MOSE_domain % bc(i) % b
                    Im  = MOSE_domain % bc(i) % i
                    Jm  = MOSE_domain % bc(i) % j 
                    Km  = MOSE_domain % bc(i) % k 
                    Fm  = MOSE_domain % bc(i) % f
                    InjID  = MOSE_domain % bc(i) % bs
                    !Connection 1D-multi-D
                    call Hybrid_Connection_Eul( Im, Jm, Km, Fm, MOSE_domain % blk(Bm), MOSKA_domain % blk(InjID))
                    
                
                
                end select
            enddo
       
        !$omp end do
        enddo
        enddo
    

    end subroutine Flux_MOSExMOSKA
    
    subroutine Hybrid_Coupling_Fill( Im, Jm, Km, Fm, MOSE_Block, MOSKA_Block)
        implicit none
        integer, intent(in) :: Im,Jm,Km,Fm
        real(8)             :: tmp_val(nprim), Area
        type(MOSKA_block_type), intent(inout) ::  MOSKA_Block
        type(MOSE_block_type), intent(inout)  ::   MOSE_Block
        integer :: Ig,Jg,Kg,Fs,g, dir, face_i,face_j,face_k,kk
        
        call Face_Index( Fm, dir, Im, Jm, Km, face_i, face_j, face_k )

        Area = MOSE_Block%dir(dir)%f(face_i,face_j,face_k)%A
        

        select case (Fm)
            
            case(1) 
                
                Fs = 2

                do g = 1, gc
                    Ig = Im - guide_(Fm,1)*(g)
                    Jg = Jm - guide_(Fm,2)*(g)
                    Kg = Km - guide_(Fm,3)*(g)
                    MOSE_Block % P(1:nu,Ig,Jg,Kg) = MOSKA_Block % P (1:nu, MOSKA_Block%dim + guide(Fs,1)*(g-1) )
                    MOSE_Block % P(nv:nw,Ig,Jg,Kg) = 0.0
                    MOSE_Block % P(np:nprim,Ig,Jg,Kg) = MOSKA_Block % P (np-2:nprim-2, MOSKA_Block%dim + guide(Fs,1)*(g-1) )
                enddo

                do g = 1, gc
                    Ig = MOSKA_Block%dim - guide_(Fs,1)*(g)
                  
                     tmp_val = MOSE_Block%P(:, Im + guide_(Fm,1)*(g-1), Jm, Km) * (Area/MOSKA_Block%tmpA)
                    do kk = 1, nu
                   
                        !$omp atomic
                        MOSKA_Block % P(kk,Ig) = MOSKA_Block % P(kk,Ig) +  tmp_val(kk) 
                    enddo
                    !The -2 is added because index nu,np etc.. are from the MOSE domain which has two additional variables: v and w
                   ! do kk = np-2, np+npass-2
                    
                   !     !$omp atomic
                !        MOSKA_Block % P(kk,Ig) = MOSKA_Block % P(kk,Ig) + tmp_val(kk+2)
               !     enddo
        
                   
                enddo
               
      

        end select
             
    end subroutine Hybrid_Coupling_Fill

    subroutine Hybrid_Connection_Eul ( Im, Jm, Km, Fm, MOSE_Blk, MOSKA_Blk)
        use MOSE_Lib_Limiters
        use FLINT_Lib_Thermodynamic
        use MOSE_Lib_Riemann

        implicit none
        integer, intent(in) :: Im, Jm, Km, Fm
        type(MOSE_block_type), intent(inout) :: MOSE_Blk
        type(MOSKA_block_type), intent(inout) :: MOSKA_Blk
        ! Local
        integer :: modfm, modfm1, modfm2, modfm3, Dirm, Face_i, Face_j, Face_k, s
        integer :: I1, J1, K1, I2, J2, K2, Fs, II1 , II2
        integer :: I3, J3, K3, I4, J4, K4, kk, npass
        real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, diff21, diff32, diff43, Slope
        real(R8), dimension(nprim) :: Prim_L, Prim_R, Flux
        real(R8) :: Normal(3), Area, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
        real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_w, F_e, Fmass
        real(R8) :: su, sel_L, sel_R, tmp(nprim)


        call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
        call Face_Index ( Fm, Dirm, Im, Jm, Km, Face_i, Face_j, Face_k )

        ! Stencil construction: 4 points for 2nd order reconstruction from both sides at boundary face
        I1 = modfm*( Im - guide_(Fm,1) * gc) + modfm1*( Im + guide_(Fm,1) )
        J1 = modFm*( Jm - guide_(Fm,2) * gc) + modFm1*( Jm + guide_(Fm,2) )
        K1 = modFm*( Km - guide_(Fm,3) * gc) + modFm1*( Km + guide_(Fm,3) )

        I2 = modFm1*Im + modFm*( Im - guide_(Fm,1))
        J2 = modFm1*Jm + modFm*( Jm - guide_(Fm,2))
        K2 = modFm1*Km + modFm*( Km - guide_(Fm,3))


         select case (Fm)
            
            case(1) 
                
                Fs = 2

                    !Using the 1D domain geometry
                    II1 = MOSKA_Blk%dim + guide_(Fs,1)*(gc-1)
                    II2 = MOSKA_Blk%dim 
                
               
        

        end select

        I3 = modFm1*( Im - guide_(Fm,1)) + modFm*Im
        J3 = modFm1*( Jm - guide_(Fm,2)) + modFm*Jm
        K3 = modFm1*( Km - guide_(Fm,3)) + modFm*Km

        I4 = modFm1*( Im - guide_(Fm,1) * gc) + modFm*(Im + guide_(Fm,1))
        J4 = modFm1*( Jm - guide_(Fm,2) * gc) + modFm*(Jm + guide_(Fm,2))
        K4 = modFm1*( Km - guide_(Fm,3) * gc) + modFm*(Km + guide_(Fm,3))

        ! Primitive variables for the stencil
        
        Prim1 = MOSE_Blk % P(:,I1,J1,K1)
        Prim2 = MOSE_Blk % P(:,I2,J2,K2)
        Prim3 = MOSE_Blk % P(:,I3,J3,K3)
        Prim4 = MOSE_Blk % P(:,I4,J4,K4)
        
        Normal = MOSE_Blk % dir(Dirm) % f(Face_i,Face_j,Face_k) % n
        Area = MOSE_Blk % dir(Dirm) % f(Face_i,Face_j,Face_k) % a
        !Using the 1D domain geometry
        Dl1 = MOSKA_Blk % Dl(II1)  * 0.5d0
        Dl2 = MOSKA_Blk % Dl(II2)  * 0.5d0
        !Using the multi-D domain geometry
        Dl3 = MOSE_Blk % Dl(I3,J3,K3) % c(Dirm) * 0.5d0
        Dl4 = MOSE_Blk % Dl(I4,J4,K4) % c(Dirm) * 0.5d0
  
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
       
        call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), a_L, rho_L, &
                    Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), a_R, rho_R, &
                    normal(1),normal(2),normal(3), F_r, F_u, F_v, F_w, F_E)

        su = sign ( 0.5d0, F_r )
        Sel_L = 0.5d0 + su
        Sel_R = su - 0.5d0

        ! Fluxes
        Fmass = F_r * area
        Flux(1:nsc) = Fmass*( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
        Flux(nu) = F_u * Area
        Flux(nv) = F_v * Area
        Flux(nw) = F_w * Area
        Flux(np) = F_E * Area

        if (RANS) &
        Flux(nt:nprim) = Fmass*(sel_L*Prim_L(nt:nprim)/rho_L - sel_R*Prim_R(nt:nprim)/rho_R)

       ! if (npass>0) &
       ! Flux(np+1:np+npass) = Fmass*(sel_L*Prim_L(np+1:np+npass)/rho_L - sel_R*Prim_R(np+1:np+npass)/rho_R)
        
        ! Residual update
        MOSE_Blk % r(:,Im,Jm,Km) = MOSE_Blk % r(:,Im,Jm,Km) + Flux * modfm2
    
        if (Fm == 1) then
             tmp = Flux * modfm2
              do kk = 1, nu
                !$omp atomic
                MOSKA_Blk % r(kk,MOSKA_Blk%dim) = MOSKA_Blk % r(kk,MOSKA_Blk%dim) - tmp(kk)!- Flux(1:nu) * modfm2
              enddo
             do kk = np-2, nprim-2
                !$omp atomic
                MOSKA_Blk % r(kk,MOSKA_Blk%dim) = MOSKA_Blk % r(kk,MOSKA_Blk%dim) - tmp(kk+2) !- Flux(np:) * modfm2
             enddo

        endif

   
    end subroutine Hybrid_Connection_Eul

    subroutine Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

        implicit none
        integer, intent(in) :: fm
        integer, intent(out) :: modfm, modfm1, modfm2, modfm3

        modFm  = mod(Fm,2)
        modfm1 = 1-modFm  
        modfm2 = 1-2*modFm
        modfm3 = 2*( 2*modFm-1 )

    end subroutine Compute_Modfm

    subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
        
        implicit none
        integer, intent(in) :: face, cell_i, cell_j, cell_k
        integer, intent(out) :: dir, face_i, face_j, face_k

        select case ( face )
        
        case(1:2)
            dir = 1
            face_i = cell_i - mod ( face, 2 )
            face_j = cell_j
            face_k = cell_k
        
        case(3:4)
            dir = 2
            face_i = cell_i
            face_j = cell_j - mod ( face, 2 )
            face_k = cell_k
        
        case(5:6)
            dir = 3
            face_i = cell_i
            face_j = cell_j
            face_k = cell_k - mod ( face, 2 )
    
        end select

    end subroutine Face_Index

  end module hydraMM_Mod_Coupling
