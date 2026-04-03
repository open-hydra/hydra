module hydraMF_Lib_BC_Coupled_Wall_Temperature
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: FUSS_BC_Coupled_Wall_Temperature, MOSE_BC_Coupled_Wall_Temperature
  
contains

  subroutine FUSS_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, Blk_FUSS, Is, Js, Ks, Fs, Blk_MOSE, Mg, Tg, new_flux, Ovar )
    use hydra_Parameters_m, only: guide_
    use FUSS_Advanced_Types_m
    use FUSS_Base_Types_m
    use FUSS_Lib_Solid
    use FUSS_Lib_BC_Fluxes, only: Compute_Wall_Properties
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    integer,intent(in)                     :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(FUSS_block_type), intent(inout)   :: Blk_FUSS
    type(MOSE_block_type), intent(inout)   :: Blk_MOSE
    type(FUSS_tensor_3D_type), intent(in)  :: Mg
    real(R8), intent(in)                   :: Tg(6)
    real(R8), intent(inout)                :: new_flux
    real(R8), optional, dimension(2), intent(out) :: Ovar

    ! Local
    integer  :: Dir, Face_im, Face_jm, Face_km, modfm, modfm1, modfm2, modfm3, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, M_loc(3,3), M_con(3,3)
    real(R8) :: T_loc, T_ip_loc, T_im_loc, T_jp_loc, T_jm_loc, T_kp_loc, T_km_loc, kappa_loc
    real(R8) :: rho, Rgas, mil
    real(R8) :: T_con, kappa_con
    real(R8) :: Num_loc, Num_con, Den_loc, Den_con
    real(R8) :: Tint, T_wall, q_wall, T, kappa, Gradient(3), Flux
  
    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )

    Ig = Im - guide_(Fm,1)
    Jg = Jm - guide_(Fm,2)
    Kg = Km - guide_(Fm,3)

    ! Metric stuff
    Normal = Blk_FUSS % dir(Dir) % f(Face_im,Face_jm,Face_km) % N
    Area   = Blk_FUSS % dir(Dir) % f(Face_im,Face_jm,Face_km) % A
    M_loc  = Blk_FUSS % M (Im,Jm,Km) % c 
    M_con  = Mg % c

    ! Stencil building
    T_loc    = Blk_FUSS % T (Im,Jm,Km)
    call co_Kappa( Blk_FUSS % matID(Im,Jm,Km), T_loc, kappa_loc )
    T_ip_loc = Blk_FUSS % T (Im+1,Jm,Km)
    T_im_loc = Blk_FUSS % T (Im-1,Jm,Km)
    T_jp_loc = Blk_FUSS % T (Im,Jm+1,Km)
    T_jm_loc = Blk_FUSS % T (Im,Jm-1,Km)
    if (ndir == 2) then
      T_kp_loc = huge(0.0)
      T_km_loc = huge(0.0)
    else
      T_kp_loc = Blk_FUSS % T (Im,Jm,Km+1)
      T_km_loc = Blk_FUSS % T (Im,Jm,Km-1)
    endif

    T_con = Tg(1)
    call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Is,Js,Ks), rho, Rgas )
    call co_k_mi_lam_Wilke ( Blk_MOSE % P(1:nsc,Is,Js,Ks), rho, T_con, mil, kappa_con )

    ! Tint computations
    select case ( Fm )
      case ( 1 : 2 )
        Den_loc = 2.0d0 * dot_product( M_loc(1,:), Normal ) * modfm2
        Num_loc = 2.0d0*T_loc * dot_product( M_loc(1,:), Normal ) * modfm2 &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(1,:), Normal ) * modfm2
        Num_con = 2.0d0*T_con * dot_product( M_con(1,:), Normal ) * modfm2 &
                + 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(2,:), Normal ) &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(3,:), Normal )
      case ( 3 : 4 )
        Den_loc = 2.0d0 * dot_product( M_loc(2,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(2,:), Normal ) * modfm2 &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(2,:), Normal ) * modfm2
        Num_con = 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(1,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(2,:), Normal ) * modfm2 &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(3,:), Normal )
      case ( 5 : 6 )
        Den_loc = 2.0d0 * dot_product( M_loc(3,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(3,:), Normal ) * modfm2
        Den_con = 2.0d0 * dot_product( M_con(3,:), Normal ) * modfm2
        Num_con = 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(1,:), Normal ) &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(2,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(3,:), Normal ) * modfm2
    end select

    Tint = (kappa_loc*Num_loc + kappa_con*Num_con) / (kappa_loc*Den_loc + kappa_con*Den_con)

    ! Calcolo dei flussi come nel caso di Temperature imposta Tint = Twall
    T_wall = Tint

    ! boundary cell variables
    T = Blk_FUSS % T (Im,Jm,Km)
    call co_Kappa( Blk_FUSS%matID(Im,Jm,Km), T, kappa )

    ! Difference in Across-face direction: symmetry; tangential direction: null
    Gradient = 0d0
    Gradient(Dir) = ( T - T_wall ) * modfm3
    if (ndir == 2) Gradient(3) = 0.0
    Gradient = matmul ( Gradient, M_loc )

    Flux = Area * kappa * dot_product ( Gradient, Normal )
    
    ! Residual update
    new_flux = - Flux * modfm2
    
    ! Wall variables
    q_wall = kappa * dot_product( Gradient, Normal )
    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine FUSS_BC_Coupled_Wall_Temperature


  subroutine MOSE_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, Blk_MOSE, Is, Js, Ks, Fs, Blk_FUSS, Mg, Pg, new_flux, Ovar )
    use hydra_Parameters_m, only: guide_
    use MOSE_Advanced_Types_m
    use MOSE_Base_Types_m
    use MOSE_Global_m
    use MOSE_Lib_Fluid
    use MOSE_Lib_RANS
    use MOSE_Lib_BC_Fluxes, only: Compute_Wall_Properties
    use FLINT_Lib_Thermodynamic
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Solid
    implicit none
    integer,intent(in) :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(MOSE_block_type), intent(inout)   :: Blk_MOSE
    type(FUSS_block_type), intent(inout)   :: Blk_FUSS
    type(MOSE_tensor_3D_type), intent(in)  :: Mg
    real(R8), intent(in)                   :: Pg(1,6)
    real(R8), intent(inout)                :: new_flux(nprim)
    real(R8), optional, dimension(8), intent(out) :: Ovar
    ! Local
    integer  :: Dir, Face_im, Face_jm, Face_km, modfm, modfm1, modfm2, modfm3, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, M_loc(3,3), M_con(3,3)
    real(R8) :: T_loc, T_ip_loc, T_im_loc, T_jp_loc, T_jm_loc, T_kp_loc, T_km_loc, kappa_loc
    real(R8) :: T_con, kappa_con
    real(R8) :: Num_loc, Num_con, Den_loc, Den_con
    real(R8) :: rho, Rgas, mil, kl, Prim(nprim), Prim_wall(nprim), dist, dl, T
    real(R8) :: Tint, T_wall, q_wall, rho_wall, Stress(3), Gradient(nprim,3), Flux(nprim)

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )

    Ig = Im - guide_(Fm,1)
    Jg = Jm - guide_(Fm,2)
    Kg = Km - guide_(Fm,3)

    ! Metric stuff
    Normal = Blk_MOSE % dir(Dir) % f(Face_im,Face_jm,Face_km) % N
    Area   = Blk_MOSE % dir(Dir) % f(Face_im,Face_jm,Face_km) % A
    M_loc  = Blk_MOSE % M (Im,Jm,Km) % c 
    M_con  = Mg % c
    
    ! Stencil building               
    call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Im,Jm,Km), rho, Rgas )
    T_loc = EOS( p=Blk_MOSE%P(np,Im,Jm,Km), rho=rho, R=Rgas )
    call co_k_mi_lam_Wilke ( Blk_MOSE%P(1:nsc,Im,Jm,Km), rho, T_loc, mil, kappa_loc )
    call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Im+1,Jm,Km), rho, Rgas )
    T_ip_loc = EOS( p=Blk_MOSE%P(np,Im+1,Jm,Km), rho=rho, R=Rgas )
    call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Im-1,Jm,Km), rho, Rgas )
    T_im_loc = EOS( p=Blk_MOSE%P(np,Im-1,Jm,Km), rho=rho, R=Rgas )
    call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Im,Jm+1,Km), rho, Rgas )
    T_jp_loc = EOS( p=Blk_MOSE%P(np,Im,Jm+1,Km), rho=rho, R=Rgas )
    call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Im,Jm-1,Km), rho, Rgas )
    T_jm_loc = EOS( p=Blk_MOSE%P(np,Im,Jm-1,Km), rho=rho, R=Rgas )
    if (ndir == 2) then
      T_kp_loc = huge(0.0)
      T_km_loc = huge(0.0)
    else
      call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Im,Jm,Km+1), rho, Rgas )
      T_kp_loc = EOS( p=Blk_MOSE%P(np,Im,Jm,Km+1), rho=rho, R=Rgas )
      call co_rotot_Rtot ( Blk_MOSE%P(1:nsc,Im,Jm,Km-1), rho, Rgas )
      T_km_loc = EOS( p=Blk_MOSE%P(np,Im,Jm,Km-1), rho=rho, R=Rgas )
    endif

    T_con    = Pg(1,1)
    call co_Kappa( Blk_FUSS % matID(Is,Js,Ks), T_con, kappa_con )

    ! Tint computations
    select case ( Fm )
      case ( 1 : 2 )
        Den_loc = 2.0d0 * dot_product( M_loc(1,:), Normal ) * modfm2
        Num_loc = 2.0d0*T_loc * dot_product( M_loc(1,:), Normal ) * modfm2 &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(1,:), Normal ) * modfm2
        Num_con = 2.0d0*T_con * dot_product( M_con(1,:), Normal ) * modfm2 &
                + 0.5d0*(Pg(1,4)-Pg(1,3)) * dot_product( M_con(2,:), Normal ) &
                + 0.5d0*(Pg(1,6)-Pg(1,5)) * dot_product( M_con(3,:), Normal )
      case ( 3 : 4 )
        Den_loc = 2.0d0 * dot_product( M_loc(2,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(2,:), Normal ) * modfm2 &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(2,:), Normal ) * modfm2
        Num_con = 0.5d0*(Pg(1,4)-Pg(1,3)) * dot_product( M_con(1,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(2,:), Normal ) * modfm2 &
                + 0.5d0*(Pg(1,6)-Pg(1,5)) * dot_product( M_con(3,:), Normal )
      case ( 5 : 6 )
        Den_loc = 2.0d0 * dot_product( M_loc(3,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(3,:), Normal ) * modfm2
        Den_con = 2.0d0 * dot_product( M_con(3,:), Normal ) * modfm2
        Num_con = 0.5d0*(Pg(1,4)-Pg(1,3)) * dot_product( M_con(1,:), Normal ) &
                + 0.5d0*(Pg(1,6)-Pg(1,5)) * dot_product( M_con(2,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(3,:), Normal ) * modfm2
    end select
    
    Tint = (kappa_loc*Num_loc + kappa_con*Num_con) / (kappa_loc*Den_loc + kappa_con*Den_con)

    ! Calcolo dei flussi come nel caso di Temperature imposta Tint = Twall
    T_wall = Tint

    Prim = Blk_MOSE % P (:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    ! Boundary layer: p(wall) = p(cell). Also, non catalytic wall, so R(wall) = R(cell) and ci(wall) = ci(cell), but not roi
    rho_wall = Prim(np) / ( Rgas * T_wall )
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    ! Laminar viscosity/conductivity at interface computed at Twall
    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )

    ! Difference in Across-face direction: symmetry; tangential direction: null
    Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3                  ! Velocity
    Gradient(np,Dir) = ( T - T_wall ) * modfm3                  ! Temperature

    if (model==2)  then
      dl = Blk_MOSE % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil, Prim_Wall(nt:nprim), dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3 
    end if

    ! Transformation
    if (ndir == 2) Gradient(:,3) = 0.0
    Gradient = matmul ( Gradient, M_loc )

    ! Stress vector
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:nprim) )

    ! Fluxes
    Flux(1:nsc) = 0d0 
    Flux(nu:nw) = Stress * Area
    Flux(np) = 0.5 * Area * kl * dot_product ( Gradient(np,:), Normal )

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dist )
      ! Ghost-cell extrapolation of RANS variables
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                    rho, rho_wall, Blk_MOSE % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    new_flux = - 2.0 * Flux * modfm2

    ! Wall variables
    q_wall = kl * dot_product( Gradient(np,:), Normal )
    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), Tw=T_wall,          &
                                                    rhow=rho_wall, mu=mil,                          &
                                                    y=Blk_MOSE%dl(im,jm,km)%c(dir)*0.5d0,            &
                                                    qw=q_wall, &
                                                    exit_array=Ovar)

  end subroutine MOSE_BC_Coupled_Wall_Temperature


  subroutine Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    implicit none
    integer, intent(in)  :: Fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1-modFm  
    modfm2 = 1-2*modFm
    modfm3 = 2*( 2*modFm-1 )

  end subroutine Compute_Modfm

  
  subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
    implicit none
    integer, intent(in)  :: face, cell_i, cell_j, cell_k
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

end module hydraMF_Lib_BC_Coupled_Wall_Temperature