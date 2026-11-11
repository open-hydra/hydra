module hydraAF_Lib_BC_Coupled_Wall_Temperature
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: FUSS_BC_Coupled_Wall_Temperature, ARES_BC_Coupled_Wall_Temperature
  
contains

  !>@brief Solid-side coupled wall.
  !>
  !> The solid only ever conducts: its conductivity is kappa_loc, from co_Kappa on the
  !> FUSS cell material, and it is what the heat flux put into the FUSS residual uses.
  !> The routine also evaluates the FLUID conductivity kappa_con (FLINT tables, from the
  !> ARES cell state) because the interface temperature Tint is a weighted average of the
  !> two cell temperatures whose weights are the two conductivities: FUSS must reproduce
  !> exactly the number ARES computes, or the wall temperature becomes discontinuous.
  !>
  !> @param[in] k_fluid_wall  TOTAL conductivity OF THE FLUID on the wall face,
  !>                     k_l(p,Twall) + k_t, evaluated by the ARES-side routine and
  !>                     mirrored here by hydraAF_Mod_CouplingExchange. It replaces the
  !>                     local estimate kappa_con on the fluid side of Tint, and never
  !>                     touches kappa_loc (solid), which conducts and nothing else.
  !>                     It must be the very number the fluid heat flux is built with:
  !>                     Tint is the root of k_s*grad(T_s).n = k_f*grad(T_f).n, so any
  !>                     difference between the k_f that enters Tint and the one that
  !>                     multiplies the fluid gradient shows up directly as a mismatch
  !>                     between the two one-sided wall fluxes. Leaving k_t out costs a
  !>                     factor (1+k_t/k_l); evaluating k_l at the cell state instead of
  !>                     at the wall costs the variation of k_l across the first cell.
  !>                     K_WALL_UNSET (first call, before the fluid side has produced
  !>                     anything) falls back to kappa_con from the local fluid state.
  subroutine FUSS_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, Blk_FUSS, Is, Js, Ks, Fs, P_con, Mg, Tg, new_flux, Ovar, k_fluid_wall )
    use hydra_Parameters_m, only: guide_
    use FUSS_Advanced_Types_m
    use FUSS_Base_Types_m
    use FUSS_Lib_Solid
    use FUSS_Lib_BC_Fluxes, only: Compute_Wall_Properties
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    integer,intent(in)                     :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(FUSS_block_type), intent(inout)   :: Blk_FUSS
    real(R8), intent(in)                   :: P_con(:)   ! P(1:np) of ARES secondary cell (Is,Js,Ks)
    type(FUSS_tensor_3D_type), intent(in)  :: Mg
    real(R8), intent(in)                   :: Tg(6)
    real(R8), intent(inout)                :: new_flux
    real(R8), optional, dimension(2), intent(out) :: Ovar
    real(R8), optional, intent(in)         :: k_fluid_wall

    ! Local
    integer  :: Dir, Face_im, Face_jm, Face_km, modfm, modfm1, modfm2, modfm3, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, M_loc(3,3), M_con(3,3)
    real(R8) :: T_loc, T_ip_loc, T_im_loc, T_jp_loc, T_jm_loc, T_kp_loc, T_km_loc, kappa_loc
    real(R8) :: T_con, kappa_con, h_con, kappa_fluid
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
    h_con = pT2h( p=P_con(np), T=T_con )
    kappa_con = ph2vars ( P_con(np), h_con, k_tab2D )

    ! Conductivity to put on the FLUID side of Tint: the wall value the fluid solver
    ! actually uses in its heat flux, mirrored here (see the routine header). The local
    ! kappa_con above, taken at the fluid CELL state and with no turbulent part, is only
    ! the fallback for the very first call.
    kappa_fluid = kappa_con
    if ( present(k_fluid_wall) ) then
      if ( k_fluid_wall >= 0d0 ) kappa_fluid = k_fluid_wall
    endif

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

    Tint = ( kappa_loc*Num_loc + kappa_fluid*Num_con ) &
         / ( kappa_loc*Den_loc + kappa_fluid*Den_con )

    ! Calcolo dei flussi come nel caso di Temperature imposta Tint = Twall
    T_wall = Tint

    ! boundary cell variables
    T = Blk_FUSS % T (Im,Jm,Km)
    call co_Kappa( Blk_FUSS%matID(Im,Jm,Km), T, kappa )

    ! Wall-normal derivative from the half-cell one-sided difference, tangential ones
    ! from the cell-centred central differences. Both are needed: Tint above is the root
    ! of the FULL flux balance k_s*grad(T_s).n = k_f*grad(T_f).n, so a flux built from the
    ! across-face term alone does not reproduce the balance Tint was built from and the
    ! two one-sided wall fluxes disagree by the tangential part wherever the grid is not
    ! orthogonal to the wall. Signs follow the interior scheme (FUSS_Lib_Diffusive).
    Gradient = 0d0
    Gradient(Dir) = ( T - T_wall ) * modfm3
    select case ( Fm )
      case ( 1 : 2 )
        Gradient(2) = 0.5d0 * ( T_jp_loc - T_jm_loc )
        Gradient(3) = 0.5d0 * ( T_kp_loc - T_km_loc )
      case ( 3 : 4 )
        Gradient(1) = 0.5d0 * ( T_ip_loc - T_im_loc )
        Gradient(3) = 0.5d0 * ( T_kp_loc - T_km_loc )
      case ( 5 : 6 )
        Gradient(1) = 0.5d0 * ( T_ip_loc - T_im_loc )
        Gradient(2) = 0.5d0 * ( T_jp_loc - T_jm_loc )
    end select
    if (ndir == 2) Gradient(3) = 0.0
    Gradient = matmul ( Gradient, M_loc )

    Flux = Area * kappa * dot_product ( Gradient, Normal )
    
    ! Residual update
    new_flux = - Flux * modfm2
    
    ! Wall variables
    q_wall = kappa * dot_product( Gradient, Normal )
    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine FUSS_BC_Coupled_Wall_Temperature


  !>@brief Fluid-side coupled wall.
  !>
  !> Here kappa_loc is the fluid conductivity and kappa_con the solid one (co_Kappa on the
  !> partner FUSS cell material), i.e. the mirror image of the solid-side routine.
  !>
  !> @param[in]  k_fluid_wall      total fluid conductivity k_l+k_t on the wall face, used
  !>                          BOTH to build the interface temperature AND to build the
  !>                          heat flux. It is the value produced by the PREVIOUS call
  !>                          (k_l and k_t both depend on Twall, which depends on them):
  !>                          in a partitioned coupling the solid side can never have
  !>                          anything fresher, and using the same lagged number on both
  !>                          sides and in the flux is what makes the two one-sided wall
  !>                          fluxes identical at every iteration, not just at the fixed
  !>                          point. K_WALL_UNSET falls back to the local estimate.
  !> @param[out] k_fluid_wall_new  k_l(p,Twall)+k_t just evaluated on this face, stored for
  !>                          the next call and mirrored to the FUSS partner entry.
  subroutine ARES_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, Blk_ARES, Is, Js, Ks, Fs, matID_con, Mg, Pg, new_flux, Ovar, k_fluid_wall, k_fluid_wall_new )
    use hydra_Parameters_m, only: guide_
    use ARES_Advanced_Types_m
    use ARES_Base_Types_m
    use ARES_Global_m
    use ARES_Lib_Fluid
    use ARES_Lib_RANS
    use ARES_Config_Types_m, only: obj_rans
    use ARES_Lib_BC_Fluxes, only: Compute_Wall_Properties
    use ARES_Lib_Prt_Correction, only: delta_Prt
    use FLINT_Lib_Thermodynamic
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Solid
    implicit none
    integer,intent(in) :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(ARES_block_type), intent(inout)   :: Blk_ARES
    real(R8), intent(in)                   :: matID_con   ! FUSS secondary cell (Is,Js,Ks) matID
    type(ARES_tensor_3D_type), intent(in)  :: Mg
    real(R8), intent(in)                   :: Pg(1,6)
    real(R8), intent(inout)                :: new_flux(nprim)
    real(R8), optional, dimension(8), intent(out) :: Ovar
    real(R8), optional, intent(in)         :: k_fluid_wall
    real(R8), optional, intent(out)        :: k_fluid_wall_new
    ! Local
    integer  :: Dir, Face_im, Face_jm, Face_km, modfm, modfm1, modfm2, modfm3, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, M_loc(3,3), M_con(3,3)
    real(R8) :: T_loc, T_ip_loc, T_im_loc, T_jp_loc, T_jm_loc, T_kp_loc, T_km_loc, kappa_loc
    real(R8) :: T_con, kappa_con, kappa_fluid
    real(R8) :: Num_loc, Num_con, Den_loc, Den_con
    real(R8) :: rho, mil, kl, Prim(nprim), Prim_wall(nprim), dl, T
    real(R8) :: Tint, T_wall, q_wall, rho_wall, Stress(3), Gradient(nprim,3), Flux(nprim)
    real(R8) :: pressure, h_wall, mit, kt, cp, hs
    real(R8) :: dPrt ! Turbulent Prandtl correction for rough wall

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )

    Ig = Im - guide_(Fm,1)
    Jg = Jm - guide_(Fm,2)
    Kg = Km - guide_(Fm,3)

    ! Metric stuff
    Normal = Blk_ARES % dir(Dir) % f(Face_im,Face_jm,Face_km) % N
    Area   = Blk_ARES % dir(Dir) % f(Face_im,Face_jm,Face_km) % A
    M_loc  = Blk_ARES % M (Im,Jm,Km) % c 
    M_con  = Mg % c
    
    ! Roughness
    if ( obj_rans%rough ) then 
      hs = Blk_ARES%ks(Im,Jm,Km)
    else  
      hs = 0d0
    end if

    ! Stencil building               
    T_loc = ph2vars( Blk_ARES%P(np,Im,Jm,Km), Blk_ARES%P(nh,Im,Jm,Km), T_tab )
    kappa_loc = ph2vars( Blk_ARES%P(np,Im,Jm,Km), Blk_ARES%P(nh,Im,Jm,Km), k_tab2D )
    T_ip_loc = ph2vars( Blk_ARES%P(np,Im+1,Jm,Km), Blk_ARES%P(nh,Im+1,Jm,Km), T_tab )
    T_im_loc = ph2vars( Blk_ARES%P(np,Im-1,Jm,Km), Blk_ARES%P(nh,Im-1,Jm,Km), T_tab )
    T_jp_loc = ph2vars( Blk_ARES%P(np,Im,Jm+1,Km), Blk_ARES%P(nh,Im,Jm+1,Km), T_tab )
    T_jm_loc = ph2vars( Blk_ARES%P(np,Im,Jm-1,Km), Blk_ARES%P(nh,Im,Jm-1,Km), T_tab )
    if (ndir == 2) then
      T_kp_loc = huge(0.0)
      T_km_loc = huge(0.0)
    else
      T_kp_loc = ph2vars( Blk_ARES%P(np,Im,Jm,Km+1), Blk_ARES%P(nh,Im,Jm,Km+1), T_tab )
      T_km_loc = ph2vars( Blk_ARES%P(np,Im,Jm,Km-1), Blk_ARES%P(nh,Im,Jm,Km-1), T_tab )
    endif

    T_con    = Pg(1,1)
    call co_Kappa( matID_con, T_con, kappa_con )

    ! Total fluid conductivity on the wall face, lagged by one call and shared with the
    ! solid side (see the routine header). kappa_loc, taken at the CELL state and without
    ! the turbulent part, is only the fallback for the first call.
    kappa_fluid = kappa_loc
    if ( present(k_fluid_wall) ) then
      if ( k_fluid_wall >= 0d0 ) kappa_fluid = k_fluid_wall
    endif

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
    
    Tint = ( kappa_fluid*Num_loc + kappa_con*Num_con ) &
         / ( kappa_fluid*Den_loc + kappa_con*Den_con )

    ! Calcolo dei flussi come nel caso di Temperature imposta Tint = Twall
    T_wall = Tint

    Prim = Blk_ARES % P (:,Im,Jm,Km)
    rho  = ph2vars( Prim(np), Prim(nh), rho_tab )  ! cell density
    T    = ph2vars( Prim(np), Prim(nh), T_tab   ) ! cell Temperature

    ! Initialization
    Gradient = 0d0

    ! Boundary layer: p(wall) = p(cell). Also, non catalytic wall, so R(wall) = R(cell) and ci(wall) = ci(cell), but not roi
    pressure = Prim(np) ! wall pressure = cell pressure
    h_wall   = pT2h(pressure, T_wall ) ! Wall entalpy
    rho_wall = ph2vars(pressure, h_wall, rho_tab )   ! Wall density

    Prim_Wall = Prim        ! -------------------------------- ! 
    Prim_Wall(nh) = h_wall ! -------------------------------- !     

    ! Laminar viscosity/conductivity at interface computed at Twall
    mil = ph2vars( pressure, h_wall, mi_tab2D )
    kl  = ph2vars( pressure, h_wall, k_tab2D  )

    ! Across-face direction: symmetry. Tangential direction: null for velocity and for
    ! the RANS variables, the temperature one is filled just below.
    Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3                  ! Velocity
    Gradient(nh,Dir) = ( T - T_wall ) * modfm3                  ! Temperature
    ! Tangential part of the temperature gradient, the counterpart of the terms Tint is
    ! built with (see the solid-side routine): without it the heat flux is only the
    ! across-face projection and the two sides of the interface exchange different powers
    ! on a grid that is not orthogonal to the wall.
    select case ( Fm )
      case ( 1 : 2 )
        Gradient(nh,2) = 0.5d0 * ( T_jp_loc - T_jm_loc )
        Gradient(nh,3) = 0.5d0 * ( T_kp_loc - T_km_loc )
      case ( 3 : 4 )
        Gradient(nh,1) = 0.5d0 * ( T_ip_loc - T_im_loc )
        Gradient(nh,3) = 0.5d0 * ( T_kp_loc - T_km_loc )
      case ( 5 : 6 )
        Gradient(nh,1) = 0.5d0 * ( T_ip_loc - T_im_loc )
        Gradient(nh,2) = 0.5d0 * ( T_jp_loc - T_jm_loc )
    end select

    if (model==2)  then
      dl = Blk_ARES % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil_wall=mil, rans_wall=Prim_Wall(nt:nprim), rans_cell=Prim(nt:nprim), &
                                 metric=dot_product( M_loc(Dir,:), normal ), dist=dl, ks=hs )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3 
      call Eddy_Viscosity ( mut=mit, rans_variables=Prim_Wall(nt:nprim), &
                            mul=mil, rho=rho_wall, vel_gradient=Gradient(nu:nw,:), &
                            walldist=0d0, ks=hs)   ! wall face: y = 0 (dl stays the cell distance for Set_Wall_Values)
    else 
      mit = 0d0
    end if

    ! Transformation
    if (ndir == 2) Gradient(:,3) = 0.0
    Gradient = matmul ( Gradient, M_loc )

    ! Stress vector
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, mit, Prim(nt:) )
    cp = ph2vars(pressure, h_wall, hT_tab )
    if ( obj_rans%Prt_correction ) then
      ! Evaluated at the wall face (y=0). The rough-wall value mitilde_w = vk*u_tau*rho*0.03*hs
      ! yields the same u_tau as the cell value taken at y=dl, but leaves the near-wall damping
      ! exp(-y/hs) equal to 1, as required for a flux computed on the wall face itself.
      dPrt = delta_Prt ( mil, cp, kl, Prim_Wall(nt), hs, 0d0 )
    else
      dPrt = 0.00d0
    end if
    kt = mit*cp/(obj_rans%Prt+dPrt)

    ! Fluxes. The heat flux uses kappa_fluid, i.e. exactly the conductivity Tint was
    ! built with on both sides: (kl+kt) freshly evaluated here would be a different
    ! number and the interface would exchange two different powers. The fresh value is
    ! published through k_fluid_wall_new and becomes kappa_fluid at the next call.
    Flux(np) = 0d0 
    Flux(nu:nw) = Stress * Area
    Flux(nh) = 0.5d0 * Area * kappa_fluid * dot_product ( Gradient(nh,:), Normal )

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dl )
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk_ARES % P(nt:nprim,Ig,Jg,Kg) )
    endif


    ! Residual update
    new_flux = - 2.0 * Flux * modfm2

    ! Wall variables
    q_wall = kappa_fluid * dot_product( Gradient(nh,:), Normal )
    if (present(k_fluid_wall_new)) k_fluid_wall_new = kl + kt
    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), into=-modfm2, Tw=T_wall,          &
                                                    rhow=rho_wall, mu=mil,                          &
                                                    y=Blk_ARES%dl(im,jm,km)%c(dir)*0.5d0,                &
                                                    qw=q_wall,                                      &
                                                    exit_array=Ovar)

  end subroutine ARES_BC_Coupled_Wall_Temperature


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

end module hydraAF_Lib_BC_Coupled_Wall_Temperature