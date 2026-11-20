module hydraMF_Lib_BC_Coupled_Wall_Temperature
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: FUSS_BC_Coupled_Wall_Temperature, MOSE_BC_Coupled_Wall_Temperature

contains

  !>@brief Solid-side coupled wall.
  !>
  !> The solid only ever conducts: its conductivity is kappa_loc, from co_Kappa on the
  !> FUSS cell material, and it is what the heat flux put into the FUSS residual uses.
  !> The routine also evaluates the FLUID conductivity kappa_con (FLINT ideal-gas mixture,
  !> from the MOSE cell state) because the interface temperature Tint is a weighted average
  !> of the two sides whose weights are the two conductivities: FUSS must reproduce exactly
  !> the number MOSE computes, or the wall temperature becomes discontinuous.
  !>
  !> Each side enters Tint through its half of the balance, (A, Teff) of hydraMF_Lib_Metrics:
  !> the own one is computed here, the partner one by the rank that owns the partner cell
  !> (103: hydraMF_Lib_Ghost / hydraMF_Mod_CouplingExchange; 104: the segments), so that the
  !> two sides build the same Tint whatever the orientation of the two meshes.
  !>
  !> @param[in] P_con     P(1:nsc) of the MOSE partner cell: its species densities
  !> @param[in] T_con     temperature of the MOSE partner cell (fallback conductivity only)
  !> @param[in] A_con     1/d of the MOSE partner cell
  !> @param[in] Teff_con  temperature of the MOSE partner cell corrected by its tangential
  !>                      terms (and, on a 104 segment, moved to the segment centroid)
  !> @param[in] k_fluid_wall  TOTAL conductivity OF THE FLUID on the wall face,
  !>                     k_l(Twall) + k_t, evaluated by the MOSE-side routine and
  !>                     mirrored here by hydraMF_Mod_CouplingExchange. It replaces the
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
  !> @param[in] dx_seg    104 only: centroid of the segment minus the centre of this face;
  !>                      Teff of this side is moved there along the wall.
  subroutine FUSS_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, Blk_FUSS, P_con, T_con, A_con, Teff_con, &
                                               new_flux, Ovar, k_fluid_wall, dx_seg )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Solid
    use FUSS_Lib_BC_Fluxes, only: Compute_Wall_Properties
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, co_k_mi_lam_Wilke
    use hydraMF_Lib_Metrics, only: FUSS_Wall_Half, Compute_Modfm, Face_Index
    implicit none
    integer,intent(in)                     :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout)   :: Blk_FUSS
    real(R8), intent(in)                   :: P_con(:)   ! P(1:nsc) of the MOSE partner cell
    real(R8), intent(in)                   :: T_con, A_con, Teff_con
    real(R8), intent(inout)                :: new_flux
    real(R8), optional, dimension(2), intent(out) :: Ovar
    real(R8), optional, intent(in)         :: k_fluid_wall
    real(R8), optional, intent(in)         :: dx_seg(3)

    ! Local
    integer  :: Dir, Face_im, Face_jm, Face_km, modfm, modfm1, modfm2, modfm3
    real(R8) :: Area, kappa_loc, A_loc, Teff_loc, gpar(3)
    real(R8) :: kappa_con, kappa_fluid, rho, Rgas, mil
    real(R8) :: T_wall, q_wall, Flux

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )

    Area = Blk_FUSS % dir(Dir) % f(Face_im,Face_jm,Face_km) % A

    ! Own half of the balance, with the conductivity of the cell; on a 104 segment (only there)
    ! Teff is moved to the centroid of the segment with the gradient along the wall
    call co_Kappa( Blk_FUSS % matID(Im,Jm,Km), Blk_FUSS % T(Im,Jm,Km), kappa_loc )
    if ( present(dx_seg) ) then
      call FUSS_Wall_Half ( Blk_FUSS, Im, Jm, Km, Fm, A_loc, Teff_loc, gpar )
      Teff_loc = Teff_loc + dot_product ( gpar, dx_seg )
    else
      call FUSS_Wall_Half ( Blk_FUSS, Im, Jm, Km, Fm, A_loc, Teff_loc )
    end if

    ! Laminar conductivity of the gas of the partner cell at its temperature (fallback)
    call co_rotot_Rtot ( P_con(1:nsc), rho, Rgas )
    call co_k_mi_lam_Wilke ( P_con(1:nsc), rho, T_con, mil, kappa_con )

    ! Conductivity to put on the FLUID side of Tint: the wall value the fluid solver
    ! actually uses in its heat flux, mirrored here (see the routine header). The local
    ! kappa_con above, taken at the fluid CELL state and with no turbulent part, is only
    ! the fallback for the very first call.
    kappa_fluid = kappa_con
    if ( present(k_fluid_wall) ) then
      if ( k_fluid_wall >= 0d0 ) kappa_fluid = k_fluid_wall
    endif

    ! Interface temperature: root of the balance of the two halves
    T_wall = ( kappa_loc*A_loc*Teff_loc + kappa_fluid*A_con*Teff_con ) &
           / ( kappa_loc*A_loc          + kappa_fluid*A_con          )

    ! Wall heat flux of this side, k grad(T).N = k modfm2 A (Tw - Teff): the same metric,
    ! tangential terms included, Tint was built with, so the two one-sided wall fluxes agree
    q_wall = kappa_loc * modfm2 * A_loc * ( T_wall - Teff_loc )
    Flux   = Area * q_wall

    ! Residual update
    new_flux = - Flux * modfm2

    ! Wall variables
    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine FUSS_BC_Coupled_Wall_Temperature


  !>@brief Fluid-side coupled wall.
  !>
  !> Here kappa_loc is the fluid conductivity and kappa_con the solid one (co_Kappa on the
  !> partner FUSS cell material, at the partner cell temperature T_con, as the FUSS side
  !> takes it), i.e. the mirror image of the solid-side routine. The gas at the wall, the
  !> wall values of the turbulence variables, the eddy viscosity, the wall stress and the
  !> turbulence flux are those of MOSE's own wall of prescribed temperature
  !> (BC_Wall_Temperature in MOSE_Lib_BC_Fluxes_Wall_Temperature), with T_wall the root of
  !> the balance of the two halves.
  !>
  !> @param[in]  matID_con, T_con  material and temperature of the FUSS partner cell
  !> @param[in]  A_con, Teff_con   half of the balance of the FUSS partner cell (see the
  !>                          solid-side routine)
  !> @param[in]  k_fluid_wall      total fluid conductivity k_l+k_t on the wall face, used
  !>                          BOTH to build the interface temperature AND to build the
  !>                          heat flux. It is the value produced by the PREVIOUS call
  !>                          (k_l and k_t both depend on Twall, which depends on them):
  !>                          in a partitioned coupling the solid side can never have
  !>                          anything fresher, and using the same lagged number on both
  !>                          sides and in the flux is what makes the two one-sided wall
  !>                          fluxes identical at every iteration, not just at the fixed
  !>                          point. K_WALL_UNSET falls back to the local estimate.
  !> @param[out] k_fluid_wall_new  k_l(Twall)+k_t just evaluated on this face, stored for
  !>                          the next call and mirrored to the FUSS partner entry.
  !> @param[in]  dx_seg       104 only: centroid of the segment minus the centre of this face.
  !> @param[in]  k_rough      sand-grain roughness of the face (MOSE bc record; absent = smooth)
  subroutine MOSE_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, Blk_MOSE, matID_con, T_con, A_con, Teff_con, &
                                               new_flux, Ovar, k_fluid_wall, k_fluid_wall_new, dx_seg, k_rough )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Lib_Fluid
    use MOSE_Lib_RANS
    use MOSE_Config_Types_m, only: obj_rans
    use MOSE_Lib_BC_Fluxes, only: Compute_Wall_Properties, Wall_State_Temperature
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, co_k_mi_lam_Wilke, f_cp
    use FUSS_Lib_Solid
    use hydraMF_Lib_Metrics, only: MOSE_Wall_Half, Compute_Modfm, Face_Index
    implicit none
    integer,intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout)   :: Blk_MOSE
    real(R8), intent(in)                   :: matID_con   ! FUSS partner cell matID
    real(R8), intent(in)                   :: T_con, A_con, Teff_con
    real(R8), intent(inout)                :: new_flux(nprim)
    real(R8), optional, dimension(8), intent(out) :: Ovar
    real(R8), optional, intent(in)         :: k_fluid_wall
    real(R8), optional, intent(out)        :: k_fluid_wall_new
    real(R8), optional, intent(in)         :: dx_seg(3)
    real(R8), optional, intent(in)         :: k_rough
    ! Local
    integer  :: Dir, Face_im, Face_jm, Face_km, modfm, modfm1, modfm2, modfm3
    real(R8) :: Normal(3), Area, Dist, Mu(3,3)
    real(R8) :: kappa_loc, A_loc, Teff_loc, gpar(3)
    real(R8) :: kappa_con, kappa_fluid
    real(R8) :: rho, Rgas, T, mil_c, mil, kl, Prim(nprim), Prim_wall(nprim), dl
    real(R8) :: T_wall, q_wall, rho_wall, Stress(3), Gradient(nprim,3), Flux(nprim)
    real(R8) :: mit, kt, kr

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )

    kr = 0d0
    if ( present(k_rough) ) kr = k_rough

    ! Metric stuff: Mu for the gradients of the velocity and of the turbulence variables
    ! (hydraMF_Lib_Metrics)
    Normal = Blk_MOSE % dir(Dir) % f(Face_im,Face_jm,Face_km) % N
    Area   = Blk_MOSE % dir(Dir) % f(Face_im,Face_jm,Face_km) % A
    Dist   = 1d-20 ! the turbulence flux is computed at the wall, as in MOSE BC_Wall_Temperature

    ! Cell state, and the laminar conductivity of the cell (the fallback of the first call)
    Prim = Blk_MOSE % P (:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )
    call co_k_mi_lam_Wilke ( Prim(1:nsc), rho, T, mil_c, kappa_loc )

    ! Own half of the balance (on a 104 segment, only there, Teff moved to its centroid)
    if ( present(dx_seg) ) then
      call MOSE_Wall_Half ( Blk_MOSE, Im, Jm, Km, Fm, A_loc, Teff_loc, gpar, Mu )
      Teff_loc = Teff_loc + dot_product ( gpar, dx_seg )
    else
      call MOSE_Wall_Half ( Blk_MOSE, Im, Jm, Km, Fm, A_loc, Teff_loc, Mu=Mu )
    end if

    call co_Kappa( matID_con, T_con, kappa_con )

    ! Total fluid conductivity on the wall face, lagged by one call and shared with the
    ! solid side (see the routine header). kappa_loc, taken at the CELL state and without
    ! the turbulent part, is only the fallback for the first call.
    kappa_fluid = kappa_loc
    if ( present(k_fluid_wall) ) then
      if ( k_fluid_wall >= 0d0 ) kappa_fluid = k_fluid_wall
    endif

    ! Interface temperature: root of the balance of the two halves
    T_wall = ( kappa_fluid*A_loc*Teff_loc + kappa_con*A_con*Teff_con ) &
           / ( kappa_fluid*A_loc          + kappa_con*A_con          )

    ! Initialization
    Gradient = 0d0

    ! Boundary layer: p(wall) = p(cell). Also, non catalytic wall, so R(wall) = R(cell) and ci(wall) = ci(cell), but
    ! not roi. Laminar viscosity/conductivity at the interface computed at T_wall: the wall state of MOSE's own walls
    Prim_wall = Prim
    call Wall_State_Temperature ( Prim, rho, Rgas, T_wall, rho_wall, Prim_wall(1:nsc), mil, kl )

    ! Across-face direction: symmetry. Tangential direction: null for velocity and for the
    ! RANS variables (wall values). The temperature gradient is not needed here: its normal
    ! component, the only one the heat flux takes, is modfm2 A (Tw - Teff).
    Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3                  ! Velocity

    ! Wall values from the cell values per unit mass at the wall density, across the wall
    ! distance of the cell, as MOSE BC_Wall_Temperature: per unit mass the same values as the
    ! ghost cells of MOSE Ghost_Wall, which take the state of the cell
    if (model==2)  then
      dl = Blk_MOSE % yn(Im,Jm,Km)
      call RANS_Set_Wall_Values( mil, Prim(nt:nprim) * rho_wall / rho, Prim_Wall(nt:nprim), dl, kr )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
    end if

    ! Transformation
    if (ndir == 2) Gradient(:,3) = 0.0
    Gradient = matmul ( Gradient, Mu )

    ! Eddy viscosity at the wall face, with the velocity gradient in physical space (after the
    ! transformation), as MOSE BC_Wall_Temperature: a smooth wall zeroes the turbulence
    ! variables there, so it is zero; a rough wall leaves nit_wall > 0 (SA-rough)
    mit = 0d0
    kt  = 0d0
    if ( model==2 .and. kr > 0d0 ) then
      call Eddy_Viscosity ( mut=mit, rans_variables=Prim_Wall(nt:nprim), mul=mil, rho=rho_wall, &
                            vel_gradient=Gradient(nu:nw,:), walldist=0d0, k_rough=kr )
      kt = mit * f_cp ( Prim_wall(1:nsc), T_wall, rho_wall ) / obj_rans%Prt
    end if

    ! Stress vector, with the turbulence variables of the wall (no Reynolds stress, no k), as
    ! in MOSE BC_Wall_Temperature
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, mit, Prim_Wall(nt:) )

    ! Fluxes. The heat flux uses kappa_fluid, i.e. exactly the conductivity Tint was
    ! built with on both sides: (kl+kt) freshly evaluated here would be a different
    ! number and the interface would exchange two different powers. The fresh value is
    ! published through k_fluid_wall_new and becomes kappa_fluid at the next call.
    q_wall = kappa_fluid * modfm2 * A_loc * ( T_wall - Teff_loc )
    Flux        = 0d0
    Flux(nu:nw) = Stress * Area
    Flux(np)    = 0.5d0 * Area * q_wall

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=Dist )
      ! No ghost value is written here: the ghost cells of the coupled walls, turbulence
      ! variables included, are those of MOSE Fill_Ghost_Cell (Ghost_Wall), one set for the
      ! fluxes and the source terms
    endif

    ! Residual update
    new_flux = - 2.0 * Flux * modfm2

    ! Wall variables
    if (present(k_fluid_wall_new)) k_fluid_wall_new = kl + kt
    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=Prim(np), into=-modfm2, Tw=T_wall,          &
                                                    rhow=rho_wall, mu=mil,                          &
                                                    y=Blk_MOSE%dl(Im,Jm,Km)%c(Dir)*0.5d0,                &
                                                    qw=q_wall,                                      &
                                                    exit_array=Ovar)

  end subroutine MOSE_BC_Coupled_Wall_Temperature

end module hydraMF_Lib_BC_Coupled_Wall_Temperature
