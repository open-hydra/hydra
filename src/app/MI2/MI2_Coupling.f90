module hydraMI2_Mod_Coupling
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: MOSExIGLOO
  public :: MOSExICE
  public :: ICExMOSE
  public :: Gas_Field_For_IGLOO
  public :: IGLOOxMOSE
  public :: Initialize_Residual_Store
  public :: Build_Family_Material_Map
  public :: Accumulate_Injected_Mass
  public :: Gas_Species_Mass
  public :: Accumulate_Injected_Energy
  public :: Accumulate_Injected_ICE
  public :: Gas_Total_Energy
  public :: Gas_Mean_Temperature
  public :: Check_ICE_Mass_Channel
  public :: Check_ICE_Source_Finite
  public :: Check_IGLOO_Source_Finite
  public :: Check_IGLOO_Euler_Finite
  public :: Gas_Momentum
  public :: ICE_To_Gas_Totals
  public :: IGLOO_To_Gas_Totals

  !> Under-relaxation store: the previous coupling step's injected residual, per block.
  type, public :: resField
    real(R8), dimension(:,:,:,:), allocatable :: resBlock
  end type resField

  !> family -> parcel material id, built once after IGLOO setup (getSourceTerms needs it).
  integer, allocatable, private :: famMatID(:)

  !> Number of (cell, species) extractions clamped by the consumed-species floor guard.
  integer, public :: nMassLimited = 0

  !> Mass MI2 has injected into each gas species over the run [kg].
  real(R8), allocatable, public :: massInjected(:)

  !> Energy MI2 has handed the gas over the run [J] (sourceEn plus datum correction), and the
  !  datum-correction part of it.
  real(R8), public :: energyInjected  = 0._R8
  real(R8), public :: energyDatumPart = 0._R8
  logical,  public :: energyAudited   = .false.

  !> Inflow mass flow MOSExIGLOO handed IGLOO at its last call [kg/s], and the number of face records.
  real(R8), public :: mdotHanded = 0._R8
  integer,  public :: nHanded    = 0

contains

  !> Mass of each gas species in the domain [kg] = SUM_cells rho_s*vol.
  subroutine Gas_Species_Mass ( MOSE_domain, mass )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: mass(nsc)
    integer :: b, i, j, k, s

    mass = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        do s = 1, nsc
          mass(s) = mass(s) + MOSE_domain%blk(b)%P(s,i,j,k) * MOSE_domain%blk(b)%vol(i,j,k)
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine Gas_Species_Mass


  !> Gas momentum [kg m/s] = SUM_cells rho*u*vol.
  subroutine Gas_Momentum ( MOSE_domain, mom )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc, nu, nw
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: mom(3)
    integer :: b, i, j, k

    mom = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        mom = mom + sum(MOSE_domain%blk(b)%P(1:nsc,i,j,k)) * MOSE_domain%blk(b)%P(nu:nw,i,j,k) &
                  * MOSE_domain%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
    enddo

  end subroutine Gas_Momentum


  !> What ICE's sources hand the gas over the interior cells and the families, as ICExMOSE applies them: force
  !  F [N], energy rate Q [W] (with the datum offset of a material that converts), mass
  !  rate M [kg/s] into the gas species, and the range of the relaxation time tau [s].
  subroutine ICE_To_Gas_Totals ( MOSE_domain, ICE_domain, F, Q, M, taumin, taumax )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    type(ICE_domain_type),  intent(in)  :: ICE_domain
    real(R8),               intent(out) :: F(3), Q, M, taumin, taumax
    integer  :: b, g, i, j, k, imat
    real(R8) :: vol, dm

    F = 0._R8; Q = 0._R8; M = 0._R8; taumin = huge(1._R8); taumax = -huge(1._R8)
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        imat = mat_of(g)
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          associate ( cp => ICE_domain%blk(b)%cond_phase(g) )
            F = F - cp%source(2:4,i,j,k) * vol
            Q = Q - cp%source(nbase(g)-1,i,j,k) * vol
            taumin = min(taumin, cp%tau(i,j,k))
            taumax = max(taumax, cp%tau(i,j,k))
            if (has_stoichiometry_ICE) then
              dm = - cp%source(1,i,j,k) * vol
              M  = M + sum(nu_ice(imat,1:nsc)) * dm
              Q  = Q + dhDatumICE(imat) * dm
            endif
          end associate
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine ICE_To_Gas_Totals


  !> What IGLOO's source-mode rates hand the gas over the interior cells, before under-relaxation:
  !  force F [N], energy rate Q [W] with the datum correction, mass rate M [kg/s] into the gas species.
  subroutine IGLOO_To_Gas_Totals ( MOSE_domain, IGLOO, F, Q, M )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nm
    use hydraMI2_Stoichiometry, only: nu_gas => nu, has_stoichiometry, dhDatum, capped, nuCap, nuRel, dhProduct
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    type(obj_IGLOO),        intent(in)  :: IGLOO
    real(R8),               intent(out) :: F(3), Q, M
    integer  :: b, i, j, k, imat
    real(R8) :: cr(4)

    F = 0._R8; Q = 0._R8; M = 0._R8
    do b = 1, size(MOSE_domain%blk)
      if (.not. allocated(IGLOO%source(b)%sourceMom)) cycle
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        F = F + IGLOO%source(b)%sourceMom(:,i,j,k)
        Q = Q + IGLOO%source(b)%sourceEn(i,j,k)
        if (.not. has_stoichiometry) cycle
        do imat = 1, nm
          if (capped(imat)) then
            cr = Cap_Rates( imat, IGLOO%source(b)%sourceMass(imat,i,j,k), IGLOO%source(b)%sourceMass(nm+2*imat-1,i,j,k), &
                            IGLOO%source(b)%sourceMass(nm+2*imat,i,j,k) )
            M = M + sum(nu_gas(imat,1:nsc))*cr(2) + sum(nuCap(imat,1:nsc))*cr(3) + sum(nuRel(imat,1:nsc))*cr(4)
            Q = Q + dhDatum(imat)*cr(1) - dhProduct(imat)*(cr(3) - cr(4))
            cycle
          endif
          M = M + sum(nu_gas(imat,1:nsc)) * IGLOO%source(b)%sourceMass(imat,i,j,k)
          Q = Q + dhDatum(imat) * IGLOO%source(b)%sourceMass(imat,i,j,k)
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine IGLOO_To_Gas_Totals


  !> The rates of a capped IGLOO material [kg/s] from its net mass column and its cap's two: (1) the material
  !  consumed, (2) the part of it burnt into gas species, (3) the oxide its cap gained, (4) the cap oxide handed
  !  to the gas.
  pure function Cap_Rates ( m, net, gained, handed ) result(cr)
    use hydraMI2_Stoichiometry, only: nuOxide
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: net, gained, handed
    real(R8) :: cr(4)
    cr(3) = gained
    cr(4) = handed
    cr(1) = net + gained - handed
    cr(2) = cr(1) - gained/nuOxide(m)
  end function Cap_Rates


  !> Accumulate the mass handed to the gas over a coupling step of duration dtBlock:
  !  SUM_m nu(m,s)*sourceMass(m) (absolute kg/s, no volume factor), a capped material through its three rows.
  !  Called once per coupling step.
  subroutine Accumulate_Injected_Mass ( IGLOO, dtBlock )
    use MOSE_Global_m,          only: nsc
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nm
    use hydraMI2_Stoichiometry, only: nu_gas => nu, has_stoichiometry, capped, nuCap, nuRel
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    real(R8),        intent(in) :: dtBlock
    integer  :: b, m, s
    real(R8) :: tot, cr(4)

    if (.not. has_stoichiometry) return
    if (.not. allocated(massInjected)) then
      allocate(massInjected(1:nsc)); massInjected = 0._R8
    endif
    if (dtBlock <= 0._R8) return

    do b = 1, size(IGLOO%source)
      if (.not. allocated(IGLOO%source(b)%sourceMass)) cycle
      do m = 1, nm
        tot = sum(IGLOO%source(b)%sourceMass(m,:,:,:))
        if (capped(m)) then
          cr = Cap_Rates( m, tot, sum(IGLOO%source(b)%sourceMass(nm+2*m-1,:,:,:)), &
                          sum(IGLOO%source(b)%sourceMass(nm+2*m,:,:,:)) )
          do s = 1, nsc
            massInjected(s) = massInjected(s) + (nu_gas(m,s)*cr(2) + nuCap(m,s)*cr(3) + nuRel(m,s)*cr(4)) * dtBlock
          enddo
          cycle
        endif
        do s = 1, nsc
          massInjected(s) = massInjected(s) + nu_gas(m,s) * tot * dtBlock
        enddo
      enddo
    enddo

  end subroutine Accumulate_Injected_Mass

  !> Total gas energy [J] = SUM_cells cons(np)*vol (kinetic energy included).
  subroutine Gas_Total_Energy ( MOSE_domain, energy )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: np
    use FLINT_Lib_Thermodynamic, only: prim2cons
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: energy
    integer  :: b, i, j, k
    real(R8) :: cons(np)

    energy = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        cons   = prim2cons( MOSE_domain%blk(b)%P(1:np,i,j,k) )
        energy = energy + cons(np) * MOSE_domain%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
    enddo

  end subroutine Gas_Total_Energy


  !> Mass-weighted mean gas temperature [K], gas mass [kg] and mass-weighted frozen cv [J/kg/K].
  subroutine Gas_Mean_Temperature ( MOSE_domain, Tbar, mgas, cvbar )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: nsc, np
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_tabT, cp_tab, Ri_tab
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: Tbar, mgas, cvbar
    integer  :: b, i, j, k, s
    real(R8) :: rho, R, T, dm, cvcell

    Tbar = 0._R8; mgas = 0._R8; cvbar = 0._R8
    do b = 1, size(MOSE_domain%blk)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        call co_rotot_Rtot ( MOSE_domain%blk(b)%P(1:nsc,i,j,k), rho, R )
        T    = EOS ( MOSE_domain%blk(b)%P(np,i,j,k), rho, R )
        dm   = rho * MOSE_domain%blk(b)%vol(i,j,k)
        cvcell = 0._R8
        do s = 1, nsc
          cvcell = cvcell + MOSE_domain%blk(b)%P(s,i,j,k) * ( f_tabT(T,s,cp_tab) - Ri_tab(s) )
        enddo
        mgas  = mgas  + dm
        Tbar  = Tbar  + dm*T
        cvbar = cvbar + cvcell * MOSE_domain%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
    enddo
    if (mgas > 0._R8) then
      Tbar  = Tbar  / mgas
      cvbar = cvbar / mgas
    endif

  end subroutine Gas_Mean_Temperature


  !> Accumulate the energy handed to the gas over a coupling step of duration dtBlock:
  !  sourceEn + SUM_m dhDatum(m)*sourceMass(m) (absolute W), a capped material with its oxide's offset. Called
  !  once per coupling step.
  subroutine Accumulate_Injected_Energy ( IGLOO, dtBlock )
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nm
    use hydraMI2_Stoichiometry, only: has_stoichiometry, dhDatum, capped, dhProduct
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    real(R8),        intent(in) :: dtBlock
    integer  :: b, m
    real(R8) :: datum, cr(4)

    if (dtBlock <= 0._R8) return

    do b = 1, size(IGLOO%source)
      if (.not. allocated(IGLOO%source(b)%sourceEn)) cycle
      energyInjected = energyInjected + sum(IGLOO%source(b)%sourceEn) * dtBlock
      energyAudited  = .true.

      if (.not. has_stoichiometry) cycle
      if (.not. allocated(IGLOO%source(b)%sourceMass)) cycle
      datum = 0._R8
      do m = 1, nm
        if (capped(m)) then
          cr = Cap_Rates( m, sum(IGLOO%source(b)%sourceMass(m,:,:,:)), sum(IGLOO%source(b)%sourceMass(nm+2*m-1,:,:,:)), &
                          sum(IGLOO%source(b)%sourceMass(nm+2*m,:,:,:)) )
          datum = datum + dhDatum(m)*cr(1) - dhProduct(m)*(cr(3) - cr(4))
          cycle
        endif
        datum = datum + dhDatum(m) * sum(IGLOO%source(b)%sourceMass(m,:,:,:))
      enddo
      energyInjected  = energyInjected  + datum*dtBlock
      energyDatumPart = energyDatumPart + datum*dtBlock
    enddo

  end subroutine Accumulate_Injected_Energy


  !> Accumulate what ICE's sources hand the gas over a coupling step of duration dtBlock, as ICExMOSE applies them:
  !  the mass of each material that converts, through its row, and the energy with its datum offset. Called once
  !  per coupling step when an ICE material converts.
  subroutine Accumulate_Injected_ICE ( MOSE_domain, ICE_domain, dtBlock )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(ICE_domain_type),  intent(in) :: ICE_domain
    real(R8),               intent(in) :: dtBlock
    integer  :: b, g, i, j, k, s, imat
    real(R8) :: vol, mtot, etot

    if (.not. has_stoichiometry_ICE) return
    if (.not. allocated(massInjected)) then
      allocate(massInjected(1:nsc)); massInjected = 0._R8
    endif
    if (dtBlock <= 0._R8) return

    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%source)) cycle
        imat = mat_of(g)
        mtot = 0._R8; etot = 0._R8
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          associate ( cp => ICE_domain%blk(b)%cond_phase(g) )
            mtot = mtot - cp%source(1,i,j,k) * vol
            etot = etot - cp%source(nbase(g)-1,i,j,k) * vol
          end associate
        enddo; enddo; enddo
        do s = 1, nsc
          massInjected(s) = massInjected(s) + nu_ice(imat,s) * mtot * dtBlock
        enddo
        energyInjected  = energyInjected  + (etot + dhDatumICE(imat)*mtot) * dtBlock
        energyDatumPart = energyDatumPart + dhDatumICE(imat)*mtot * dtBlock
        energyAudited   = .true.
      enddo
    enddo

  end subroutine Accumulate_Injected_ICE


  !> Resolve family -> material once, after IGLOO%setup has assigned famID.
  subroutine Build_Family_Material_Map ( IGLOO )
    use IGLOO_module,    only: obj_IGLOO
    use IGLOO_variables, only: nfam
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    integer :: m, g, fam

    if (allocated(famMatID)) deallocate(famMatID)
    allocate(famMatID(1:max(nfam,1)))
    famMatID = 1

    do m = 1, size(IGLOO%material)
      do g = 1, IGLOO%material(m)%ngroups
        fam = IGLOO%material(m)%group(g)%famID
        if (fam >= 1 .and. fam <= nfam) famMatID(fam) = m
      enddo
    enddo

  end subroutine Build_Family_Material_Map


  !> Allocate the under-relaxation store, sized from the MOSE domain.
  subroutine Initialize_Residual_Store ( MOSE_domain, oldResidual )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nprim
    implicit none
    type(MOSE_domain_type), intent(in)               :: MOSE_domain
    type(resField), allocatable, intent(out)         :: oldResidual(:)
    integer :: b

    allocate(oldResidual(1:size(MOSE_domain%blk)))
    do b = 1, size(oldResidual)
      allocate(oldResidual(b)%resBlock(1:nprim, 1:MOSE_domain%blk(b)%dim(1), &
                                                1:MOSE_domain%blk(b)%dim(2), &
                                                1:MOSE_domain%blk(b)%dim(3)))
      oldResidual(b)%resBlock = 0._R8
    enddo

  end subroutine Initialize_Residual_Store


  !> Inject IGLOO's reaction on the gas into MOSE's residual (R = R - S), from MOSE's per-RK-stage
  !  callback: runs inside MOSE's parallel region, orphaned !$OMP DO only.
  !    source branch : source{Mass,Mom,En} are absolute rates (no *vol); mass through nu(m,s),
  !                    plus the enthalpy-datum correction and the consumed-species floor guard
  !    euler  branch : getSourceTerms is per unit volume (*vol); momentum and energy only
  subroutine IGLOOxMOSE ( MOSE_domain, IGLOO, eulCoupling, relaxation, iter, iMax, relCoeff, &
                          oldResidual )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: nsc, nu, nw, np, nprim
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    use IGLOO_module,            only: obj_IGLOO
    use IGLOO_variables,         only: nfam, nm
    use hydraMI2_Stoichiometry,  only: nu_gas => nu, has_stoichiometry, dhDatum, capped, nuCap, nuRel, dhProduct
    use MOSE_Config_Types_m,     only: obj_chemistry
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: eulCoupling, relaxation
    integer,                intent(in)    :: iter, iMax
    real(R8),               intent(in)    :: relCoeff
    type(resField),         intent(inout) :: oldResidual(:)
    ! Local
    real(R8), parameter :: toll = 1.d-15
    logical  :: useRelax
    integer  :: b, i, j, k, m
    ! mirrors blk%R(1:nprim): species 1:nsc, momentum nu:nw, energy np, then soot/passive/RANS (left at zero)
    real(R8) :: partial_residuals(nprim), sources(4)
    real(R8) :: rhoi(nsc), rho, u(3), T, p, R, gam, kappa, mu, coeff
    integer  :: imat, isp
    real(R8) :: dm, avail, dtStage, cr(4)

    ! relaxation ramp
    useRelax = (relCoeff > toll .and. relCoeff < (1.d0-toll)) .and. relaxation
    if (iter < iMax*(1.d0-relCoeff) .and. useRelax) then
      coeff = real(iter,R8) / real(iMax,R8)
      coeff = coeff*(1.d0/(1.d0-relCoeff) - exp(-(coeff/relCoeff))**2)
      coeff = max(0.d0, min(1.d0, coeff))
    elseif ((iMax==1) .and. useRelax) then
      coeff = relCoeff
    else
      coeff = 1.d0
      useRelax = .false.
    endif

    do b = 1, size(MOSE_domain%blk)

      !$OMP DO COLLAPSE(3) PRIVATE(partial_residuals,sources,rhoi,rho,u,p,T,R,gam,kappa,mu,i,j,k,m, &
      !$OMP&                       imat,isp,dm,avail,dtStage,cr)
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)

        partial_residuals(:) = 0._R8

        if (eulCoupling) then

          rhoi = MOSE_domain%blk(b)%P(1:nsc,i,j,k)
          u    = MOSE_domain%blk(b)%P(nu:nw,i,j,k)
          p    = MOSE_domain%blk(b)%P(np,i,j,k)
          call co_rotot_Rtot ( rhoi, rho, R )
          gam = f_gamma ( rhoi, p, rho, R )
          T   = EOS ( p, rho, R )
          call co_k_mi_lam_Wilke ( rhoi, rho, T, mu, kappa )

          sources = 0._R8
          do m = 1, nfam
            sources = sources + IGLOO%getSourceTerms( [rho,u,T,gam,R,mu,kappa],      &
                                    IGLOO%euler(b,m)%velocity(:,i,j,k),              &
                                    IGLOO%euler(b,m)%temperature(i,j,k),             &
                                    IGLOO%euler(b,m)%density(i,j,k),                 &
                                    IGLOO%euler(b,m)%np(i,j,k), famMatID(m) )
          enddo
          partial_residuals(nu:nw) = - sources(1:3)
          partial_residuals(np)    = - sources(4)
          partial_residuals = partial_residuals * MOSE_domain%blk(b)%vol(i,j,k)

        else

          partial_residuals(nu:nw) = IGLOO%source(b)%sourceMom(:,i,j,k)
          partial_residuals(np)    = IGLOO%source(b)%sourceEn(i,j,k)

          ! mass coupling: sourceMass per material, distributed over the gas species through nu
          if (has_stoichiometry) then
            do imat = 1, nm
              dm = IGLOO%source(b)%sourceMass(imat,i,j,k)
              ! an oxide cap: the material burnt into gas species, the oxide the cap gained, the cap handed to the gas
              if (capped(imat)) then
                cr = Cap_Rates( imat, dm, IGLOO%source(b)%sourceMass(nm+2*imat-1,i,j,k), &
                                IGLOO%source(b)%sourceMass(nm+2*imat,i,j,k) )
                do isp = 1, nsc
                  partial_residuals(isp) = partial_residuals(isp) + nu_gas(imat,isp)*cr(2) + nuCap(imat,isp)*cr(3) &
                                                                  + nuRel(imat,isp)*cr(4)
                enddo
                partial_residuals(np) = partial_residuals(np) + dhDatum(imat)*cr(1) - dhProduct(imat)*(cr(3) - cr(4))
                cycle
              endif
              if (dm == 0._R8) cycle
              do isp = 1, nsc
                if (nu_gas(imat,isp) == 0._R8) cycle
                partial_residuals(isp) = partial_residuals(isp) + nu_gas(imat,isp)*dm
              enddo

              ! enthalpy-datum correction, on the same dm as the mass term
              partial_residuals(np) = partial_residuals(np) + dhDatum(imat)*dm
            enddo

            ! floor guard: cap the extraction RATE of a consumed species at the rate that empties the
            ! cell over this stage's step, dtlocal*strangcoeff (the step MOSE's Newstate applies)
            dtStage = MOSE_domain%blk(b)%dtlocal(i,j,k) * obj_chemistry%strangcoeff
            if (dtStage > 0._R8) then
              do isp = 1, nsc
                if (partial_residuals(isp) >= 0._R8) cycle
                avail = MOSE_domain%blk(b)%P(isp,i,j,k) * MOSE_domain%blk(b)%vol(i,j,k) / dtStage
                if (-partial_residuals(isp) > avail) then
                  partial_residuals(isp) = -avail
                  !$OMP ATOMIC UPDATE
                  nMassLimited = nMassLimited + 1
                endif
              enddo
            endif
          endif

        endif

        if (useRelax) then
          partial_residuals(1:nsc) = coeff*partial_residuals(1:nsc) &
                                   + (1._R8-coeff)*oldResidual(b)%resBlock(1:nsc,i,j,k)
          partial_residuals(nu:nw) = coeff*partial_residuals(nu:nw) &
                                   + (1._R8-coeff)*oldResidual(b)%resBlock(nu:nw,i,j,k)
          partial_residuals(np)    = coeff*partial_residuals(np) &
                                   + (1._R8-coeff)*oldResidual(b)%resBlock(np,i,j,k)
          oldResidual(b)%resBlock(:,i,j,k) = partial_residuals
        endif

        MOSE_domain%blk(b)%R(:,i,j,k) = MOSE_domain%blk(b)%R(:,i,j,k) - partial_residuals

      enddo ; enddo ; enddo
      !$OMP END DO

    enddo

  end subroutine IGLOOxMOSE

  !> Derived gas state (rhoi, vel, rho, R, gam, T, mu, kappa) from one MOSE primitive column.
  subroutine Gas_State_From_Prim ( prim, rhoi, vel, rho, R, gam, T, mu, kappa )
    use MOSE_Global_m,           only: nsc, nu, nw, np
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    implicit none
    real(R8), intent(in)  :: prim(:)
    real(R8), intent(out) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa

    rhoi = prim(1:nsc)
    vel  = prim(nu:nw)
    call co_rotot_Rtot ( rhoi, rho, R )
    gam = f_gamma ( rhoi, prim(np), rho, R )
    T   = EOS ( prim(np), rho, R )
    call co_k_mi_lam_Wilke ( rhoi, rho, T, mu, kappa )

  end subroutine Gas_State_From_Prim


  !> Build the ORION gas field IGLOO ingests at setup from MOSE's domain: x, y, z, rho(1..nsc),
  !  u, v, w, T, mil, kl, gam, R (import_gas maps varname v to vars(v-3)).
  subroutine Gas_Field_For_IGLOO ( MOSE_domain, MOSE_IOfield, gas_for_igloo )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: nsc
    use Lib_ORION_data,          only: orion_data, copyORION
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    type(orion_data),       intent(in)  :: MOSE_IOfield
    type(orion_data),       intent(out) :: gas_for_igloo
    ! Local
    integer  :: b, i, j, k, v, nv, ni, nj, nk
    real(R8) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa

    ! start from MOSE's field: mesh, block names and format options
    call copyORION ( MOSE_IOfield, gas_for_igloo )

    nv = nsc + 8   ! rho(1..nsc), u, v, w, T, mil, kl, gam, R

    if (allocated(gas_for_igloo%varnames)) deallocate(gas_for_igloo%varnames)
    allocate(gas_for_igloo%varnames(3+nv))
    gas_for_igloo%varnames(1) = 'x'
    gas_for_igloo%varnames(2) = 'y'
    gas_for_igloo%varnames(3) = 'z'
    do v = 1, nsc
      write(gas_for_igloo%varnames(3+v),'(A,I0,A)') 'rho(', v, ')'
    enddo
    gas_for_igloo%varnames(3+nsc+1) = 'u'
    gas_for_igloo%varnames(3+nsc+2) = 'v'
    gas_for_igloo%varnames(3+nsc+3) = 'w'
    gas_for_igloo%varnames(3+nsc+4) = 'T'
    gas_for_igloo%varnames(3+nsc+5) = 'mil'
    gas_for_igloo%varnames(3+nsc+6) = 'kl'
    gas_for_igloo%varnames(3+nsc+7) = 'gam'
    gas_for_igloo%varnames(3+nsc+8) = 'R'

    do b = 1, size(MOSE_domain%blk)
      ni = MOSE_domain%blk(b)%dim(1)
      nj = MOSE_domain%blk(b)%dim(2)
      nk = MOSE_domain%blk(b)%dim(3)

      if (allocated(gas_for_igloo%block(b)%vars)) deallocate(gas_for_igloo%block(b)%vars)
      allocate(gas_for_igloo%block(b)%vars(1:nv,1:ni,1:nj,1:nk))

      !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(rhoi,vel,rho,kappa,mu,R,gam,T,i,j,k)
      !$OMP DO COLLAPSE(3)
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni

            call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), &
                                       rhoi, vel, rho, R, gam, T, mu, kappa )

            gas_for_igloo%block(b)%vars(1:nsc      ,i,j,k) = rhoi
            gas_for_igloo%block(b)%vars(nsc+1:nsc+3,i,j,k) = vel
            gas_for_igloo%block(b)%vars(nsc+4      ,i,j,k) = T
            gas_for_igloo%block(b)%vars(nsc+5      ,i,j,k) = mu
            gas_for_igloo%block(b)%vars(nsc+6      ,i,j,k) = kappa
            gas_for_igloo%block(b)%vars(nsc+7      ,i,j,k) = gam
            gas_for_igloo%block(b)%vars(nsc+8      ,i,j,k) = R

          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo

    write(*,'(A,I0,A,I0,A)') ' [MI2] built gas field for IGLOO: ', nsc, &
      ' species, ', nv, ' variables (rho(i), u, v, w, T, mil, kl, gam, R)'

  end subroutine Gas_Field_For_IGLOO

  !> Push the MOSE gas state (interior only: IGLOO rebuilds its own ghost ring in solve) into
  !  IGLOO's gasblock and re-impose the gas inflow mass flow mdotGas [kg/s] on the inlet faces, every
  !  sweep. Opens its own parallel region.
  subroutine MOSExIGLOO ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc, gc, MOSE_phase_prefix
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(obj_IGLOO),        intent(inout) :: IGLOO
    ! Local
    integer  :: b, i, j, k, f, m, n, lower, upper
    integer  :: i0, i1, j0, j1, k0, k1, ni, nj, nk
    real(R8) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa
    real(R8) :: mSum, mdot
    integer  :: nSum, noFlow

    upper = 0
    mdotHanded = 0._R8
    nHanded    = 0
    do b = 1, size(MOSE_domain%blk)

      ! mesh-consistency check (not the copy range): the gasblock, ghost ring included under
      ! gas-order = 2, must lie within MOSE's ghosted extent
      associate ( gb => IGLOO%gasblock(b) )
      i0 = lbound(gb%temperature,1); i1 = ubound(gb%temperature,1)
      j0 = lbound(gb%temperature,2); j1 = ubound(gb%temperature,2)
      k0 = lbound(gb%temperature,3); k1 = ubound(gb%temperature,3)
      if ( i0 < 1-gc .or. i1 > MOSE_domain%blk(b)%dim(1)+gc .or. &
           j0 < 1-gc .or. j1 > MOSE_domain%blk(b)%dim(2)+gc .or. &
           k0 < 1-gc .or. k1 > MOSE_domain%blk(b)%dim(3)+gc ) then
        write(*,'(A,I0,A)') ' [MI2] MOSExIGLOO: IGLOO gasblock ', b, &
          ' extends beyond MOSE''s ghosted range -- the two meshes disagree.'
        error stop 1
      endif

      ! interior only: IGLOO rebuilds its own ghost ring in solve (fillGhostGradient) before reading it
      ni = MOSE_domain%blk(b)%dim(1); nj = MOSE_domain%blk(b)%dim(2); nk = MOSE_domain%blk(b)%dim(3)
      !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(rhoi,vel,rho,kappa,mu,R,gam,T,i,j,k)
      !$OMP DO COLLAPSE(3)
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni

            call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), &
                                       rhoi, vel, rho, R, gam, T, mu, kappa )

            gb%density    (1:nsc,i,j,k) = rhoi
            gb%velocity   (1:3  ,i,j,k) = vel
            gb%temperature(      i,j,k) = T
            gb%gam        (      i,j,k) = gam
            gb%R          (      i,j,k) = R
            gb%mil        (      i,j,k) = mu
            gb%kl         (      i,j,k) = kappa

          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      end associate

      ! gas inflow mass flow on the inlet faces (MOSE_domain%bc is flat over blocks and faces);
      ! other records carry no mass flow and leave IGLOO's value alone. A solid-propellant grain in its
      ! inflow mode (ATLAS 502, grain-wall = F: a 404 with srm_inflow) hands the gas mass flux MOSE last
      ! computed at the face, bc%g_srm = mdot_ref (p/p_ref)^n at the pressure of the boundary cell (its
      ! mdot_ref already carries 1 - krho_tot, ATLAS), not the flux at p_ref: IGLOO then injects
      ! krho/(1 - krho_tot) times it.
      do f = 1, 6
        lower = upper + 1
        upper = upper + MOSE_domain%n_bf(b,f)
        ! the walk takes the next n_bf(b,f) records as those of block b, face f: MOSE keeps them in the order of its file
        do i = lower, upper
          if (MOSE_domain%bc(i)%b == b .and. MOSE_domain%bc(i)%f == f) cycle
          write(*,'(A,I0,3A,I0,A,I0,A,I0,A,I0,A)') ' [ERROR] MI2: record ', i, ' of INPUT/', trim(MOSE_phase_prefix), &
            'bc.txt is of block ', MOSE_domain%bc(i)%b, ' face ', MOSE_domain%bc(i)%f, ' where the coupler expects block ', &
            b, ' face ', f, ': the gas boundary records must be grouped by block, then by face (1 to 6), in that order'
          error stop 1
        enddo
        mSum = 0._R8; nSum = 0; noFlow = 0
        !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i,m,n,mdot)
        !$OMP DO REDUCTION(+:mSum,nSum) REDUCTION(max:noFlow)
        do i = lower, upper
          call IGLOO%geoblock(b)%ijk2fmn ( MOSE_domain%bc(i)%i, MOSE_domain%bc(i)%j, &
                                           MOSE_domain%bc(i)%k, f, m, n )
          if (m<=0 .or. n<=0) cycle
          select case (MOSE_domain%bc(i)%type)
          case (401, 402, 405)             ! mdot is the face mass flow MOSE computes [kg/s] along the solver normal:
            mdot = abs(MOSE_domain%bc(i)%mdot)   ! negative for an inflow on faces 2, 4 and 6
          case (406, 407)
            mdot = MOSE_domain%bc(i)%mdot
          case (403, 404, 420)             ! mdot is the prescribed mass flux [kg/(s m2)], an inflow whatever its sign,
            if (MOSE_domain%bc(i)%srm_inflow) then   ! or for a grain the burn law's at the cell pressure (g_srm >= 0)
              mdot = abs(MOSE_domain%bc(i)%g_srm) * IGLOO%geoblock(b)%face(f)%cell(m,n)%area
            else
              mdot = abs(MOSE_domain%bc(i)%mdot) * IGLOO%geoblock(b)%face(f)%cell(m,n)%area
            end if
          case default                     ! no mass flow in the record: a particle inlet that scales with it is refused below
            if (IGLOO%geoblock(b)%face(f)%cell(m,n)%bcdef == 401) noFlow = max(noFlow, i)
            cycle
          end select
          IGLOO%geoblock(b)%face(f)%cell(m,n)%mdotGas = mdot
          mSum = mSum + mdot
          nSum = nSum + 1
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        if (noFlow > 0) then
          write(*,'(A,I0,A,I0,A,I0,A)') ' [ERROR] MI2: block ', b, ' face ', f, ': particles injected in proportion to '// &
            'the gas mass flow (particle BC 401) on a gas boundary of kind ', MOSE_domain%bc(noFlow)%type, &
            ', whose record holds no mass flow'
          error stop 1
        endif
        mdotHanded = mdotHanded + mSum
        nHanded    = nHanded + nSum
      enddo

    enddo

  end subroutine MOSExIGLOO

  !> Push the MOSE gas state into ICE's gas_phase: prim = (rho, u, v, w, T), interior only, plus
  !  R, gam, mu, k. Opens its own parallel region.
  subroutine MOSExICE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    ! Local
    integer  :: b, i, j, k
    real(R8) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa

    do b = 1, size(MOSE_domain%blk)

      !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(rhoi,vel,rho,kappa,mu,R,gam,T,i,j,k)
      !$OMP DO COLLAPSE(3)
      do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
          do i = 1, MOSE_domain%blk(b)%dim(1)

            call Gas_State_From_Prim ( MOSE_domain%blk(b)%P(:,i,j,k), &
                                       rhoi, vel, rho, R, gam, T, mu, kappa )

            ICE_domain%blk(b)%gas_phase%prim(1  ,i,j,k) = rho
            ICE_domain%blk(b)%gas_phase%prim(2:4,i,j,k) = vel
            ICE_domain%blk(b)%gas_phase%prim(5  ,i,j,k) = T
            ICE_domain%blk(b)%gas_phase%R   (     i,j,k) = R
            ICE_domain%blk(b)%gas_phase%gam (     i,j,k) = gam
            ICE_domain%blk(b)%gas_phase%mu  (     i,j,k) = mu
            ICE_domain%blk(b)%gas_phase%k   (     i,j,k) = kappa

          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    enddo

  end subroutine MOSExICE


  !> Stop when a family whose material has no record in the products file loses mass (cond_phase%source(1) not
  !  zero): ICExMOSE would not hand that mass to the gas. Called once per coupling step, outside any parallel region.
  subroutine Check_ICE_Mass_Channel ( ICE_domain )
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of, ICE_phase_prefix
    use hydraMI2_Stoichiometry, only: has_stoichiometry_ICE, mappedICE
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer  :: b, g
    real(R8) :: peak

    do g = 1, ngroups
      if (has_stoichiometry_ICE) then
        if (mappedICE(mat_of(g))) cycle
      endif
      peak = 0._R8
      do b = 1, size(ICE_domain%blk)
        if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%source)) cycle
        peak = max( peak, maxval(abs(ICE_domain%blk(b)%cond_phase(g)%source(1,:,:,:))) )
      enddo
      if (peak > 0._R8) then
        write(*,'(A,I0,A,I0,3A,ES12.4,A)') ' [ERROR] MI2: ICE family ', g, ' (material ', mat_of(g), &
          ') loses mass and INPUT/', trim(ICE_phase_prefix), 'products.txt has no record for its material (peak of '// &
          'cond_phase%source(1) ', peak, ' kg/(m3 s)): the gas would not receive it'
        error stop 1
      endif
    enddo

  end subroutine Check_ICE_Mass_Channel


  !> Stop on a value that is not finite in what ICE hands the gas, the slots ICExMOSE reads (the momentum and the
  !  energy of every family; the mass of a material that converts): MOSE's conservative update
  !  keeps the old temperature of a cell whose energy residual is NaN, so such a value would act unseen. Called once per
  !  coupling step, outside any parallel region.
  subroutine Check_ICE_Source_Finite ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Finite_m,          only: Not_Finite
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: has_stoichiometry_ICE, mappedICE
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(ICE_domain_type),  intent(in) :: ICE_domain
    integer :: b, g, i, j, k, n, s
    logical :: mapped
    character(len=16) :: val

    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%source)) cycle
        mapped = .false.
        if (has_stoichiometry_ICE) mapped = mappedICE(mat_of(g))
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          do n = 1, 4                        ! momentum source(2:4), energy source(nbase(g)-1)
            s = merge(n+1, nbase(g)-1, n < 4)
            if (.not. Not_Finite(ICE_domain%blk(b)%cond_phase(g)%source(s,i,j,k))) cycle
            write(val,'(ES12.5)') ICE_domain%blk(b)%cond_phase(g)%source(s,i,j,k)
            write(*,'(A,I0,A,I0,A,3(1X,I0),A,I0,4A)') ' [ERROR] MI2: ICE handed a non-finite source (block ', b, &
              ', group ', g, ', cell', i, j, k, ', source(', s, '), the ', trim(merge('momentum', 'energy  ', n < 4)), &
              ') = ', trim(adjustl(val))
            error stop 1
          enddo
          if (.not. mapped) cycle
          if (Not_Finite(ICE_domain%blk(b)%cond_phase(g)%source(1,i,j,k))) then
            write(val,'(ES12.5)') ICE_domain%blk(b)%cond_phase(g)%source(1,i,j,k)
            write(*,'(A,I0,A,I0,A,3(1X,I0),2A)') ' [ERROR] MI2: ICE handed a non-finite source (block ', b, &
              ', group ', g, ', cell', i, j, k, ', source(1), the mass) = ', trim(adjustl(val))
            error stop 1
          endif
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine Check_ICE_Source_Finite


  !> Stop on a value that is not finite in the source-mode rates IGLOO hands the gas, those IGLOOxMOSE reads (sourceMom,
  !  sourceEn, and sourceMass under a stoichiometry), for the reason of Check_ICE_Source_Finite. Called once per coupling
  !  step, after IGLOO's solve, outside any parallel region.
  subroutine Check_IGLOO_Source_Finite ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Finite_m,          only: Not_Finite
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nsm
    use hydraMI2_Stoichiometry, only: has_stoichiometry
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    integer :: b, i, j, k, n

    do b = 1, size(MOSE_domain%blk)
      associate ( src => IGLOO%source(b) )
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        if (allocated(src%sourceMom)) then
          do n = 1, 3
            if (Not_Finite(src%sourceMom(n,i,j,k))) call Refuse('sourceMom', n, src%sourceMom(n,i,j,k))
          enddo
        endif
        if (allocated(src%sourceEn)) then
          if (Not_Finite(src%sourceEn(i,j,k))) call Refuse('sourceEn', 0, src%sourceEn(i,j,k))
        endif
        if (has_stoichiometry .and. allocated(src%sourceMass)) then
          do n = 1, nsm()
            if (Not_Finite(src%sourceMass(n,i,j,k))) call Refuse('sourceMass', n, src%sourceMass(n,i,j,k))
          enddo
        endif
      enddo; enddo; enddo
      end associate
    enddo

  contains

    subroutine Refuse ( name, n, v )
      character(len=*), intent(in) :: name
      integer,          intent(in) :: n
      real(R8),         intent(in) :: v
      character(len=16) :: what, val
      what = name
      if (n > 0) write(what,'(2A,I0,A)') name, '(', n, ')'
      write(val,'(ES12.5)') v
      write(*,'(A,I0,A,3(1X,I0),4A)') ' [ERROR] MI2: IGLOO handed a non-finite source (block ', b, ', cell', i, j, k, &
        ', ', trim(what), ') = ', trim(adjustl(val))
      error stop 1
    end subroutine Refuse

  end subroutine Check_IGLOO_Source_Finite


  !> Stop on a value that is not finite in the Eulerian fields IGLOO hands the gas under IGLOO-iter = 2, those IGLOOxMOSE
  !  turns into momentum and energy sources (velocity, temperature, density, np of every family), for the reason of
  !  Check_ICE_Source_Finite. Called once per coupling step, after IGLOO's solve, outside any parallel region.
  subroutine Check_IGLOO_Euler_Finite ( MOSE_domain, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Finite_m,         only: Not_Finite
    use IGLOO_module,          only: obj_IGLOO
    use IGLOO_variables,       only: nfam
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    integer :: b, m, i, j, k, n

    do b = 1, size(MOSE_domain%blk)
      do m = 1, nfam
        associate ( e => IGLOO%euler(b,m) )
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          if (allocated(e%velocity)) then
            do n = 1, 3
              if (Not_Finite(e%velocity(n,i,j,k))) call Refuse('velocity', n, e%velocity(n,i,j,k))
            enddo
          endif
          if (allocated(e%temperature)) then
            if (Not_Finite(e%temperature(i,j,k))) call Refuse('temperature', 0, e%temperature(i,j,k))
          endif
          if (allocated(e%density)) then
            if (Not_Finite(e%density(i,j,k))) call Refuse('density', 0, e%density(i,j,k))
          endif
          if (allocated(e%np)) then
            if (Not_Finite(e%np(i,j,k))) call Refuse('np', 0, e%np(i,j,k))
          endif
        enddo; enddo; enddo
        end associate
      enddo
    enddo

  contains

    subroutine Refuse ( name, n, v )
      character(len=*), intent(in) :: name
      integer,          intent(in) :: n
      real(R8),         intent(in) :: v
      character(len=16) :: what, val
      what = name
      if (n > 0) write(what,'(2A,I0,A)') name, '(', n, ')'
      write(val,'(ES12.5)') v
      write(*,'(A,I0,A,I0,A,3(1X,I0),4A)') ' [ERROR] MI2: IGLOO handed a non-finite Eulerian field (block ', b, &
        ', family ', m, ', cell', i, j, k, ', ', trim(what), ') = ', trim(adjustl(val))
      error stop 1
    end subroutine Refuse

  end subroutine Check_IGLOO_Euler_Finite


  !> Inject ICE's reaction on the gas into MOSE's residual, from MOSE's per-RK-stage callback:
  !  runs inside MOSE's parallel region, orphaned !$OMP DO only. cond_phase%source is the per-unit-
  !  volume rate on the condensed phase over ICE's last step (its RK stages weighted as its update
  !  weighs them), so the gas receives -source (times vol): source(2:4) the
  !  momentum, source(nbase-1) the energy, and for a material that converts source(1), the mass,
  !  through its row and with its datum offset.
  subroutine ICExMOSE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nu, nw, np, nprim
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE, mappedICE
    use MOSE_Config_Types_m,    only: obj_chemistry
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(ICE_domain_type),  intent(in)    :: ICE_domain
    ! Local
    integer  :: b, g, i, j, k, imat, isp, nshare
    real(R8) :: partial_residuals(nprim), dm, avail, dtStage
    logical  :: mapped, consumes

    ! the families whose material consumes a gas species share what a cell holds of it
    nshare = 0
    if (has_stoichiometry_ICE) then
      do g = 1, ngroups
        if (mappedICE(mat_of(g))) then
          if (any(nu_ice(mat_of(g),1:nsc) < 0._R8)) nshare = nshare + 1
        endif
      enddo
    endif

    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups

        imat   = mat_of(g)
        mapped = .false.; consumes = .false.
        if (has_stoichiometry_ICE) then
          mapped = mappedICE(imat)
          if (mapped) consumes = any(nu_ice(imat,1:nsc) < 0._R8)
        endif

        !$OMP DO COLLAPSE(3) PRIVATE(partial_residuals,i,j,k,isp,dm,avail,dtStage)
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)

          partial_residuals(:) = 0._R8

          partial_residuals(nu:nw) = - ICE_domain%blk(b)%cond_phase(g)%source(2:4,i,j,k)
          partial_residuals(np)    = - ICE_domain%blk(b)%cond_phase(g)%source(nbase(g)-1,i,j,k)

          ! a material that converts: its mass through its row, the datum offset on the same mass
          if (mapped) then
            dm = - ICE_domain%blk(b)%cond_phase(g)%source(1,i,j,k)
            do isp = 1, nsc
              partial_residuals(isp) = nu_ice(imat,isp)*dm
            enddo
            partial_residuals(np) = partial_residuals(np) + dhDatumICE(imat)*dm
            ! floor guard: a consumed species leaves at most its share of what the cell holds over this stage's
            ! step, dtlocal*strangcoeff (the step MOSE's Newstate applies), as for IGLOO's rates
            if (consumes) then
              dtStage = MOSE_domain%blk(b)%dtlocal(i,j,k) * obj_chemistry%strangcoeff
              if (dtStage > 0._R8) then
                do isp = 1, nsc
                  if (partial_residuals(isp) >= 0._R8) cycle
                  avail = MOSE_domain%blk(b)%P(isp,i,j,k) / (dtStage * real(nshare, R8))
                  if (-partial_residuals(isp) > avail) then
                    partial_residuals(isp) = -avail
                    !$OMP ATOMIC UPDATE
                    nMassLimited = nMassLimited + 1
                  endif
                enddo
              endif
            endif
          endif

          MOSE_domain%blk(b)%R(:,i,j,k) = MOSE_domain%blk(b)%R(:,i,j,k) &
                                        - partial_residuals * MOSE_domain%blk(b)%vol(i,j,k)

        enddo ; enddo ; enddo
        !$OMP END DO

      enddo
    enddo

  end subroutine ICExMOSE


end module hydraMI2_Mod_Coupling
