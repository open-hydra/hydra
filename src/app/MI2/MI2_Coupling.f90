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
  public :: Gas_Total_Energy
  public :: Gas_Mean_Temperature
  public :: Check_ICE_Mass_Channel
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

  !> Latch: the Check_ICE_Mass_Channel warning is printed once.
  logical, private :: ICE_mass_channel_reported = .false.

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


  !> What ICE's sources hand the gas over the interior cells and the families: force F [N] and energy
  !  rate Q [W] (-source*vol, as ICExMOSE applies them), and the range of the relaxation time tau [s].
  subroutine ICE_To_Gas_Totals ( MOSE_domain, ICE_domain, F, Q, taumin, taumax )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use ICE_Global_m,          only: ngroups, nbase
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    type(ICE_domain_type),  intent(in)  :: ICE_domain
    real(R8),               intent(out) :: F(3), Q, taumin, taumax
    integer  :: b, g, i, j, k
    real(R8) :: vol

    F = 0._R8; Q = 0._R8; taumin = huge(1._R8); taumax = -huge(1._R8)
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          associate ( cp => ICE_domain%blk(b)%cond_phase(g) )
            F = F - cp%source(2:4,i,j,k) * vol
            Q = Q - cp%source(nbase(g)-1,i,j,k) * vol
            taumin = min(taumin, cp%tau(i,j,k))
            taumax = max(taumax, cp%tau(i,j,k))
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
    use hydraMI2_Stoichiometry, only: nu_gas => nu, has_stoichiometry, dhDatum
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    type(obj_IGLOO),        intent(in)  :: IGLOO
    real(R8),               intent(out) :: F(3), Q, M
    integer :: b, i, j, k, imat

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
          M = M + sum(nu_gas(imat,1:nsc)) * IGLOO%source(b)%sourceMass(imat,i,j,k)
          Q = Q + dhDatum(imat) * IGLOO%source(b)%sourceMass(imat,i,j,k)
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine IGLOO_To_Gas_Totals


  !> Accumulate the mass handed to the gas over a coupling step of duration dtBlock:
  !  SUM_m nu(m,s)*sourceMass(m) (absolute kg/s, no volume factor). Called once per coupling step.
  subroutine Accumulate_Injected_Mass ( IGLOO, dtBlock )
    use MOSE_Global_m,          only: nsc
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nm
    use hydraMI2_Stoichiometry, only: nu_gas => nu, has_stoichiometry
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    real(R8),        intent(in) :: dtBlock
    integer  :: b, m, s
    real(R8) :: tot

    if (.not. has_stoichiometry) return
    if (.not. allocated(massInjected)) then
      allocate(massInjected(1:nsc)); massInjected = 0._R8
    endif
    if (dtBlock <= 0._R8) return

    do b = 1, size(IGLOO%source)
      if (.not. allocated(IGLOO%source(b)%sourceMass)) cycle
      do m = 1, nm
        tot = sum(IGLOO%source(b)%sourceMass(m,:,:,:))
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
  !  sourceEn + SUM_m dhDatum(m)*sourceMass(m) (absolute W). Called once per coupling step.
  subroutine Accumulate_Injected_Energy ( IGLOO, dtBlock )
    use IGLOO_module,           only: obj_IGLOO
    use IGLOO_variables,        only: nm
    use hydraMI2_Stoichiometry, only: has_stoichiometry, dhDatum
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    real(R8),        intent(in) :: dtBlock
    integer  :: b, m
    real(R8) :: datum

    if (dtBlock <= 0._R8) return

    do b = 1, size(IGLOO%source)
      if (.not. allocated(IGLOO%source(b)%sourceEn)) cycle
      energyInjected = energyInjected + sum(IGLOO%source(b)%sourceEn) * dtBlock
      energyAudited  = .true.

      if (.not. has_stoichiometry) cycle
      if (.not. allocated(IGLOO%source(b)%sourceMass)) cycle
      datum = 0._R8
      do m = 1, nm
        datum = datum + dhDatum(m) * sum(IGLOO%source(b)%sourceMass(m,:,:,:))
      enddo
      energyInjected  = energyInjected  + datum*dtBlock
      energyDatumPart = energyDatumPart + datum*dtBlock
    enddo

  end subroutine Accumulate_Injected_Energy


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
    use hydraMI2_Stoichiometry,  only: nu_gas => nu, has_stoichiometry, dhDatum
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
    real(R8) :: dm, avail, dtStage

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
      !$OMP&                       imat,isp,dm,avail,dtStage)
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
    use MOSE_Global_m,         only: nsc, gc
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(obj_IGLOO),        intent(inout) :: IGLOO
    ! Local
    integer  :: b, i, j, k, f, m, n, lower, upper
    integer  :: i0, i1, j0, j1, k0, k1, ni, nj, nk
    real(R8) :: rhoi(nsc), vel(3), rho, R, gam, T, mu, kappa
    real(R8) :: mSum
    integer  :: nSum

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
      ! other records carry no mass flow and leave IGLOO's value alone
      do f = 1, 6
        lower = upper + 1
        upper = upper + MOSE_domain%n_bf(b,f)
        mSum = 0._R8; nSum = 0
        !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i,m,n)
        !$OMP DO REDUCTION(+:mSum,nSum)
        do i = lower, upper
          call IGLOO%geoblock(b)%ijk2fmn ( MOSE_domain%bc(i)%i, MOSE_domain%bc(i)%j, &
                                           MOSE_domain%bc(i)%k, f, m, n )
          if (m<=0 .or. n<=0) cycle
          associate ( cell => IGLOO%geoblock(b)%face(f)%cell(m,n) )
            select case (MOSE_domain%bc(i)%type)
            case (401, 402, 405, 406, 407)   ! mdot is the face mass flow MOSE computes [kg/s]
              cell%mdotGas = MOSE_domain%bc(i)%mdot
            case (403, 404, 420)             ! mdot is the prescribed mass flux [kg/(s m2)]
              cell%mdotGas = MOSE_domain%bc(i)%mdot * cell%area
            case default
              cycle
            end select
            mSum = mSum + cell%mdotGas
            nSum = nSum + 1
          end associate
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
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


  !> Warn once if ICE's mass channel cond_phase%source(1), which ICExMOSE drops, is no longer
  !  zero. Called once per coupling step, outside any parallel region.
  subroutine Check_ICE_Mass_Channel ( ICE_domain )
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: ngroups
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer  :: b, g
    real(R8) :: peak

    if (ICE_mass_channel_reported) return

    peak = 0._R8
    do b = 1, size(ICE_domain%blk)
      do g = 1, ngroups
        if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%source)) cycle
        peak = max( peak, maxval(abs(ICE_domain%blk(b)%cond_phase(g)%source(1,:,:,:))) )
      enddo
    enddo

    if (peak > 0._R8) then
      ICE_mass_channel_reported = .true.
      write(*,*)
      write(*,'(A,ES12.4)') ' [MI2] WARNING: ICE cond_phase%source(1) is NO LONGER ZERO -- peak ', peak
      write(*,'(A)') '       ICExMOSE drops that channel, so this mass is NOT reaching the gas.'
      write(*,'(A)') '       Before trusting a combusting/evaporating ICE run, give ICExMOSE an'
      write(*,'(A)') '       ICE-group -> gas-species map and route source(1) into R(1:nsc).'
      write(*,'(A)') '       hydraMI2_Stoichiometry''s nu(m,s) is keyed on IGLOO MATERIALS and must'
      write(*,'(A)') '       NOT be reused here -- it would attribute the mass to the wrong species.'
      write(*,*)
    endif

  end subroutine Check_ICE_Mass_Channel


  !> Inject ICE's reaction on the gas into MOSE's residual, from MOSE's per-RK-stage callback:
  !  runs inside MOSE's parallel region, orphaned !$OMP DO only. cond_phase%source is the per-unit-
  !  volume rate on the condensed phase, so the gas receives -source (times vol); source(1) is the
  !  mass (not coupled), source(2:4) the momentum, source(nbase-1) the energy.
  subroutine ICExMOSE ( MOSE_domain, ICE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc, nu, nw, np, nprim
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use ICE_Global_m,          only: ngroups, nbase
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(ICE_domain_type),  intent(in)    :: ICE_domain
    ! Local
    integer  :: b, g, i, j, k
    real(R8) :: partial_residuals(nprim)

    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups

        !$OMP DO COLLAPSE(3) PRIVATE(partial_residuals,i,j,k)
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)

          partial_residuals(:) = 0._R8

          ! mass (source(1)) is not coupled: no ICE-group -> gas-species map
          partial_residuals(nu:nw) = - ICE_domain%blk(b)%cond_phase(g)%source(2:4,i,j,k)
          partial_residuals(np)    = - ICE_domain%blk(b)%cond_phase(g)%source(nbase(g)-1,i,j,k)

          MOSE_domain%blk(b)%R(:,i,j,k) = MOSE_domain%blk(b)%R(:,i,j,k) &
                                        - partial_residuals * MOSE_domain%blk(b)%vol(i,j,k)

        enddo ; enddo ; enddo
        !$OMP END DO

      enddo
    enddo

  end subroutine ICExMOSE


end module hydraMI2_Mod_Coupling
