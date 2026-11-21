module hydraMI2_Mod_Coupling
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use IGLOO_data_block,     only: obj_sourceblock, obj_eulerblock
  use ICE_Advanced_Types_m, only: ICE_domain_type

  implicit none
  private

  interface
    !> libc rename(2): the hand-off's temporary file onto the hand-off, in one step
    function c_rename ( from, to ) bind(C, name='rename') result( rc )
      use, intrinsic :: iso_c_binding, only: c_char, c_int
      character(kind=c_char), intent(in) :: from(*), to(*)
      integer(c_int) :: rc
    end function c_rename
  end interface

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
  public :: Gone_Delivered
  public :: Gone_Stored
  public :: IGLOO_Smoke_Accumulate
  public :: IGLOO_Smoke_Deliver
  public :: Gone_Clock
  public :: Source_Clock, Source_Mean
  public :: Gone_To_Rate
  public :: ICE_Exchange_File, Write_ICE_Exchange, Read_ICE_Exchange
  public :: Restrict_IGLOO_Sources, Restrict_ICE_Sources

  !> Under-relaxation store: the previous coupling step's injected residual, per block.
  type, public :: resField
    real(R8), dimension(:,:,:,:), allocatable :: resBlock
  end type resField

  !> A value per cell of a block (interior cells).
  type :: vol_t
    real(R8), allocatable :: v(:,:,:)
  end type vol_t

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

  !> The share of the removals ICE stores (cond_phase%gone) that each MOSE iteration of the current block delivers to
  !  the gas: 1/MOSE-iter of the cycle while a block runs, 0 outside one.
  real(R8), public :: goneShare = 0._R8

  !> The smoke of IGLOO's burning materials: per material, for the smoke of its combustion record (1) and the smoke
  !  its cap oxide becomes when handed at burnout (2), the ICE family it is born in (0: none), kg per kg, the smoke's
  !  datum offset; the material's d-smoke. smokeIGOn: one at least. accIG: the amounts a MOSE block gives the smoke
  !  families, per block (6, ICE family, cell), added to ICE's state after the block.
  type, public :: igloo_smoke_t
    integer  :: fam(2) = 0
    real(R8) :: nu(2) = 0._R8, dh(2) = 0._R8, dSmoke = 0._R8
  end type igloo_smoke_t
  type(igloo_smoke_t), allocatable, public :: smokeIG(:)
  logical, public :: smokeIGOn = .false.
  type :: smoke_acc_t
    real(R8), allocatable :: a(:,:,:,:,:)
  end type smoke_acc_t
  type(smoke_acc_t), allocatable :: accIG(:)
  real(R8), allocatable :: dhFamIG(:)
  !> Cells where a smoke was born and no burning parcel's columns gave its state (it took u = 0, T = T_datum_ref).
  integer, public :: nSmokeNoState = 0
  !> What the coupler hands between ICE and the gas or IGLOO, set once after ICE's set-up (hydra_MI2): .true. (MOSE or
  !  ICE steps in pseudo-time) rates, so that a converged steady field holds per unit time whatever the local steps:
  !  IGLOO's smoke as a source each smoke family integrates over ICE's own steps (cond_phase%ext), ICE's removals as
  !  the rate of their store (goneRate); .false. (both time-accurate) amounts: IGLOO's smoke of each MOSE block added to
  !  ICE's state after it, ICE's removals delivered whole over the next MOSE block (goneShare).
  logical, public :: handRates = .false.
  !> The removals of a steady run as rates (handRates), per block: a(5, ICE family, cell) the store ICE made over its
  !  last block divided by the time it was made over, per unit volume and time (mass, momentum (3), energy), kept while
  !  ICE does not step, as ICE's source is; t(ICE family, cell) that time, the sum of the family's dt over the ICE
  !  iterations since the store was last turned into a rate.
  type :: gone_rate_t
    real(R8), allocatable :: a(:,:,:,:,:), t(:,:,:,:)
  end type gone_rate_t
  type(gone_rate_t), allocatable :: goneRate(:)
  !> A time-accurate ICE block of more than one step: per block and family, the sums over its steps of dt*source,
  !  dt*heat and dt, cell by cell; turned into the block's rates after it (Source_Mean) and cleared.
  type :: src_sum_t
    real(R8), allocatable :: s(:,:,:,:), h(:,:,:), t(:,:,:), b(:,:,:,:)
  end type src_sum_t
  type(src_sum_t), allocatable :: srcSum(:,:)
  integer :: srcSteps = 0

  !> The condensed phases' sources restricted from the finest grid to the coarser level MOSE steps, with the sub-cells
  !  per coarse cell along i, j, k (ratio) and their volume (vol).
  type, public :: coarse_src_t
    private
    integer :: level = 0
    integer,               allocatable :: ratio(:,:)
    type(vol_t),           allocatable :: vol(:)
    type(obj_sourceblock), allocatable :: src(:)
    type(obj_eulerblock),  allocatable :: eul(:,:)
    type(ICE_domain_type)              :: ice
    type(gone_rate_t),     allocatable :: gr(:)
  end type coarse_src_t

  !> Inflow mass flow MOSExIGLOO handed IGLOO at its last call [kg/s], and the number of face records.
  real(R8), public :: mdotHanded = 0._R8
  integer,  public :: nHanded    = 0

  !> A variable of the file a run with ICE ends with (ICE_Exchange_File): its name, what it holds (XS source, XH heat,
  !  XB birth, XG the removals' store, XU the relaxation time, XR a smoke family's srate, XA and XT goneRate's a and t),
  !  family, component.
  type :: xvar_t
    character(len=32) :: name = ''
    integer :: what = 0, g = 0, n = 0
  end type xvar_t
  integer, parameter :: XS = 1, XH = 2, XB = 3, XG = 4, XR = 5, XA = 6, XT = 7, XU = 8
  character(len=24), parameter :: XAUX(2) = [character(len=24) :: 'MI2_ICE_HAS_SOURCES', 'MI2_MOSE_ITERATION']

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
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE, dhSmokeICE
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
            F = F - (cp%source(2:4,i,j,k) + cp%birth(2:4,i,j,k)) * vol
            Q = Q - (cp%source(nbase(g)-1,i,j,k) + cp%birth(5,i,j,k)) * vol
            taumin = min(taumin, cp%tau(i,j,k))
            taumax = max(taumax, cp%tau(i,j,k))
            if (has_stoichiometry_ICE) then
              dm = - cp%source(1,i,j,k) * vol
              M  = M + sum(nu_ice(imat,1:nsc)) * dm
              Q  = Q + dhDatumICE(imat) * dm - dhSmokeICE(imat) * cp%birth(1,i,j,k) * vol
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
    real(R8) :: cr(4), srate(6,2)

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
        if (.not. smokeIGOn) cycle
        do imat = 1, nm
          call IGLOO_Smoke_Rates( IGLOO%source(b)%sourceMass(:,i,j,k), imat, srate )
          F = F - srate(2:4,1) - srate(2:4,2)
          Q = Q - srate(5,1) - srate(5,2) - smokeIG(imat)%dh(1)*srate(1,1) - smokeIG(imat)%dh(2)*srate(1,2)
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


  !> The smoke an IGLOO material's burning makes in a cell, as rates [per s, the cell's]: (:,1) the smoke of its
  !  combustion record, (:,2) the smoke its cap oxide handed at burnout becomes; mass, momentum (3), energy on ICE's
  !  datum with the kinetic part, number (each drop at d-smoke); born with the parcels' velocity and temperature
  !  weighted by what they burnt there (IGLOO's columns). src: the cell's column of IGLOO's mass source.
  subroutine IGLOO_Smoke_Rates ( src, imat, rate, nostate )
    use IGLOO_variables,        only: nm, smokeRow
    use ICE_Config_Types_m,     only: ICE_condensed => obj_condensed
    use ICE_Global_m,           only: ICE_mat_of => mat_of
    use ICE_Lib_Properties,     only: mat_rho, mat_e, mat_cp_const
    use hydraMI2_Stoichiometry, only: capped, T_datum_ref
    implicit none
    real(R8), intent(in)  :: src(:)
    integer,  intent(in)  :: imat
    real(R8), intent(out) :: rate(6,2)
    logical,  intent(out), optional :: nostate      !> a smoke born where no parcel's columns gave its state
    real(R8), parameter :: pi = 3.14159265358979323846_R8
    real(R8) :: md(2), w, u(3), T, es, cr(4)
    integer  :: n

    rate = 0._R8
    if (present(nostate)) nostate = .false.
    if (all(smokeIG(imat)%fam == 0)) return
    if (capped(imat)) then
      cr = Cap_Rates( imat, src(imat), src(nm+2*imat-1), src(nm+2*imat) )
      md = [smokeIG(imat)%nu(1)*cr(2), smokeIG(imat)%nu(2)*cr(4)]
    else
      md = [smokeIG(imat)%nu(1)*src(imat), 0._R8]
    endif
    if (all(md == 0._R8)) return
    w = src(smokeRow(imat,1))
    if (w > 0._R8) then
      u = src(smokeRow(imat,2):smokeRow(imat,4))/w
      T = src(smokeRow(imat,5))/w
    else
      u = 0._R8
      T = T_datum_ref
      if (present(nostate)) nostate = .true.
    endif
    do n = 1, 2
      if (md(n) == 0._R8 .or. smokeIG(imat)%fam(n) == 0) cycle
      associate (smat => ICE_condensed(ICE_mat_of(smokeIG(imat)%fam(n))))
        if (smat%cs_varies) then
          es = mat_e(smat, T)
        else
          es = mat_cp_const(smat)*T
        endif
        rate(1,n)   = md(n)
        rate(2:4,n) = md(n)*u
        rate(5,n)   = md(n)*(es + 0.5_R8*dot_product(u, u))
        rate(6,n)   = md(n)/(pi/6._R8*mat_rho(smat, T)*smokeIG(imat)%dSmoke**3)
      end associate
    enddo

  end subroutine IGLOO_Smoke_Rates


  !> After a MOSE iteration on the level it stepped: IGLOO's smoke rates times dtlocal, added to the block's amounts (a
  !  coarse cell's to its first sub-cell: only totals are read in a steady run). Outside any parallel region.
  subroutine IGLOO_Smoke_Accumulate ( MOSE_domain, IGLOO, fine, cs )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain, fine
    type(obj_IGLOO),        intent(in) :: IGLOO
    type(coarse_src_t),     intent(in) :: cs
    integer :: one(3, size(MOSE_domain%blk))

    if (.not. smokeIGOn) return
    call Smoke_Acc_Setup( fine )
    if (MOSE_domain%mg_level == 1) then
      one = 1
      call Smoke_Accumulate_Cells( MOSE_domain, IGLOO%source, one )
    else
      call Smoke_Accumulate_Cells( MOSE_domain, cs%src, cs%ratio )
    endif

  end subroutine IGLOO_Smoke_Accumulate


  !> The cells of IGLOO_Smoke_Accumulate, from the source columns src on the cells of MOSE_domain; a cell's amounts go to
  !  the finest grid's cell (i-1)*ratio+1 along each direction.
  subroutine Smoke_Accumulate_Cells ( MOSE_domain, src, ratio )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_variables,       only: nm
    implicit none
    type(MOSE_domain_type),             intent(in) :: MOSE_domain
    type(obj_sourceblock), allocatable, intent(in) :: src(:)
    integer,                            intent(in) :: ratio(:,:)
    integer  :: b, i, j, k, imat, n, fi, fj, fk
    real(R8) :: srate(6,2), dt
    logical  :: nostate

    do b = 1, size(MOSE_domain%blk)
      if (.not. allocated(src(b)%sourceMass)) cycle
      do k = 1, MOSE_domain%blk(b)%dim(3)
      do j = 1, MOSE_domain%blk(b)%dim(2)
      do i = 1, MOSE_domain%blk(b)%dim(1)
        fi = (i-1)*ratio(1,b) + 1
        fj = (j-1)*ratio(2,b) + 1
        fk = (k-1)*ratio(3,b) + 1
        dt = MOSE_domain%blk(b)%dtlocal(i,j,k)
        do imat = 1, nm
          call IGLOO_Smoke_Rates( src(b)%sourceMass(:,i,j,k), imat, srate, nostate )
          if (nostate) nSmokeNoState = nSmokeNoState + 1
          do n = 1, 2
            if (srate(1,n) == 0._R8) cycle
            accIG(b)%a(:,smokeIG(imat)%fam(n),fi,fj,fk) = accIG(b)%a(:,smokeIG(imat)%fam(n),fi,fj,fk) + srate(:,n)*dt
          enddo
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine Smoke_Accumulate_Cells


  !> The amounts of the MOSE blocks, per block (6, ICE family, cell), and each family's datum offset: allocated at the
  !  first call, emptied.
  subroutine Smoke_Acc_Setup ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Global_m,          only: ICE_ngroups => ngroups
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer :: b, imat, n

    if (allocated(accIG)) return
    allocate(accIG(size(MOSE_domain%blk)), dhFamIG(ICE_ngroups))
    dhFamIG = 0._R8
    do imat = 1, size(smokeIG)
      do n = 1, 2
        if (smokeIG(imat)%fam(n) > 0) dhFamIG(smokeIG(imat)%fam(n)) = smokeIG(imat)%dh(n)
      enddo
    enddo
    do b = 1, size(MOSE_domain%blk)
      allocate(accIG(b)%a(6, ICE_ngroups, MOSE_domain%blk(b)%dim(1), MOSE_domain%blk(b)%dim(2), &
                          MOSE_domain%blk(b)%dim(3)))
      accIG(b)%a = 0._R8
    enddo

  end subroutine Smoke_Acc_Setup


  !> After a MOSE block: the amounts of smoke it gave each family (IGLOO_Smoke_Accumulate), which the gas did not get;
  !  M, P, E their totals (E on the gas's formation scale), taken out of what the energy audit says the gas received;
  !  the amounts emptied. The hand-off, in its form (handRates): the amounts added to ICE's state cell by cell (per
  !  volume); or the rates of IGLOO's sources as they stand (srcOn; else 0), per volume, set as the source of each
  !  smoke family (cond_phase%ext), R their totals (mass, momentum (3), energy on the formation scale, number; per
  !  second). Outside any parallel region.
  subroutine IGLOO_Smoke_Deliver ( MOSE_domain, ICE_domain, IGLOO, srcOn, M, P, E, R )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    use IGLOO_variables,       only: nm
    use ICE_Global_m,          only: ngroups, mat_of, nbase, ncond
    use ICE_Config_Types_m,    only: ICE_condensed => obj_condensed
    use ICE_Lib_Model,         only: prim_2_cons, cons_2_prim, assign_prim_2_cons, assign_cons_2_prim
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: srcOn
    real(R8),               intent(out)   :: M, P(3), E, R(6)
    integer  :: b, g, i, j, k, ns, imat, n
    real(R8) :: vol, ac(6), Ed, srate(6,2), sr(6)
    real(R8), allocatable :: c(:)

    M = 0._R8; P = 0._R8; E = 0._R8; Ed = 0._R8; R = 0._R8
    if (.not. smokeIGOn) return
    call Smoke_Acc_Setup( MOSE_domain )
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        if (.not. any(smokeIG(:)%fam(1) == g .or. smokeIG(:)%fam(2) == g)) cycle
        ! the conversions of family g's closure (ICE points them at the family it last stepped)
        call assign_prim_2_cons(g)
        call assign_cons_2_prim(g)
        ns = nbase(g)
        if (allocated(c)) deallocate(c)
        allocate(c(ncond(g)))
        if (handRates) then
          if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%ext)) &
            allocate(ICE_domain%blk(b)%cond_phase(g)%ext(ncond(g), ICE_domain%blk(b)%dim(1), &
                     ICE_domain%blk(b)%dim(2), ICE_domain%blk(b)%dim(3)))
          ICE_domain%blk(b)%cond_phase(g)%ext = 0._R8
        endif
        associate (smat => ICE_condensed(mat_of(g)))
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          if (handRates .and. srcOn) then
            if (allocated(IGLOO%source(b)%sourceMass)) then
              sr = 0._R8
              do imat = 1, nm
                call IGLOO_Smoke_Rates( IGLOO%source(b)%sourceMass(:,i,j,k), imat, srate )
                do n = 1, 2
                  if (smokeIG(imat)%fam(n) == g) sr = sr + srate(:,n)
                enddo
              enddo
              c = 0._R8
              c(1:4)  = sr(1:4)/vol
              c(ns-1) = sr(5)/vol
              c(ns)   = sr(6)/vol
              ICE_domain%blk(b)%cond_phase(g)%ext(:,i,j,k) = c
              R(1:4) = R(1:4) + sr(1:4)
              R(5)   = R(5) + sr(5) + dhFamIG(g)*sr(1)
              R(6)   = R(6) + sr(6)
            endif
          endif
          ac = accIG(b)%a(:,g,i,j,k)
          if (all(ac == 0._R8)) cycle
          if (.not. handRates) then
            c = prim_2_cons(ICE_domain%blk(b)%cond_phase(g)%prim(:,i,j,k), smat)
            c(1:4)  = c(1:4)  + ac(1:4)/vol
            c(ns-1) = c(ns-1) + ac(5)/vol
            c(ns)   = c(ns)   + ac(6)/vol
            ICE_domain%blk(b)%cond_phase(g)%prim(:,i,j,k) = cons_2_prim(c, smat)
          endif
          M = M + ac(1)
          P = P + ac(2:4)
          E = E + ac(5) + dhFamIG(g)*ac(1)
          Ed = Ed + dhFamIG(g)*ac(1)
          accIG(b)%a(:,g,i,j,k) = 0._R8
        enddo; enddo; enddo
        end associate
      enddo
    enddo
    energyInjected  = energyInjected  - E
    energyDatumPart = energyDatumPart - Ed

  end subroutine IGLOO_Smoke_Deliver


  !> After a MOSE block that delivered the share frac of the removals ICE stores: the books take what was delivered
  !  (mass through the rows, energy with the datum offset), the stores keep the rest; M, P, E are the totals delivered
  !  on the gas's scale. Outside any parallel region.
  subroutine Gone_Delivered ( MOSE_domain, ICE_domain, frac, M, P, E )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE, mappedICE, dhSmokeICE, &
                                      nuSmokeICE
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    real(R8),               intent(in)    :: frac
    real(R8),               intent(out)   :: M, P(3), E
    integer  :: b, g, i, j, k, s, imat
    real(R8) :: vol, dm, gn(5)

    M = 0._R8; P = 0._R8; E = 0._R8
    if (.not. has_stoichiometry_ICE) return
    if (.not. allocated(massInjected)) then
      allocate(massInjected(1:nsc)); massInjected = 0._R8
    endif
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        imat = mat_of(g)
        if (.not. mappedICE(imat)) cycle
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          gn  = ICE_domain%blk(b)%cond_phase(g)%gone(:,i,j,k)
          if (all(gn == 0._R8)) cycle
          dm = gn(1)*frac*vol
          do s = 1, nsc
            massInjected(s) = massInjected(s) + nu_ice(imat,s)*dm
          enddo
          energyInjected  = energyInjected  + gn(5)*frac*vol + (dhDatumICE(imat) - nuSmokeICE(imat)*dhSmokeICE(imat))*dm
          energyDatumPart = energyDatumPart + (dhDatumICE(imat) - nuSmokeICE(imat)*dhSmokeICE(imat))*dm
          energyAudited   = .true.
          M = M + sum(nu_ice(imat,1:nsc))*dm
          P = P + gn(2:4)*frac*vol
          E = E + gn(5)*frac*vol + (dhDatumICE(imat) - nuSmokeICE(imat)*dhSmokeICE(imat))*dm
          ICE_domain%blk(b)%cond_phase(g)%gone(:,i,j,k) = gn*(1._R8 - frac)
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine Gone_Delivered


  !> The arrays of goneRate, allocated at the first call and zero.
  subroutine Gone_Rate_Setup ( ICE_domain )
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: ngroups
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer :: b

    if (allocated(goneRate)) return
    allocate(goneRate(ICE_domain%nb))
    do b = 1, ICE_domain%nb
      allocate(goneRate(b)%a(5, ngroups, ICE_domain%blk(b)%dim(1), ICE_domain%blk(b)%dim(2), ICE_domain%blk(b)%dim(3)), &
               goneRate(b)%t(ngroups, ICE_domain%blk(b)%dim(1), ICE_domain%blk(b)%dim(2), ICE_domain%blk(b)%dim(3)))
      goneRate(b)%a = 0._R8
      goneRate(b)%t = 0._R8
    enddo

  end subroutine Gone_Rate_Setup


  !> After an ICE iteration of a steady run: each converting family's dt added, cell by cell, to the time its store of
  !  removals is being made over. Outside any parallel region.
  subroutine Gone_Clock ( ICE_domain )
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Stoichiometry, only: has_stoichiometry_ICE, mappedICE
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer :: b, g, i, j, k

    if (.not. has_stoichiometry_ICE) return
    call Gone_Rate_Setup( ICE_domain )
    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        if (.not. mappedICE(mat_of(g))) cycle
        if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%gone)) cycle
        do k = 1, ICE_domain%blk(b)%dim(3)
        do j = 1, ICE_domain%blk(b)%dim(2)
        do i = 1, ICE_domain%blk(b)%dim(1)
          goneRate(b)%t(g,i,j,k) = goneRate(b)%t(g,i,j,k) + ICE_domain%blk(b)%cond_phase(g)%dt(i,j,k)
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine Gone_Clock


  !> After a time-accurate ICE step: each family's dt*source, dt*heat, dt*birth (a family that makes a smoke) and dt
  !  added to the block's sums, cell by cell. Outside any parallel region.
  subroutine Source_Clock ( ICE_domain )
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer :: b, g, i, j, k

    if (.not. allocated(srcSum)) allocate(srcSum(ICE_domain%nb, ngroups))
    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        associate (cp => ICE_domain%blk(b)%cond_phase(g), sm => srcSum(b,g))
          if (.not. (allocated(cp%source) .and. allocated(cp%heat))) cycle
          if (.not. allocated(sm%s)) then
            allocate(sm%s, mold=cp%source); allocate(sm%h, mold=cp%heat); allocate(sm%t, mold=cp%heat)
            sm%s = 0._R8; sm%h = 0._R8; sm%t = 0._R8
            if (allocated(cp%birth)) then
              allocate(sm%b, mold=cp%birth); sm%b = 0._R8
            endif
          endif
          do k = 1, ICE_domain%blk(b)%dim(3)
          do j = 1, ICE_domain%blk(b)%dim(2)
          do i = 1, ICE_domain%blk(b)%dim(1)
            sm%s(:,i,j,k) = sm%s(:,i,j,k) + cp%dt(i,j,k)*cp%source(:,i,j,k)
            sm%h(i,j,k)   = sm%h(i,j,k)   + cp%dt(i,j,k)*cp%heat(i,j,k)
            sm%t(i,j,k)   = sm%t(i,j,k)   + cp%dt(i,j,k)
            if (allocated(sm%b)) sm%b(:,i,j,k) = sm%b(:,i,j,k) + cp%dt(i,j,k)*cp%birth(:,i,j,k)
          enddo; enddo; enddo
        end associate
      enddo
    enddo
    srcSteps = srcSteps + 1

  end subroutine Source_Clock


  !> After a time-accurate ICE block: where it made more than one step, each family's source, heat and birth become
  !  the block's rates, the sums of Source_Clock divided by the time; the sums cleared. Outside any parallel region.
  subroutine Source_Mean ( ICE_domain )
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups
    implicit none
    type(ICE_domain_type), intent(inout) :: ICE_domain
    integer :: b, g, i, j, k

    if (.not. allocated(srcSum)) return
    do b = 1, ICE_domain%nb
      do g = 1, ngroups
        associate (cp => ICE_domain%blk(b)%cond_phase(g), sm => srcSum(b,g))
          if (.not. allocated(sm%s)) cycle
          if (srcSteps > 1) then
            do k = 1, ICE_domain%blk(b)%dim(3)
            do j = 1, ICE_domain%blk(b)%dim(2)
            do i = 1, ICE_domain%blk(b)%dim(1)
              if (sm%t(i,j,k) > 0._R8) then
                cp%source(:,i,j,k) = sm%s(:,i,j,k)/sm%t(i,j,k)
                cp%heat(i,j,k)     = sm%h(i,j,k)/sm%t(i,j,k)
                if (allocated(sm%b)) cp%birth(:,i,j,k) = sm%b(:,i,j,k)/sm%t(i,j,k)
              endif
            enddo; enddo; enddo
          endif
          sm%s = 0._R8; sm%h = 0._R8; sm%t = 0._R8
          if (allocated(sm%b)) sm%b = 0._R8
        end associate
      enddo
    enddo
    srcSteps = 0

  end subroutine Source_Mean


  !> At a MOSE block's start in a steady run: where ICE stepped since the last call, the store of removals divided by
  !  the time it was made over becomes the rate the gas receives (goneRate, per unit volume and time), the store and
  !  its time emptied; elsewhere the rate stands. M, P, E: the rates' totals on the gas's scale (mass through the
  !  rows, momentum, energy with the datum offset), per second. Outside any parallel region.
  subroutine Gone_To_Rate ( MOSE_domain, ICE_domain, M, P, E )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE, mappedICE, dhSmokeICE, &
                                      nuSmokeICE
    implicit none
    type(MOSE_domain_type), intent(in)    :: MOSE_domain
    type(ICE_domain_type),  intent(inout) :: ICE_domain
    real(R8),               intent(out)   :: M, P(3), E
    integer  :: b, g, i, j, k, imat
    real(R8) :: vol, dm, tg

    M = 0._R8; P = 0._R8; E = 0._R8
    if (.not. has_stoichiometry_ICE) return
    call Gone_Rate_Setup( ICE_domain )
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        imat = mat_of(g)
        if (.not. mappedICE(imat)) cycle
        if (.not. allocated(ICE_domain%blk(b)%cond_phase(g)%gone)) cycle
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          tg = goneRate(b)%t(g,i,j,k)
          if (tg > 0._R8) then
            goneRate(b)%a(:,g,i,j,k) = ICE_domain%blk(b)%cond_phase(g)%gone(1:5,i,j,k)/tg
            ICE_domain%blk(b)%cond_phase(g)%gone(1:5,i,j,k) = 0._R8
            goneRate(b)%t(g,i,j,k) = 0._R8
          endif
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          dm  = goneRate(b)%a(1,g,i,j,k)*vol
          M   = M + sum(nu_ice(imat,1:nsc))*dm
          P   = P + goneRate(b)%a(2:4,g,i,j,k)*vol
          E   = E + goneRate(b)%a(5,g,i,j,k)*vol + (dhDatumICE(imat) - nuSmokeICE(imat)*dhSmokeICE(imat))*dm
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine Gone_To_Rate


  !> What the removals ICE stores still hold for the gas, on the gas's scale: M the mass through the rows, P the
  !  momentum, E the energy with the datum offset. Outside any parallel region.
  subroutine Gone_Stored ( MOSE_domain, ICE_domain, M, P, E )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE, mappedICE, dhSmokeICE, &
                                      nuSmokeICE
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    type(ICE_domain_type),  intent(in)  :: ICE_domain
    real(R8),               intent(out) :: M, P(3), E
    integer  :: b, g, i, j, k, imat
    real(R8) :: vol, dm, gn(5)

    M = 0._R8; P = 0._R8; E = 0._R8
    if (.not. has_stoichiometry_ICE) return
    do b = 1, size(MOSE_domain%blk)
      do g = 1, ngroups
        imat = mat_of(g)
        if (.not. mappedICE(imat)) cycle
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          gn  = ICE_domain%blk(b)%cond_phase(g)%gone(:,i,j,k)
          dm  = gn(1)*vol
          M   = M + sum(nu_ice(imat,1:nsc))*dm
          P   = P + gn(2:4)*vol
          E   = E + gn(5)*vol + (dhDatumICE(imat) - nuSmokeICE(imat)*dhSmokeICE(imat))*dm
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine Gone_Stored


  !> Accumulate what ICE's sources hand the gas over a coupling step of duration dtBlock, as ICExMOSE applies them:
  !  the mass of each material that converts, through its row, and the energy with its datum offset. Called once
  !  per coupling step when an ICE material converts.
  subroutine Accumulate_Injected_ICE ( MOSE_domain, ICE_domain, dtBlock )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE, dhSmokeICE
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
            etot = etot - (cp%source(nbase(g)-1,i,j,k) + cp%birth(5,i,j,k) + dhSmokeICE(imat)*cp%birth(1,i,j,k)) * vol
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
  !  callback: runs inside MOSE's parallel region, orphaned !$OMP DO only. On MOSE's finest grid IGLOO's own
  !  arrays; on a coarser level the restricted ones of cs.
  subroutine IGLOOxMOSE ( MOSE_domain, IGLOO, eulCoupling, relaxation, iter, iMax, relCoeff, &
                          oldResidual, cs )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use IGLOO_module,            only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: eulCoupling, relaxation
    integer,                intent(in)    :: iter, iMax
    real(R8),               intent(in)    :: relCoeff
    type(resField),         intent(inout) :: oldResidual(:)
    type(coarse_src_t),     intent(in)    :: cs

    if (MOSE_domain%mg_level == 1) then
      call IGLOO_Residual( MOSE_domain, IGLOO, IGLOO%source, IGLOO%euler, eulCoupling, relaxation, iter, iMax, &
                           relCoeff, oldResidual )
    else
      call IGLOO_Residual( MOSE_domain, IGLOO, cs%src, cs%eul, eulCoupling, relaxation, iter, iMax, &
                           relCoeff, oldResidual, cs%vol )
    endif

  end subroutine IGLOOxMOSE


  !> IGLOO's reaction on the gas into MOSE's residual, from the arrays given on the cells of MOSE_domain; volE, when
  !  given, the volume the Eulerian branch's per-volume source is taken over, else the cell's own.
  !    source branch : source{Mass,Mom,En} are absolute rates (no *vol); mass through nu(m,s),
  !                    plus the enthalpy-datum correction and the consumed-species floor guard
  !    euler  branch : getSourceTerms is per unit volume (*vol); momentum and energy only
  subroutine IGLOO_Residual ( MOSE_domain, IGLOO, src, eul, eulCoupling, relaxation, iter, iMax, relCoeff, &
                              oldResidual, volE )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: nsc, nu, nw, np, nprim
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    use IGLOO_module,            only: obj_IGLOO
    use IGLOO_variables,         only: nfam, nm
    use hydraMI2_Stoichiometry,  only: nu_gas => nu, has_stoichiometry, dhDatum, capped, nuCap, nuRel, dhProduct
    use MOSE_Config_Types_m,     only: obj_chemistry
    implicit none
    type(MOSE_domain_type),             intent(inout) :: MOSE_domain
    type(obj_IGLOO),                    intent(in)    :: IGLOO
    type(obj_sourceblock), allocatable, intent(in)    :: src(:)
    type(obj_eulerblock),  allocatable, intent(in)    :: eul(:,:)
    logical,                            intent(in)    :: eulCoupling, relaxation
    integer,                            intent(in)    :: iter, iMax
    real(R8),                           intent(in)    :: relCoeff
    type(resField),                     intent(inout) :: oldResidual(:)
    type(vol_t),          optional,     intent(in)    :: volE(:)
    ! Local
    real(R8), parameter :: toll = 1.d-15
    logical  :: useRelax
    integer  :: b, i, j, k, m
    ! mirrors blk%R(1:nprim): species 1:nsc, momentum nu:nw, energy np, then soot/passive/RANS (left at zero)
    real(R8) :: partial_residuals(nprim), sources(4)
    real(R8) :: rhoi(nsc), rho, u(3), T, p, R, gam, kappa, mu, coeff
    integer  :: imat, isp
    real(R8) :: dm, avail, dtStage, cr(4), srate(6,2), vE

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
      !$OMP&                       imat,isp,dm,avail,dtStage,cr,srate,vE)
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
                                    eul(b,m)%velocity(:,i,j,k),                  &
                                    eul(b,m)%temperature(i,j,k),                 &
                                    eul(b,m)%density(i,j,k),                     &
                                    eul(b,m)%np(i,j,k), famMatID(m) )
          enddo
          partial_residuals(nu:nw) = - sources(1:3)
          partial_residuals(np)    = - sources(4)
          vE = MOSE_domain%blk(b)%vol(i,j,k)
          if (present(volE)) vE = volE(b)%v(i,j,k)
          partial_residuals = partial_residuals * vE

        else

          partial_residuals(nu:nw) = src(b)%sourceMom(:,i,j,k)
          partial_residuals(np)    = src(b)%sourceEn(i,j,k)

          ! mass coupling: sourceMass per material, distributed over the gas species through nu
          if (has_stoichiometry) then
            do imat = 1, nm
              dm = src(b)%sourceMass(imat,i,j,k)
              ! an oxide cap: the material burnt into gas species, the oxide the cap gained, the cap handed to the gas
              if (capped(imat)) then
                cr = Cap_Rates( imat, dm, src(b)%sourceMass(nm+2*imat-1,i,j,k), &
                                src(b)%sourceMass(nm+2*imat,i,j,k) )
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
            ! the smoke of a material that makes one is born in its ICE family: the gas receives IGLOO's books less it
            if (smokeIGOn) then
              do imat = 1, nm
                call IGLOO_Smoke_Rates( src(b)%sourceMass(:,i,j,k), imat, srate )
                partial_residuals(nu:nw) = partial_residuals(nu:nw) - srate(2:4,1) - srate(2:4,2)
                partial_residuals(np)    = partial_residuals(np) - srate(5,1) - srate(5,2) &
                                           - smokeIG(imat)%dh(1)*srate(1,1) - smokeIG(imat)%dh(2)*srate(1,2)
              enddo
            endif

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

  end subroutine IGLOO_Residual

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
  !  runs inside MOSE's parallel region, orphaned !$OMP DO only. On MOSE's finest grid ICE's own arrays; on a
  !  coarser level the restricted ones of cs.
  subroutine ICExMOSE ( MOSE_domain, ICE_domain, cs )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(ICE_domain_type),  intent(in)    :: ICE_domain
    type(coarse_src_t),     intent(in)    :: cs

    if (MOSE_domain%mg_level == 1) then
      call ICE_Residual( MOSE_domain, ICE_domain, goneRate )
    else
      call ICE_Residual( MOSE_domain, cs%ice, cs%gr )
    endif

  end subroutine ICExMOSE


  !> ICE's reaction on the gas into MOSE's residual, from the arrays given, on the cells of MOSE_domain: ICE_domain's
  !  cond_phase and gr, the removals' rates of a steady run (goneRate). cond_phase%source is the per-unit-
  !  volume rate on the condensed phase over ICE's last step (its RK stages weighted as its update
  !  weighs them; over a time-accurate block of more than one step, the steps' weighted by their dt,
  !  Source_Mean), so the gas receives -source (times vol): source(2:4) the
  !  momentum, source(nbase-1) the energy, and for a material that converts source(1), the mass,
  !  through its row and with its datum offset.
  subroutine ICE_Residual ( MOSE_domain, ICE_domain, gr )
    use MOSE_Advanced_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m,          only: nsc, nu, nw, np, nprim
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Global_m,           only: ngroups, nbase, mat_of
    use hydraMI2_Stoichiometry, only: nu_ice => nuICE, has_stoichiometry_ICE, dhDatumICE, mappedICE, dhSmokeICE, &
                                      nuSmokeICE
    use MOSE_Config_Types_m,    only: obj_chemistry
    implicit none
    type(MOSE_domain_type),         intent(inout) :: MOSE_domain
    type(ICE_domain_type),          intent(in)    :: ICE_domain
    type(gone_rate_t), allocatable, intent(in)    :: gr(:)
    ! Local
    integer  :: b, g, i, j, k, imat, isp, nshare
    real(R8) :: partial_residuals(nprim), dm, avail, dtStage, dtIter, dmg, gn(5)
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

        !$OMP DO COLLAPSE(3) PRIVATE(partial_residuals,i,j,k,isp,dm,avail,dtStage,dtIter,dmg,gn)
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
            ! the smoke the family makes is not the gas's: its momentum, its energy and its datum offset stay out
            partial_residuals(nu:nw) = partial_residuals(nu:nw) - ICE_domain%blk(b)%cond_phase(g)%birth(2:4,i,j,k)
            partial_residuals(np)    = partial_residuals(np) - ICE_domain%blk(b)%cond_phase(g)%birth(5,i,j,k) &
                                     - dhSmokeICE(imat)*ICE_domain%blk(b)%cond_phase(g)%birth(1,i,j,k)
            ! the removals ICE stores: a steady run, the rate of the store (goneRate); a time-accurate one, this MOSE
            ! iteration's share of the amounts, as a rate over its step (each stage of it and both Strang halves
            ! integrate a constant rate to rate*dtlocal)
            if (handRates) then
              if (allocated(gr)) then
                gn  = gr(b)%a(:,g,i,j,k)
                dmg = gn(1)
                do isp = 1, nsc
                  partial_residuals(isp) = partial_residuals(isp) + nu_ice(imat,isp)*dmg
                enddo
                partial_residuals(nu:nw) = partial_residuals(nu:nw) + gn(2:4)
                partial_residuals(np)    = partial_residuals(np) + gn(5) &
                                         + (dhDatumICE(imat) - nuSmokeICE(imat)*dhSmokeICE(imat))*dmg
              endif
            elseif (goneShare > 0._R8) then
              dtIter = MOSE_domain%blk(b)%dtlocal(i,j,k)
              if (dtIter > 0._R8) then
                dmg = ICE_domain%blk(b)%cond_phase(g)%gone(1,i,j,k)*goneShare/dtIter
                do isp = 1, nsc
                  partial_residuals(isp) = partial_residuals(isp) + nu_ice(imat,isp)*dmg
                enddo
                partial_residuals(nu:nw) = partial_residuals(nu:nw) &
                                         + ICE_domain%blk(b)%cond_phase(g)%gone(2:4,i,j,k)*goneShare/dtIter
                partial_residuals(np)    = partial_residuals(np) &
                                         + ICE_domain%blk(b)%cond_phase(g)%gone(5,i,j,k)*goneShare/dtIter &
                                         + (dhDatumICE(imat) - nuSmokeICE(imat)*dhSmokeICE(imat))*dmg
              endif
            endif
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

  end subroutine ICE_Residual


  !> cs set to MOSE's level coarse: per block the sub-cells per coarse cell (finest-grid cells over coarse cells) and
  !  their volume. Stops on a level whose cells do not hold a whole number of the finest grid's.
  subroutine Coarse_Level_Setup ( cs, fine, coarse )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    type(coarse_src_t),     intent(inout) :: cs
    type(MOSE_domain_type), intent(in)    :: fine, coarse
    integer :: b, n(3), r(3), lo(3), hi(3), ci, cj, ck

    if (cs%level == coarse%mg_level) return
    if (allocated(cs%ratio)) deallocate(cs%ratio)
    if (allocated(cs%vol))   deallocate(cs%vol)
    allocate(cs%ratio(3, size(coarse%blk)), cs%vol(size(coarse%blk)))
    do b = 1, size(coarse%blk)
      n = coarse%blk(b)%dim
      if (any(mod(fine%blk(b)%dim, n) /= 0)) then
        write(*,'(A,I0,A,I0,A,3(1X,I0),A,3(1X,I0),A)') ' [ERROR] MI2: block ', b, ' of MOSE''s grid level ', &
          coarse%mg_level, ' has', n, ' cells against', fine%blk(b)%dim, &
          ' on the finest grid: not a whole number of sub-cells per cell'
        error stop 1
      endif
      r = fine%blk(b)%dim / n
      cs%ratio(:,b) = r
      allocate(cs%vol(b)%v(n(1), n(2), n(3)))
      do ck = 1, n(3)
      do cj = 1, n(2)
      do ci = 1, n(1)
        lo = ([ci, cj, ck] - 1)*r + 1
        hi = [ci, cj, ck]*r
        cs%vol(b)%v(ci,cj,ck) = sum(fine%blk(b)%vol(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
      enddo; enddo; enddo
    enddo
    cs%level = coarse%mg_level

  end subroutine Coarse_Level_Setup


  !> Stop when an array of a solver does not hold the cells of block b of MOSE's finest grid (fine): ext, its cell extents.
  subroutine Check_Cells ( what, ext, fine, b )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    implicit none
    character(len=*),       intent(in) :: what
    integer,                intent(in) :: ext(3), b
    type(MOSE_domain_type), intent(in) :: fine

    if (all(ext == fine%blk(b)%dim)) return
    write(*,'(A,I0,A,3(1X,I0),A,3(1X,I0))') ' [ERROR] MI2: '//what//' of block ', b, ' has', ext, &
      ' cells, MOSE''s finest grid', fine%blk(b)%dim
    error stop 1

  end subroutine Check_Cells


  !> IGLOO's sources restricted to MOSE's level coarse: rates summed over the sub-cells, Eulerian fields weighted as
  !  IGLOO's normalizeEUL defines them (density, number by volume; velocity, temperature by mass). Outside parallel regions.
  subroutine Restrict_IGLOO_Sources ( cs, fine, coarse, IGLOO, report )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use IGLOO_module,          only: obj_IGLOO
    use IGLOO_variables,       only: toll
    implicit none
    type(coarse_src_t),     intent(inout) :: cs
    type(MOSE_domain_type), intent(in)    :: fine, coarse
    type(obj_IGLOO),        intent(in)    :: IGLOO
    logical,                intent(in)    :: report
    character(len=6), parameter :: EUL_NAME(6) = [character(len=6) :: 'mass', 'number', 'rhou', 'rhov', 'rhow', 'rhoT']
    integer  :: b, m, c, n(3), r(3), i, j, k, ci, cj, ck
    real(R8) :: v, sv, sm, sn, smu(3), smT, ft(6), ct(6), at(6)
    character(len=24) :: nm

    call Coarse_Level_Setup( cs, fine, coarse )
    if (allocated(cs%src)) deallocate(cs%src)
    if (allocated(cs%eul)) deallocate(cs%eul)
    allocate(cs%src(size(coarse%blk)))
    if (allocated(IGLOO%euler)) allocate(cs%eul(size(coarse%blk), size(IGLOO%euler, 2)))
    if (report) write(*,'(A,I0)') ' [MI2-report] restricted IGLOO to level ', cs%level
    do b = 1, size(coarse%blk)
      n = coarse%blk(b)%dim
      r = cs%ratio(:,b)
      ! the source branch: sums
      if (allocated(IGLOO%source(b)%sourceMass)) then
        call Check_Cells( 'IGLOO''s sourceMass', shape3(IGLOO%source(b)%sourceMass), fine, b )
        allocate(cs%src(b)%sourceMass(size(IGLOO%source(b)%sourceMass, 1), n(1), n(2), n(3)))
        call Sum_Sub_Cells( IGLOO%source(b)%sourceMass, cs%src(b)%sourceMass )
        do c = 1, size(IGLOO%source(b)%sourceMass, 1)
          write(nm,'(A,I0,A)') 'igloo-sourceMass(', c, ')'
          call Tell( nm, sum(IGLOO%source(b)%sourceMass(c,:,:,:)), sum(cs%src(b)%sourceMass(c,:,:,:)), &
                     sum(abs(IGLOO%source(b)%sourceMass(c,:,:,:))) )
        enddo
      endif
      if (allocated(IGLOO%source(b)%sourceMom)) then
        call Check_Cells( 'IGLOO''s sourceMom', shape3(IGLOO%source(b)%sourceMom), fine, b )
        allocate(cs%src(b)%sourceMom(3, n(1), n(2), n(3)))
        call Sum_Sub_Cells( IGLOO%source(b)%sourceMom, cs%src(b)%sourceMom )
        do c = 1, 3
          write(nm,'(A,I0,A)') 'igloo-sourceMom(', c, ')'
          call Tell( nm, sum(IGLOO%source(b)%sourceMom(c,:,:,:)), sum(cs%src(b)%sourceMom(c,:,:,:)), &
                     sum(abs(IGLOO%source(b)%sourceMom(c,:,:,:))) )
        enddo
      endif
      if (allocated(IGLOO%source(b)%sourceEn)) then
        call Check_Cells( 'IGLOO''s sourceEn', shape(IGLOO%source(b)%sourceEn), fine, b )
        allocate(cs%src(b)%sourceEn(n(1), n(2), n(3)))
        do ck = 1, n(3)
        do cj = 1, n(2)
        do ci = 1, n(1)
          cs%src(b)%sourceEn(ci,cj,ck) = sum(IGLOO%source(b)%sourceEn((ci-1)*r(1)+1:ci*r(1), (cj-1)*r(2)+1:cj*r(2), &
                                                                      (ck-1)*r(3)+1:ck*r(3)))
        enddo; enddo; enddo
        call Tell( 'igloo-sourceEn', sum(IGLOO%source(b)%sourceEn), sum(cs%src(b)%sourceEn), &
                   sum(abs(IGLOO%source(b)%sourceEn)) )
      endif
      ! the Eulerian fields: density and number by volume, velocity and temperature by mass
      if (.not. allocated(cs%eul)) cycle
      do m = 1, size(cs%eul, 2)
        if (.not. allocated(IGLOO%euler(b,m)%density)) cycle
        associate (fe => IGLOO%euler(b,m))
        call Check_Cells( 'IGLOO''s Eulerian density', shape(fe%density), fine, b )
        allocate(cs%eul(b,m)%density(n(1), n(2), n(3)), cs%eul(b,m)%np(n(1), n(2), n(3)), &
                 cs%eul(b,m)%velocity(3, n(1), n(2), n(3)), cs%eul(b,m)%temperature(n(1), n(2), n(3)))
        ! the totals of the finest grid, cell by cell
        ft = 0._R8; ct = 0._R8; at = 0._R8
        do k = 1, fine%blk(b)%dim(3)
        do j = 1, fine%blk(b)%dim(2)
        do i = 1, fine%blk(b)%dim(1)
          v = fine%blk(b)%vol(i,j,k)
          ft = ft + [fe%density(i,j,k), fe%np(i,j,k), fe%density(i,j,k)*fe%velocity(:,i,j,k), &
                     fe%density(i,j,k)*fe%temperature(i,j,k)]*v
          at = at + abs([fe%density(i,j,k), fe%np(i,j,k), fe%density(i,j,k)*fe%velocity(:,i,j,k), &
                         fe%density(i,j,k)*fe%temperature(i,j,k)])*v
        enddo; enddo; enddo
        do ck = 1, n(3)
        do cj = 1, n(2)
        do ci = 1, n(1)
          sm = 0._R8; sn = 0._R8; smu = 0._R8; smT = 0._R8
          do k = (ck-1)*r(3)+1, ck*r(3)
          do j = (cj-1)*r(2)+1, cj*r(2)
          do i = (ci-1)*r(1)+1, ci*r(1)
            v   = fine%blk(b)%vol(i,j,k)
            sm  = sm  + fe%density(i,j,k)*v
            sn  = sn  + fe%np(i,j,k)*v
            smu = smu + fe%density(i,j,k)*fe%velocity(:,i,j,k)*v
            smT = smT + fe%density(i,j,k)*fe%temperature(i,j,k)*v
          enddo; enddo; enddo
          sv = cs%vol(b)%v(ci,cj,ck)
          cs%eul(b,m)%density(ci,cj,ck) = sm/sv
          cs%eul(b,m)%np(ci,cj,ck)      = sn/sv
          if (cs%eul(b,m)%density(ci,cj,ck) > toll) then
            cs%eul(b,m)%velocity(:,ci,cj,ck)  = smu/sm
            cs%eul(b,m)%temperature(ci,cj,ck) = smT/sm
          else
            cs%eul(b,m)%velocity(:,ci,cj,ck)  = 0._R8
            cs%eul(b,m)%temperature(ci,cj,ck) = 0._R8
          endif
          ct = ct + [cs%eul(b,m)%density(ci,cj,ck), cs%eul(b,m)%np(ci,cj,ck), &
                     cs%eul(b,m)%density(ci,cj,ck)*cs%eul(b,m)%velocity(:,ci,cj,ck), &
                     cs%eul(b,m)%density(ci,cj,ck)*cs%eul(b,m)%temperature(ci,cj,ck)]*sv
        enddo; enddo; enddo
        end associate
        do c = 1, 6
          write(nm,'(A,I0,A,A)') 'igloo-euler(', m, ')-', trim(EUL_NAME(c))
          call Tell( nm, ft(c), ct(c), at(c) )
        enddo
      enddo
    enddo

  contains

    !> y, the sum over the sub-cells of each coarse cell of block b of x (component first)
    subroutine Sum_Sub_Cells ( x, y )
      real(R8), intent(in)  :: x(:,:,:,:)
      real(R8), intent(out) :: y(:,:,:,:)
      y = 0._R8
      do ck = 1, n(3)
      do cj = 1, n(2)
      do ci = 1, n(1)
        do k = (ck-1)*r(3)+1, ck*r(3)
        do j = (cj-1)*r(2)+1, cj*r(2)
        do i = (ci-1)*r(1)+1, ci*r(1)
          y(:,ci,cj,ck) = y(:,ci,cj,ck) + x(:,i,j,k)
        enddo; enddo; enddo
      enddo; enddo; enddo
    end subroutine Sum_Sub_Cells

    subroutine Tell ( what, f, c, a )
      character(len=*), intent(in) :: what
      real(R8),         intent(in) :: f, c, a
      if (report) write(*,'(A,I0,1X,A,3ES24.16)') ' [MI2-report] restricted block ', b, what, f, c, a
    end subroutine Tell

  end subroutine Restrict_IGLOO_Sources


  !> ICE's per-volume rates restricted to MOSE's level coarse: the sub-cells' rate times volume summed, per volume of the
  !  coarse cell, which ICE_Residual multiplies back. Outside any parallel region.
  subroutine Restrict_ICE_Sources ( cs, fine, coarse, ICE_domain, report )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Global_m,          only: ngroups
    implicit none
    type(coarse_src_t),     intent(inout) :: cs
    type(MOSE_domain_type), intent(in)    :: fine, coarse
    type(ICE_domain_type),  intent(in)    :: ICE_domain
    logical,                intent(in)    :: report
    integer  :: b, g, n(3), r(3), i, j, k, ci, cj, ck
    character(len=24) :: nm

    call Coarse_Level_Setup( cs, fine, coarse )
    if (allocated(cs%ice%blk)) deallocate(cs%ice%blk)
    allocate(cs%ice%blk(size(coarse%blk)))
    if (allocated(cs%gr)) deallocate(cs%gr)
    if (allocated(goneRate)) allocate(cs%gr(size(coarse%blk)))
    if (report) write(*,'(A,I0)') ' [MI2-report] restricted ICE to level ', cs%level
    do b = 1, size(coarse%blk)
      n = coarse%blk(b)%dim
      r = cs%ratio(:,b)
      allocate(cs%ice%blk(b)%cond_phase(ngroups))
      do g = 1, ngroups
        associate (fc => ICE_domain%blk(b)%cond_phase(g))
        if (allocated(fc%source)) then
          call Check_Cells( 'ICE''s source', shape3(fc%source), fine, b )
          allocate(cs%ice%blk(b)%cond_phase(g)%source(size(fc%source, 1), n(1), n(2), n(3)))
          call Per_Coarse_Volume( fc%source, cs%ice%blk(b)%cond_phase(g)%source, 'ice-source' )
        endif
        if (allocated(fc%birth)) then
          call Check_Cells( 'ICE''s birth', shape3(fc%birth), fine, b )
          allocate(cs%ice%blk(b)%cond_phase(g)%birth(size(fc%birth, 1), n(1), n(2), n(3)))
          call Per_Coarse_Volume( fc%birth, cs%ice%blk(b)%cond_phase(g)%birth, 'ice-birth' )
        endif
        if (allocated(fc%gone)) then
          call Check_Cells( 'ICE''s gone', shape3(fc%gone), fine, b )
          allocate(cs%ice%blk(b)%cond_phase(g)%gone(size(fc%gone, 1), n(1), n(2), n(3)))
          call Per_Coarse_Volume( fc%gone, cs%ice%blk(b)%cond_phase(g)%gone, 'ice-gone' )
        endif
        end associate
      enddo
      if (allocated(cs%gr)) then
        allocate(cs%gr(b)%a(5, ngroups, n(1), n(2), n(3)))
        do g = 1, ngroups
          call Per_Coarse_Volume( goneRate(b)%a(:,g,:,:,:), cs%gr(b)%a(:,g,:,:,:), 'ice-goneRate' )
        enddo
      endif
    enddo

  contains

    !> y, per coarse cell of block b, the sum over its sub-cells of x times their volume, over the coarse cell's volume;
    !  one report line per component
    subroutine Per_Coarse_Volume ( x, y, what )
      real(R8),         intent(in)  :: x(:,:,:,:)
      real(R8),         intent(out) :: y(:,:,:,:)
      character(len=*), intent(in)  :: what
      real(R8) :: ft(size(x,1)), ct(size(x,1)), at(size(x,1))
      integer  :: c
      ! the totals of the finest grid, cell by cell
      ft = 0._R8; ct = 0._R8; at = 0._R8
      do k = 1, fine%blk(b)%dim(3)
      do j = 1, fine%blk(b)%dim(2)
      do i = 1, fine%blk(b)%dim(1)
        ft = ft + x(:,i,j,k)*fine%blk(b)%vol(i,j,k)
        at = at + abs(x(:,i,j,k))*fine%blk(b)%vol(i,j,k)
      enddo; enddo; enddo
      do ck = 1, n(3)
      do cj = 1, n(2)
      do ci = 1, n(1)
        y(:,ci,cj,ck) = 0._R8
        do k = (ck-1)*r(3)+1, ck*r(3)
        do j = (cj-1)*r(2)+1, cj*r(2)
        do i = (ci-1)*r(1)+1, ci*r(1)
          y(:,ci,cj,ck) = y(:,ci,cj,ck) + x(:,i,j,k)*fine%blk(b)%vol(i,j,k)
        enddo; enddo; enddo
        y(:,ci,cj,ck) = y(:,ci,cj,ck)/coarse%blk(b)%vol(ci,cj,ck)
        ct = ct + y(:,ci,cj,ck)*coarse%blk(b)%vol(ci,cj,ck)
      enddo; enddo; enddo
      if (.not. report) return
      do c = 1, size(x,1)
        write(nm,'(A,I0,A,I0,A)') what//'(', g, ',', c, ')'
        write(*,'(A,I0,1X,A,3ES24.16)') ' [MI2-report] restricted block ', b, nm, ft(c), ct(c), at(c)
      enddo
    end subroutine Per_Coarse_Volume

  end subroutine Restrict_ICE_Sources


  !> The last three extents of a rank-4 array: its cells.
  pure function shape3 ( x ) result( e )
    real(R8), intent(in) :: x(:,:,:,:)
    integer :: e(3)
    e = [size(x,2), size(x,3), size(x,4)]
  end function shape3


  !> What ICE's last coupling step hands the gas and leaves for its next step, kept for a restart: the source of its last
  !  step's stages (or block), and the relaxation times and smoke rates its next step starts from, are not functions of
  !  the state ICE writes, so its field file cannot give them back.
  function ICE_Exchange_File () result(f)
    use ICE_Global_m, only: ICE_phase_prefix
    implicit none
    character(len=:), allocatable :: f

    f = 'OUTPUT/'//trim(ICE_phase_prefix)//'coupling.tec'

  end function ICE_Exchange_File


  !> The variables of ICE_Exchange_File, family by family: source, heat, birth, the removals' store, the relaxation time
  !  that bounds the next step; srate where a smoke family holds one (srate(g)), goneRate where the coupler keeps it.
  subroutine Exchange_Vars ( srate, rate, xv )
    use ICE_Global_m, only: ngroups, ncond
    implicit none
    logical,                   intent(in)  :: srate(:), rate
    type(xvar_t), allocatable, intent(out) :: xv(:)
    integer :: g, n

    allocate(xv(0))
    do g = 1, ngroups
      do n = 1, ncond(g)
        call Add( 'S', XS, n )
      enddo
      call Add( 'H', XH, 0 )
      do n = 1, 5
        call Add( 'B', XB, n )
      enddo
      do n = 1, 5
        call Add( 'G', XG, n )
      enddo
      call Add( 'TAU', XU, 0 )
      if (srate(g)) then
        do n = 1, ncond(g)
          call Add( 'R', XR, n )
        enddo
      endif
      if (rate) then
        do n = 1, 5
          call Add( 'GA', XA, n )
        enddo
        call Add( 'GT', XT, 0 )
      endif
    enddo

  contains

    subroutine Add ( tag, what, n )
      character(len=*), intent(in) :: tag
      integer,          intent(in) :: what, n
      character(len=32) :: nm

      if (n > 0) then
        write(nm,'(A,I0,A,I0)') tag, n, '_p', g
      else
        write(nm,'(A,A,I0)') tag, '_p', g
      endif
      xv = [xv, xvar_t(nm, what, g, n)]

    end subroutine Add

  end subroutine Exchange_Vars


  !> The value of the exchange variable v in cell (i,j,k) of block b; a smoke family without srate there holds none, an
  !  infinite relaxation time (no drag) is written -1.
  function Exchange_Value ( ICE_domain, b, v, i, j, k ) result(x)
    use ICE_Advanced_Types_m, only: ICE_domain_type
    implicit none
    type(ICE_domain_type), intent(in) :: ICE_domain
    integer,               intent(in) :: b, i, j, k
    type(xvar_t),          intent(in) :: v
    real(R8) :: x

    x = 0._R8
    associate (cp => ICE_domain%blk(b)%cond_phase(v%g))
      select case (v%what)
      case (XS); x = cp%source(v%n,i,j,k)
      case (XH); x = cp%heat(i,j,k)
      case (XB); x = cp%birth(v%n,i,j,k)
      case (XG); x = cp%gone(v%n,i,j,k)
      case (XU); x = cp%tau(i,j,k)
                 if (x > huge(x)) x = -1._R8
      case (XR); if (allocated(cp%srate)) x = cp%srate(v%n,i,j,k)
      case (XA); x = goneRate(b)%a(v%n,v%g,i,j,k)
      case (XT); x = goneRate(b)%t(v%g,i,j,k)
      end select
    end associate

  end function Exchange_Value


  !> At the end of a coupling step that left both solvers' fields on disk, and at the end of a run with ICE: the arrays
  !  ICExMOSE and Gone_To_Rate read and those ICE's next step starts from, to ICE_Exchange_File, with ICE's iteration
  !  (CYCLE), whether ICE has produced sources and MOSE's iteration. Tecplot ASCII through ORION: 17 digits, read back
  !  exactly. quiet: no line in the log (the writes of the coupling steps). Written to <file>.part, then renamed: a run
  !  killed during the write keeps the hand-off it had.
  subroutine Write_ICE_Exchange ( ICE_sim, hasSources, iterMOSE, quiet )
    use, intrinsic :: iso_c_binding, only: c_null_char
    use ICE_Advanced_Types_m, only: ICE_simulation_type
    use ICE_Global_m,         only: ngroups
    use ICE_Config_Types_m,   only: ice_sim_param => obj_sim_param
    use Lib_ORION_data,       only: orion_data, copyORION
    use Lib_Tecplot,          only: tec_write_structured_multiblock
    implicit none
    type(ICE_simulation_type), intent(in) :: ICE_sim
    logical,                   intent(in) :: hasSources
    integer,                   intent(in) :: iterMOSE
    logical, optional,         intent(in) :: quiet
    type(orion_data)              :: io
    type(xvar_t),     allocatable :: xv(:)
    character(len=:), allocatable :: names, f
    logical :: srate(ngroups)
    integer :: b, g, v, i, j, k, err

    do g = 1, ngroups
      srate(g) = any([(allocated(ICE_sim%domain(1)%blk(b)%cond_phase(g)%srate), b = 1, ICE_sim%domain(1)%nb)])
    enddo
    call Exchange_Vars( srate, allocated(goneRate), xv )
    call copyORION( ICE_sim%ODP(1), io )
    io%tec%format = 'ascii'; io%tec%node = .false.; io%tec%bc = .false.
    names = ''
    do v = 1, size(xv)
      names = names//' "'//trim(xv(v)%name)//'"'
    enddo
    do b = 1, size(io%block)
      if (allocated(io%block(b)%vars)) deallocate(io%block(b)%vars)
      allocate(io%block(b)%vars(size(xv), io%block(b)%Ni, io%block(b)%Nj, io%block(b)%Nk))
      do k = 1, io%block(b)%Nk
      do j = 1, io%block(b)%Nj
      do i = 1, io%block(b)%Ni
        do v = 1, size(xv)
          io%block(b)%vars(v,i,j,k) = Exchange_Value( ICE_sim%domain(1), b, xv(v), i, j, k )
        enddo
      enddo; enddo; enddo
    enddo
    f = ICE_Exchange_File()
    err = tec_write_structured_multiblock( orion=io, varnames=names, filename=f//'.part', &
                                           cycle=int(ice_sim_param%iter_general, I4), fldnames=XAUX, &
                                           fldvalues=[merge(1._R8, 0._R8, hasSources), real(iterMOSE, R8)] )
    if (err /= 0) then
      write(*,'(A,I0,A)') ' [ERROR] MI2: cannot write '//f//'.part (error ', err, ' of the writer)'
      error stop 1
    endif
    err = c_rename( f//'.part'//c_null_char, f//c_null_char )
    if (err /= 0) then
      write(*,'(A,I0,A)') ' [ERROR] MI2: cannot rename '//f//'.part to '//f//' (error ', err, ' of rename)'
      error stop 1
    endif
    if (present(quiet)) then
      if (quiet) return
    endif
    write(*,'(A)') ' [MI2] what ICE''s last coupling step hands the gas, for a restart: '//f

  end subroutine Write_ICE_Exchange


  !> An ICE restart whose solvers resumed their counters: the arrays of Write_ICE_Exchange read back from
  !  ICE_Exchange_File (used) when it was written at ICE's restored iteration and, when MOSE restarts too (checkMOSE), at
  !  MOSE's (iterMOSE). Missing, incomplete (a file MI2 did not write in full), of another layout or mesh, or written at
  !  another iteration (the run stopped between its outputs), it is not used and why says so: the caller decides from
  !  the schedule.
  subroutine Read_ICE_Exchange ( ICE_sim, iterMOSE, checkMOSE, used, hasSources, why )
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use ICE_Advanced_Types_m, only: ICE_simulation_type
    use ICE_Global_m,         only: ngroups, ncond
    use ICE_Config_Types_m,   only: ice_sim_param => obj_sim_param
    use Lib_ORION_data,       only: orion_data
    use Lib_Tecplot,          only: tec_read_structured_multiblock
    implicit none
    type(ICE_simulation_type), intent(inout) :: ICE_sim
    integer,                   intent(in)    :: iterMOSE
    logical,                   intent(in)    :: checkMOSE
    logical,                   intent(out)   :: used
    logical,                   intent(inout) :: hasSources
    character(len=:), allocatable, intent(out) :: why
    type(orion_data)              :: io
    type(xvar_t),     allocatable :: xv(:)
    character(len=:), allocatable :: f
    integer,          allocatable :: col(:)
    character(len=32) :: nm
    character(len=96) :: at
    integer(I4) :: cyc
    real(R8)    :: aux(2)
    logical     :: found(2), there, srate(ngroups), rate
    integer     :: b, g, v, i, j, k, err, nv, ndir

    used = .false.
    why  = ''
    f = ICE_Exchange_File()
    inquire(file=f, exist=there)
    if (.not. there) then
      why = f//' is missing (MI2 writes it with the fields at the end of each coupling step that leaves both solvers'' '// &
        'fields on disk and at the end of a run with ICE; an older binary did not)'
      return
    endif
    err = tec_read_structured_multiblock( orion=io, filename=f, cycle=cyc, fldnames=XAUX, fldvalues=aux, fldfound=found )
    if (err /= 0 .or. .not. all(found)) then
      write(at,'(A,I0,A)') ' (read error ', err, ', or no records '
      why = f//' is not a file MI2 wrote in full'//trim(at)//' '//trim(XAUX(1))//' and '//trim(XAUX(2))//'): a run '// &
        'killed while it wrote it, or another file'
      return
    endif
    if (cyc /= ice_sim_param%iter_general .or. (checkMOSE .and. nint(aux(2)) /= iterMOSE)) then
      write(at,'(A,I0,A,I0,A,I0)') ' is of ICE iteration ', cyc, ' and MOSE iteration ', nint(aux(2)), &
        ', the fields restored of ICE iteration ', ice_sim_param%iter_general
      why = f//trim(at)
      if (checkMOSE) then
        write(at,'(A,I0)') ' and MOSE iteration ', iterMOSE
      else
        at = ' (MOSE starts anew)'
      endif
      why = why//trim(at)//': the run that wrote them stopped between two of its outputs'
      return
    endif
    if (size(io%block) /= ICE_sim%domain(1)%nb) then
      why = Refusal( 'its number of blocks' )
      return
    endif
    do b = 1, size(io%block)
      if (any([io%block(b)%Ni, io%block(b)%Nj, io%block(b)%Nk] /= ICE_sim%domain(1)%blk(b)%dim(1:3))) then
        why = Refusal( 'the dimensions of a block' )
        return
      endif
      ! the mesh: written from ICE's ORION data and read back exactly
      if (.not. (allocated(io%block(b)%mesh) .and. allocated(ICE_sim%ODP(1)%block(b)%mesh))) then
        why = Refusal( 'the mesh of a block' )
        return
      endif
      if (any(shape(io%block(b)%mesh) /= shape(ICE_sim%ODP(1)%block(b)%mesh))) then
        why = Refusal( 'the mesh of a block' )
        return
      endif
      if (any(io%block(b)%mesh /= ICE_sim%ODP(1)%block(b)%mesh)) then
        why = Refusal( 'the mesh of a block' )
        return
      endif
    enddo
    nv   = size(io%block(1)%vars, 1)
    ndir = size(io%varnames) - nv
    rate = Column('GA1_p1') > 0
    do g = 1, ngroups
      write(nm,'(A,I0)') 'R1_p', g
      srate(g) = Column(nm) > 0
    enddo
    call Exchange_Vars( srate, rate, xv )
    if (size(xv) /= nv) then
      why = Refusal( 'its variables' )
      return
    endif
    allocate(col(size(xv)))
    do v = 1, size(xv)
      col(v) = Column(xv(v)%name)
      if (col(v) == 0) then
        why = Refusal( 'its variables ('//trim(xv(v)%name)//' missing)' )
        return
      endif
    enddo

    if (rate) call Gone_Rate_Setup( ICE_sim%domain(1) )
    do b = 1, ICE_sim%domain(1)%nb
      do g = 1, ngroups
        if (srate(g) .and. .not. allocated(ICE_sim%domain(1)%blk(b)%cond_phase(g)%srate)) &
          allocate(ICE_sim%domain(1)%blk(b)%cond_phase(g)%srate(ncond(g), ICE_sim%domain(1)%blk(b)%dim(1), &
                   ICE_sim%domain(1)%blk(b)%dim(2), ICE_sim%domain(1)%blk(b)%dim(3)))
      enddo
      do k = 1, ICE_sim%domain(1)%blk(b)%dim(3)
      do j = 1, ICE_sim%domain(1)%blk(b)%dim(2)
      do i = 1, ICE_sim%domain(1)%blk(b)%dim(1)
        do v = 1, size(xv)
          associate (cp => ICE_sim%domain(1)%blk(b)%cond_phase(xv(v)%g), x => io%block(b)%vars(col(v),i,j,k), &
                     n => xv(v)%n, g_ => xv(v)%g)
            select case (xv(v)%what)
            case (XS); cp%source(n,i,j,k)       = x
            case (XH); cp%heat(i,j,k)           = x
            case (XB); cp%birth(n,i,j,k)        = x
            case (XG); cp%gone(n,i,j,k)         = x
            case (XU); cp%tau(i,j,k)            = merge(ieee_value(1._R8, ieee_positive_inf), x, x == -1._R8)
            case (XR); cp%srate(n,i,j,k)        = x
            case (XA); goneRate(b)%a(n,g_,i,j,k) = x
            case (XT); goneRate(b)%t(g_,i,j,k)   = x
            end select
          end associate
        enddo
      enddo; enddo; enddo
    enddo
    hasSources = aux(1) > 0.5_R8
    used       = .true.

  contains

    !> The column of variable name among the file's variables, 0 if absent.
    integer function Column ( name )
      character(len=*), intent(in) :: name
      integer :: c

      Column = 0
      do c = 1, nv
        if (trim(io%varnames(ndir+c)) == trim(name)) then
          Column = c
          return
        endif
      enddo

    end function Column

    function Refusal ( what ) result( w )
      character(len=*), intent(in)  :: what
      character(len=:), allocatable :: w

      w = f//' does not match the ICE set-up of this run: '//what

    end function Refusal

  end subroutine Read_ICE_Exchange


end module hydraMI2_Mod_Coupling
