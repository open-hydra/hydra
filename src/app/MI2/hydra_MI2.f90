!>
!> hydra-MI2 -- MOSE (Eulerian gas) + ICE (Eulerian condensed) + IGLOO (Lagrangian) coupler.
!> OpenMP only: the root CMakeLists refuses USE_MPI for MI2.
!>
program hydra_MI2_program
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use hydra_Parameters_m
  use hydra_Config_Types_m
  use hydraMI2_Read_Ini,       only: Read_Inifile_MI2, Read_Inifile_MI2_Runtime
  use hydraMI2_Config_Types_m, only: hydra_MI2_obj
  use hydraMI2_Mod_Import,     only: Assign_Phases
  use hydraMI2_Stoichiometry,  only: Load_Species_Map, Report_Stoichiometry, MAP_IGLOO, MAP_ICE, &
                                     has_stoichiometry, has_stoichiometry_ICE, Build_Datum_Offsets, T_datum_ref
  use hydraMI2_Mod_Coupling,   only: MOSExIGLOO, Gas_Field_For_IGLOO, IGLOOxMOSE, &
                                     MOSExICE, ICExMOSE, &
                                     resField, Initialize_Residual_Store, Build_Family_Material_Map, &
                                     Accumulate_Injected_Mass, massInjected, Gas_Species_Mass, &
                                     nMassLimited, Accumulate_Injected_Energy, Accumulate_Injected_ICE, &
                                     Gas_Total_Energy, &
                                     Gas_Mean_Temperature, energyInjected, energyDatumPart, &
                                     energyAudited, Check_ICE_Mass_Channel, &
                                     Check_ICE_Source_Finite, Check_IGLOO_Source_Finite, &
                                     Check_IGLOO_Euler_Finite, &
                                     Gas_Momentum, ICE_To_Gas_Totals, IGLOO_To_Gas_Totals, &
                                     mdotHanded, nHanded
  use hydraMI2_Build_Stamp,    only: Print_Build_Stamp
  use Lib_ORION_data,          only: orion_data

  use MOSE_Advanced_Types_m, only: MOSE_simulation_type, MOSE_domain_type
  use MOSE_Config_Types_m,   only: MOSE_obj_sim_param => obj_sim_param, MOSE_obj_io => obj_io, &
                                   MOSE_obj_time_scheme => obj_time_scheme
  use MOSE_Procedures_m,     only: MOSE_type

  use ICE_Advanced_Types_m,  only: ICE_simulation_type
  use ICE_Config_Types_m,    only: ICE_obj_sim_param => obj_sim_param, ICE_obj_time_scheme => obj_time_scheme, &
                                   ICE_obj_condensed => obj_condensed
  use ICE_Procedures_m,      only: ICE_type

  use IGLOO_module,          only: obj_IGLOO

  implicit none

  type(MOSE_type)            :: MOSE
  type(MOSE_simulation_type) :: MOSE_simulation
  type(ICE_type)             :: ICE
  type(ICE_simulation_type)  :: ICE_simulation
  type(obj_IGLOO)            :: IGLOO
  ! Local
  integer :: nthreads
  integer :: n_cycle, n_couple, N_cycles
  integer :: TODO, i_MOSE, i_ICE
  logical :: IGLOO_first_sweep = .true.
  real(R8) :: t_block0 = 0._R8
  real(R8), allocatable :: massGas0(:), massGas1(:)
  real(R8) :: enGas0 = 0._R8, enGas1 = 0._R8
  real(R8) :: Tgas0  = 0._R8, Tgas1  = 0._R8, mGas0 = 0._R8, mGas1 = 0._R8
  real(R8) :: cvGas0 = 0._R8, cvGas1 = 0._R8
  type(orion_data) :: gas_for_igloo
  type(resField), allocatable :: oldResidual(:)
  !> Set once a solver has produced sources and never cleared: frozen sources keep acting on the gas.
  logical :: IGLOO_has_sources = .false.  ! true once IGLOO has produced sources MOSE may read
  logical :: ICE_has_sources   = .false.  ! ditto for ICE: cond_phase%source is only valid after a solve
  integer :: IGLOO_src_cycle = 0, IGLOO_src_step = 0   ! (cycle, coupling step) that produced IGLOO's sources
  integer :: ICE_src_cycle   = 0, ICE_src_step   = 0   ! same for ICE
  integer :: frozenPrintIter = -1                      ! cumulative MOSE iteration of the last frozen notice
  logical :: eulCoupling       = .false.   ! true when this cycle uses the Eulerian coupling mode
  !> Coupling mode in force when the sources MOSE is about to consume were produced.
  logical :: IGLOO_src_eulerian = .false.
  integer :: nAuditSkipped      = 0        ! coupling steps whose sources were Eulerian, so not audited
  real(R8) :: momRep0(3), momRep1(3), enRep0 = 0._R8, enRep1 = 0._R8   ! gas totals around a MOSE block

  call Print_Build_Stamp()

#if defined (_OPENMP)
  ! the team size of a parallel region, written by one thread
  !$omp parallel
  !$omp single
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end single
  !$omp end parallel
  write(*,'(A)')    ' Parallel execution'
  write(*,'(A)')    ' OpenMP:'
  write(*,'(A,I4)') ' -  Number of threads --> ', nthreads
#else
  nthreads = 1
  write(*,'(A)')    ' Serial execution'
#endif
  MOSE_obj_sim_param%nthreads = nthreads
  ICE_obj_sim_param%nthreads  = nthreads

  ! Populate hydra_obj_* AND hydra_MI2_obj from input.ini
  call Read_Inifile_MI2 ()

  ! Report the parsed coupling schedule before any solver setup
  call Report_MI2_Config()

  ! Multigrid is not supported: IGLOO has no levels
  if (hydra_obj_multigrid%MGL > 1) then
    write(*,'(A,I0,A)') '[ERROR] hydra-MI2 does not support multigrid (HYDRA-Multigrid levels = ', &
                        hydra_obj_multigrid%MGL, '). IGLOO has no multigrid levels. Use levels = 1.'
    error stop 1
  endif

  ! Route each *phase.txt of INPUT/ to its owning solver (ICE and IGLOO each solve their own phase)
  call Assign_Phases()

  ! ---------------------------------------------------------------------------
  ! Solver setup
  ! ---------------------------------------------------------------------------
  ! ICE must know it is two-way coupled before setup (it then fills cond_phase%source)
  ICE_obj_sim_param%twcoupled = hydra_MI2_obj%run_ICE

  call MOSE%setup( MOSE_simulation )

  ! ICE and IGLOO ingest the gas field as ORION data; each is set up only when enabled
  if (hydra_MI2_obj%run_ICE)   call ICE%setup( ICE_simulation, MOSE_simulation%IOfield(1) )
  ! material -> gas species of ICE's phase, from the file ATLAS writes: a material that converts needs its record
  if (hydra_MI2_obj%run_ICE) then
    block
      use ICE_Global_m,       only: ICE_phase_prefix
      use ICE_Lib_Properties, only: mat_e, mat_cp_const, mat_lv
      character(len=64), allocatable :: matnames(:)
      logical,           allocatable :: evap(:), burn(:), beck(:)
      real(R8),          allocatable :: eICEref(:), Lvmat(:), qmat(:), Thot(:), eICEhot(:)
      logical,           allocatable :: Lvconst(:)
      integer :: im, nm
      nm = size(ICE_obj_condensed)
      allocate(matnames(nm), evap(nm), burn(nm), eICEref(nm), Lvmat(nm), qmat(nm), Thot(nm), eICEhot(nm), Lvconst(nm))
      allocate(beck(nm))
      do im = 1, nm
        matnames(im) = ICE_obj_condensed(im)%name
        evap(im)     = ICE_obj_condensed(im)%evapSelect > 0
        burn(im)     = ICE_obj_condensed(im)%combSelect > 0
        beck(im)     = ICE_obj_condensed(im)%combSelect == 2
        ! where the material converts: its boiling temperature, or the one it burns from (Beckstead: T-ign)
        Thot(im) = ICE_obj_condensed(im)%Tboil
        if ( ICE_obj_condensed(im)%combSelect == 2 ) Thot(im) = ICE_obj_condensed(im)%Tign
        ! the energy ICE credits the material at T_ref and there on its own datum: its table's, or cp*T
        if ( ICE_obj_condensed(im)%cs_varies ) then
          eICEref(im) = mat_e( ICE_obj_condensed(im), T_datum_ref )
          eICEhot(im) = mat_e( ICE_obj_condensed(im), Thot(im) )
        else
          eICEref(im) = mat_cp_const( ICE_obj_condensed(im) ) * T_datum_ref
          eICEhot(im) = mat_cp_const( ICE_obj_condensed(im) ) * Thot(im)
        endif
        Lvmat(im)   = mat_lv( ICE_obj_condensed(im), T_datum_ref )
        Lvconst(im) = .not. ICE_obj_condensed(im)%use_lv
        qmat(im)    = ICE_obj_condensed(im)%qComb
      enddo
      call Load_Species_Map( MAP_ICE, ICE_phase_prefix, matnames, evap, burn )
      call Report_Stoichiometry( MAP_ICE, matnames, nm )
      call Build_Datum_Offsets( MAP_ICE, matnames, nm, beck, eICEref, Lvmat, qmat, Thot, eICEhot, Lvconst )
      deallocate(matnames, evap, burn, beck, eICEref, Lvmat, qmat, Thot, eICEhot, Lvconst)
    end block
  endif
  ! IGLOO%setup includes the reset_state of the first solve; every later solve needs its own
  if (hydra_MI2_obj%run_IGLOO) then
    call Gas_Field_For_IGLOO( MOSE_simulation%domain(1), MOSE_simulation%IOfield(1), gas_for_igloo )
    call IGLOO%setup( gas_for_igloo )
    call Build_Family_Material_Map( IGLOO )
    ! material -> gas species of IGLOO's phase, from the file ATLAS writes: a material that converts needs its record
    block
      use IGLOO_variables,      only: IGLOO_phase_prefix
      use IGLOO_Lib_Combustion, only: nuOx
      character(len=20), allocatable :: matnames(:)
      logical,           allocatable :: evap(:), burn(:), handsCap(:)
      integer,           allocatable :: oxide(:)
      integer :: im
      allocate(matnames(size(IGLOO%material)), evap(size(IGLOO%material)), burn(size(IGLOO%material)), &
               handsCap(size(IGLOO%material)), oxide(size(IGLOO%material)))
      do im = 1, size(IGLOO%material)
        matnames(im) = IGLOO%material(im)%matName
        evap(im)     = IGLOO%material(im)%evapSelect > 0
        burn(im)     = IGLOO%material(im)%combSelect > 0
        ! a material that burns with an oxide cap: the material its oxide is, and whether the cap ends in the gas
        oxide(im)    = 0
        if ( burn(im) ) oxide(im) = IGLOO%material(im)%prodID
        handsCap(im) = oxide(im) > 0 .and. .not. IGLOO%material(im)%resOn
      enddo
      call Load_Species_Map( MAP_IGLOO, IGLOO_phase_prefix, matnames, evap, burn, oxide, handsCap, nuOx )
      ! the Eulerian coupling hands the gas no mass, no species and no heat of reaction
      if ( has_stoichiometry .and. any(hydra_MI2_obj%IGLOO_iter == 2) ) then
        write(*,'(A)') ' [ERROR] MI2: an IGLOO material evaporates or burns and a cycle has IGLOO-iter = 2: the '// &
          'Eulerian coupling hands the gas momentum and energy only, without the mass, its species and its '// &
          'heat of reaction. Use IGLOO-iter = 1 (source coupling) in every cycle that runs IGLOO'
        error stop 1
      endif
      call Report_Stoichiometry( MAP_IGLOO, matnames, size(IGLOO%material) )

      ! enthalpy-datum offsets; hIGLOOref is the enthalpy IGLOO credits for the material at T_ref on
      ! the datum of its table: cp*T_ref + hOff for constant cp, hTab interpolated at T_ref otherwise
      block
        real(R8), allocatable :: hIGLOOref(:), Lvmat(:), qmat(:), Thot(:), hIGLOOhot(:)
        logical,  allocatable :: beck(:)
        logical,  allocatable :: Lvconst(:)
        integer :: it
        allocate(hIGLOOref(size(IGLOO%material)), Lvmat(size(IGLOO%material)), qmat(size(IGLOO%material)), &
                 Thot(size(IGLOO%material)), hIGLOOhot(size(IGLOO%material)), Lvconst(size(IGLOO%material)))
        allocate(beck(size(IGLOO%material)))
        do im = 1, size(IGLOO%material)
          Lvmat(im)   = IGLOO%material(im)%Lv
          Lvconst(im) = .not. IGLOO%material(im)%LvVariable
          qmat(im)    = IGLOO%material(im)%qComb
          ! where the material converts: its boiling temperature, or the one it burns from (Beckstead: T-ign)
          Thot(im) = IGLOO%material(im)%Tboil
          if ( IGLOO%material(im)%combSelect == 1 ) Thot(im) = IGLOO%material(im)%Tign
          beck(im) = IGLOO%material(im)%combSelect == 1
          if ( IGLOO%material(im)%cpVariable ) then
            associate ( ht => IGLOO%material(im)%hTab )
              it = max(lbound(ht,1), min(int(T_datum_ref), ubound(ht,1)-1))
              hIGLOOref(im) = ht(it) + (ht(it+1)-ht(it))*(T_datum_ref - real(it,R8))
              it = max(lbound(ht,1), min(int(Thot(im)), ubound(ht,1)-1))
              hIGLOOhot(im) = ht(it) + (ht(it+1)-ht(it))*(Thot(im) - real(it,R8))
            end associate
          else
            hIGLOOref(im) = IGLOO%material(im)%cp * T_datum_ref + IGLOO%material(im)%hOff
            hIGLOOhot(im) = IGLOO%material(im)%cp * Thot(im) + IGLOO%material(im)%hOff
          endif
          write(*,'(A,I0,A,ES14.7,A)') ' [MI2] IGLOO material ', im, ': table datum '// &
            trim(IGLOO%material(im)%hDatum)//', hOff = ', IGLOO%material(im)%hOff, ' J/kg'
        enddo
        call Build_Datum_Offsets( MAP_IGLOO, matnames, size(IGLOO%material), beck, hIGLOOref, Lvmat, qmat, Thot, &
                                  hIGLOOhot, Lvconst )
        deallocate(hIGLOOref, Lvmat, qmat, Thot, hIGLOOhot, Lvconst)
        deallocate(beck)
      end block

      deallocate(matnames, evap, burn, handsCap, oxide)
    end block
    ! under-relaxation store, allocated whether or not relaxation is on
    call Initialize_Residual_Store( MOSE_simulation%domain(1), oldResidual )
  endif

  ! baseline gas inventory for the closed-domain mass and energy gates reported at the end
  block
    use MOSE_Global_m, only: nsc
    allocate(massGas0(nsc), massGas1(nsc))
    call Gas_Species_Mass  ( MOSE_simulation%domain(1), massGas0 )
    call Gas_Total_Energy  ( MOSE_simulation%domain(1), enGas0 )
    call Gas_Mean_Temperature( MOSE_simulation%domain(1), Tgas0, mGas0, cvGas0 )
  end block

  call Report_Setup()

  call Check_Time_Accuracy()

  call Check_Burning_Defaults()

  ! ---------------------------------------------------------------------------
  ! Coupling cycles
  ! ---------------------------------------------------------------------------
  ! Outer loop over cycles (one entry of the MOSE-iter / ICE-iter / IGLOO-iter lists), inner loop
  ! of N_Coupling(n_cycle) coupling steps. Solve order: MOSE, IGLOO, ICE.
  N_cycles = size(hydra_obj_sim_param%MOSE_iter)
  n_cycle  = 0

  do while ( n_cycle < N_cycles )
    n_cycle = n_cycle + 1
    write(*,*)
    write(*,'(A,I0,A,I0)') ' HYDRA-MI2: CYCLE ', n_cycle, ' / ', N_cycles

    n_couple = 0
    TODO = 1
    MOSE_obj_sim_param%TODO = TODO
    if (hydra_MI2_obj%run_ICE) ICE_obj_sim_param%TODO = TODO

    do while ( TODO <= 2 )

      n_couple = n_couple + 1

      ! Hot re-read of input.ini, so a long run can be steered without restarting
      if ( hydra_obj_io%ini_diter > 0 ) then
        if ( mod(n_couple, hydra_obj_io%ini_diter) == 0 ) call Read_Inifile_MI2_Runtime()
      endif

      ! --- MOSE ---------------------------------------------------------------
      ! announce sources frozen this cycle before the block that applies them
      call Report_Frozen_Sources()

      ! domain clock at block start, for the mass and energy audit
      t_block0 = MOSE_simulation%domain(1)%time
      if ( hydra_MI2_obj%report ) then
        call Gas_Momentum    ( MOSE_simulation%domain(1), momRep0 )
        call Gas_Total_Energy( MOSE_simulation%domain(1), enRep0 )
      endif

      do i_MOSE = 1, hydra_obj_sim_param%MOSE_iter(n_cycle)
        if ( MOSE_obj_sim_param%TODO <= 2 ) then
          call MOSE%solve( MOSE_simulation, MOSE_external_routine )
          if ( MOSE_obj_sim_param%TODO <= 2 ) call MOSE%postprocess( MOSE_simulation )
        endif
      enddo

      ! audit the mass and energy this block injected (source branch only: the Eulerian branch
      ! is mass-inert)
      if ( IGLOO_has_sources .and. .not. IGLOO_src_eulerian ) then
        call Accumulate_Injected_Mass  ( IGLOO, MOSE_simulation%domain(1)%time - t_block0 )
        call Accumulate_Injected_Energy( IGLOO, MOSE_simulation%domain(1)%time - t_block0 )
      elseif ( IGLOO_has_sources ) then
        nAuditSkipped = nAuditSkipped + 1
      endif
      if ( hydra_MI2_obj%run_ICE .and. ICE_has_sources .and. has_stoichiometry_ICE ) then
        call Accumulate_Injected_ICE( MOSE_simulation%domain(1), ICE_simulation%domain(1), &
                                      MOSE_simulation%domain(1)%time - t_block0 )
      endif

      if ( hydra_MI2_obj%report ) call Report_Coupling_Step()

      ! --- IGLOO --------------------------------------------------------------
      if ( hydra_MI2_obj%run_IGLOO ) then
        if ( hydra_MI2_obj%IGLOO_iter(n_cycle) /= 0 ) then
          call Set_IGLOO_Switches(n_cycle)

          ! reset_state is mandatory before every solve except the first (IGLOO%setup ran it);
          ! it drops the sources MOSE has just consumed
          if ( IGLOO_first_sweep ) then
            IGLOO_first_sweep = .false.
          else
            call IGLOO%reset_state()
          endif

          ! gas state in, then solve
          call MOSExIGLOO( MOSE_simulation%domain(1), IGLOO )
          if ( hydra_MI2_obj%report ) write(*,'(A,ES24.16,I10)') ' [MI2-report] igloo-inlet-mdot', &
                                                                 mdotHanded, nHanded

          call IGLOO%solve()

          ! IGLOO%source / %euler are valid from here until the next reset_state
          IGLOO_has_sources  = .true.
          IGLOO_src_eulerian = eulCoupling
          IGLOO_src_cycle    = n_cycle
          IGLOO_src_step     = n_couple

          ! refuse a non-finite value before MOSE applies it: the source-mode rates, or the Eulerian fields
          if ( IGLOO_src_eulerian ) then
            call Check_IGLOO_Euler_Finite( MOSE_simulation%domain(1), IGLOO )
          else
            call Check_IGLOO_Source_Finite( MOSE_simulation%domain(1), IGLOO )
          endif
        endif
      endif

      ! --- ICE ----------------------------------------------------------------
      if ( hydra_MI2_obj%run_ICE ) then
        ! gas state in, then solve (ICE computes its sources from gas_phase internally)
        if ( hydra_MI2_obj%ICE_iter(n_cycle) /= 0 ) then
          call MOSExICE( MOSE_simulation%domain(1), ICE_simulation%domain(1) )
        endif

        do i_ICE = 1, hydra_MI2_obj%ICE_iter(n_cycle)
          if ( ICE_obj_sim_param%TODO <= 2 ) then
            call ICE%solve( ICE_simulation, Dummy_Function )
            if ( ICE_obj_sim_param%TODO <= 2 ) call ICE%postprocess( ICE_simulation )
            ! after a solve with twcoupled set, cond_phase%source holds the rates of ICE's step
            ICE_has_sources = .true.
            ICE_src_cycle   = n_cycle
            ICE_src_step    = n_couple
          endif
        enddo

        ! refuse a mass source of a material without a record
        if ( ICE_has_sources ) call Check_ICE_Mass_Channel( ICE_simulation%domain(1) )
        ! refuse a non-finite source before MOSE applies it
        if ( ICE_has_sources ) call Check_ICE_Source_Finite( MOSE_simulation%domain(1), ICE_simulation%domain(1) )
      endif

      ! --- Exit conditions ----------------------------------------------------
      ! the coupling step count is the only thing that ends a cycle: a solver's TODO = 3 is
      ! reset to 1, not propagated
      if ( n_couple >= hydra_obj_sim_param%N_Coupling(n_cycle) ) TODO = 3

      if ( TODO <= 2 ) then
        MOSE_obj_sim_param%TODO = 1
        if (hydra_MI2_obj%run_ICE) ICE_obj_sim_param%TODO = 1
      endif

    enddo

  enddo

  ! ---------------------------------------------------------------------------
  ! Final postprocess
  ! ---------------------------------------------------------------------------
  write(*,*)
  write(*,'(A)') ' HYDRA-MI2: FINAL POSTPROCESS'
  MOSE_obj_sim_param%TODO = 2
  call MOSE%postprocess( MOSE_simulation )
  if (hydra_MI2_obj%run_ICE) then
    ICE_obj_sim_param%TODO = 2
    call ICE%postprocess( ICE_simulation )
  endif
  if (hydra_MI2_obj%run_IGLOO) call IGLOO%writeout()

  if ( nMassLimited > 0 ) then
    write(*,*)
    write(*,'(A,I0,A)') ' [MI2] WARNING: the consumed-species floor guard clamped ', nMassLimited, &
      ' (cell,species) extractions.'
    write(*,'(A)') '       A gas species ran out where a reaction wanted to consume it, so the'
    write(*,'(A)') '       reaction was MASS-LIMITED: this is not the stoichiometry you asked for.'
  endif

  write(*,*)
  ! Mass-coupling audit: injected mass per species and the closed-domain balance against the
  ! change in gas mass. Meaningful only for a time-accurate run.
  if ( allocated(massInjected) ) then
    write(*,*)
    write(*,'(A)') ' [MI2] mass injected into the gas over the run [kg]'
    if ( nAuditSkipped > 0 ) then
      write(*,'(A,I0,A)') '   [WARNING] ', nAuditSkipped, ' coupling step(s) used the EULERIAN branch, '// &
                          'which is mass-inert, and are excluded.'
    endif
    if ( .not. hydra_obj_sim_param%time_accurate ) then
      write(*,'(A)') '   (time-accurate is OFF -- these are pseudo-time integrals, not physical mass)'
    endif
    block
      integer :: is
      do is = 1, size(massInjected)
        if ( abs(massInjected(is)) > 0._R8 ) &
          write(*,'(A,I0,A,ES14.6)') '   species ', is, ' : ', massInjected(is)
      enddo
      write(*,'(A,ES14.6)') '   TOTAL              : ', sum(massInjected)
    end block

    ! closed-domain balance: injected vs gas gained, both measured at the same instant
    call Gas_Species_Mass( MOSE_simulation%domain(1), massGas1 )
    write(*,*)
    write(*,'(A)') ' [MI2] closed-domain mass balance  (valid only if NO boundary flux)'
    write(*,'(A)') '   species        injected           gas gained         rel.err'
    block
      integer  :: is
      real(R8) :: dgas, den, floor_
      ! ignore species whose change is round-off on the total inventory
      floor_ = 1.e-12_R8 * max(sum(abs(massGas0)), tiny(1._R8))
      do is = 1, size(massGas1)
        dgas = massGas1(is) - massGas0(is)
        den  = max(abs(massInjected(is)), abs(dgas))
        if ( den <= floor_ ) then
          write(*,'(A,I0,A)') '   ', is, '   (untouched; change is round-off on the inventory)'
        else
          write(*,'(A,I0,A,ES16.6,A,ES16.6,A,ES11.3)') '   ', is, '   ', massInjected(is), &
                '   ', dgas, '   ', abs(dgas-massInjected(is))/den
        endif
      enddo
    end block
    write(*,*)
  endif

  ! Energy bookkeeping (datum-blind) and the enthalpy-datum diagnostic
  if ( energyAudited ) then
    call Gas_Total_Energy    ( MOSE_simulation%domain(1), enGas1 )
    call Gas_Mean_Temperature( MOSE_simulation%domain(1), Tgas1, mGas1, cvGas1 )
    write(*,*)
    write(*,'(A)') ' [MI2] closed-domain ENERGY bookkeeping  (valid only if NO boundary flux)'
    if ( .not. hydra_obj_sim_param%time_accurate ) then
      write(*,'(A)') '   (time-accurate is OFF -- pseudo-time integral, not physical energy)'
    endif
    if ( nAuditSkipped > 0 ) then
      write(*,'(A,I0,A)') '   [WARNING] ', nAuditSkipped, ' coupling step(s) used the EULERIAN branch '// &
                          'and are NOT in these totals.'
      write(*,'(A)')      '             The balance below is therefore incomplete -- it is only a gate '// &
                          'for a pure source-mode run.'
    endif
    write(*,'(A,ES16.8)') '   handed to gas  (sourceEn + datum) [J] : ', energyInjected
    write(*,'(A,ES16.8)') '     of which the datum correction   [J] : ', energyDatumPart
    write(*,'(A,ES16.8)') '   gas total energy gained          [J] : ', enGas1 - enGas0
    block
      real(R8) :: den
      den = max(abs(energyInjected), abs(enGas1-enGas0))
      if (den > 0._R8) &
        write(*,'(A,ES11.3)') '   rel.err (bookkeeping only, datum-blind) : ', &
                              abs((enGas1-enGas0) - energyInjected)/den
    end block

    write(*,*)
    write(*,'(A)') ' [MI2] enthalpy-datum diagnostic'
    write(*,'(A,F12.6,A,F12.6,A,ES11.3)') '   mass-weighted gas T:  ', Tgas0, ' -> ', Tgas1, &
                                          '   delta = ', Tgas1-Tgas0
    if ( abs(energyDatumPart) > 0._R8 .and. mGas1*cvGas1 > 0._R8 ) then
      write(*,'(A,F12.4,A,ES12.4)') '   mass-weighted frozen cv [J/kg/K] : ', cvGas1, &
                                    '   gas mass [kg] : ', mGas1
      write(*,'(A,ES16.8)') '   energy the gas would have received WITHOUT it [J] : ', &
                            energyInjected - energyDatumPart
      write(*,'(A,F12.6)')  '   => predicted delta-T without the correction  [K] : ', &
                            (Tgas1-Tgas0) - energyDatumPart/(mGas1*cvGas1)
      write(*,'(A)') '   That prediction holds the evaporated mass FIXED. Re-running with'
      write(*,'(A)') '   [MI2-Coupling] datum-correction = false must flip the sign of delta-T, but'
      write(*,'(A)') '   will not reproduce the value: the gas state feeds back into the evaporation'
      write(*,'(A)') '   rate, so the two runs evaporate different amounts. The falsifiable invariant'
      write(*,'(A)') '   is the per-kg offset printed at setup, not this linearised delta-T.'
      write(*,'(A)') '   Scope: this validates the COUPLER''s datum bookkeeping, not the evaporation'
      write(*,'(A)') '   model -- IGLOO''s latent heat cancels on both sides of the comparison.'
    endif
    write(*,*)
  endif

  write(*,'(A)') ' GAME OVER'
  write(*,*)


contains

  !> Source-term injection point, called by MOSE once per Runge-Kutta stage with the current
  !  domain, inside its parallel region (orphaned !$OMP DO only). Adds IGLOOxMOSE then ICExMOSE,
  !  each gated on its solver being active and having produced sources.
  subroutine MOSE_external_routine(MOSE_domain)
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain

    if ( hydra_MI2_obj%run_IGLOO .and. IGLOO_has_sources ) then
      call IGLOOxMOSE( MOSE_domain, IGLOO, eulCoupling, hydra_MI2_obj%relaxation,          &
                       i_MOSE, hydra_obj_sim_param%MOSE_iter(n_cycle),                     &
                       hydra_MI2_obj%relCoeff(n_cycle), oldResidual )
    endif

    ! MOSE_domain is the callback's argument, not the module-level simulation
    if ( hydra_MI2_obj%run_ICE .and. ICE_has_sources ) then
      call ICExMOSE( MOSE_domain, ICE_simulation%domain(1) )
    endif

  end subroutine MOSE_external_routine


  !> ICE accepts an external routine but never calls it; this satisfies the argument list.
  subroutine Dummy_Function
  end subroutine Dummy_Function


  !> Announce sources frozen this cycle (solver enabled, X-iter(n_cycle) = 0, sources already
  !  produced): at the first coupling step of the cycle, then at [MOSE-IO] shell-diter cadence.
  subroutine Report_Frozen_Sources()
    implicit none
    integer :: iterNow, cadence
    logical :: iglooFrozen, iceFrozen, due

    ! X_iter has size 1 when its solver is off
    iglooFrozen = .false.
    if ( hydra_MI2_obj%run_IGLOO .and. IGLOO_has_sources ) iglooFrozen = ( hydra_MI2_obj%IGLOO_iter(n_cycle) == 0 )
    iceFrozen = .false.
    if ( hydra_MI2_obj%run_ICE .and. ICE_has_sources )     iceFrozen   = ( hydra_MI2_obj%ICE_iter(n_cycle)   == 0 )
    if ( .not. (iglooFrozen .or. iceFrozen) ) return

    iterNow = MOSE_obj_sim_param%iter_general
    cadence = max(MOSE_obj_io%shell_diter, 1)
    due = ( n_couple == 1 ) .or. ( frozenPrintIter < 0 ) .or. ( iterNow/cadence > frozenPrintIter/cadence )
    if ( .not. due ) return
    frozenPrintIter = iterNow

    if ( iglooFrozen ) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      ' [MI2] IGLOO source FROZEN this cycle (IGLOO-iter = 0): the gas keeps the sources produced in cycle ', &
      IGLOO_src_cycle, ', coupling step ', IGLOO_src_step, '  [cycle ', n_cycle, ', step ', n_couple, &
      ', MOSE iter ', iterNow, ']'
    if ( iceFrozen ) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      ' [MI2] ICE source FROZEN this cycle (ICE-iter = 0): the gas keeps the sources produced in cycle ', &
      ICE_src_cycle, ', coupling step ', ICE_src_step, '  [cycle ', n_cycle, ', step ', n_couple, &
      ', MOSE iter ', iterNow, ']'

  end subroutine Report_Frozen_Sources


  !> Select IGLOO's coupling mode for this cycle: IGLOO-iter = 2 -> Eulerian (eulSwitch),
  !  1 -> source (srcSwitch). On the last cycle that runs IGLOO both switches are on for the writeout.
  subroutine Set_IGLOO_Switches(cycle_now)
    implicit none
    integer, intent(in) :: cycle_now
    integer :: i
    logical :: is_last

    if ( hydra_MI2_obj%IGLOO_iter(cycle_now) == 2 ) then
      eulCoupling      = .true.
      IGLOO%eulSwitch  = .true.
      IGLOO%srcSwitch  = .false.
    else
      eulCoupling      = .false.
      IGLOO%eulSwitch  = .false.
      IGLOO%srcSwitch  = .true.
    endif

    ! Is this the last cycle that runs IGLOO at all?
    is_last = .true.
    do i = cycle_now + 1, N_cycles
      if ( hydra_MI2_obj%IGLOO_iter(i) /= 0 ) then
        is_last = .false.
        exit
      endif
    enddo
    if ( is_last ) then
      IGLOO%eulSwitch = .true.
      IGLOO%srcSwitch = .true.
    endif

  end subroutine Set_IGLOO_Switches


  !> Echo the coupling schedule exactly as MI2 resolved it, after clamping and defaulting.
  subroutine Report_MI2_Config()
    implicit none
    integer :: i

    write(*,*)
    write(*,'(A)') ' HYDRA-MI2 coupling schedule'
    write(*,'(A,I0)')  ' -  cycles (from MOSE-iter)   --> ', size(hydra_obj_sim_param%MOSE_iter)
    write(*,'(A,L1)')  ' -  run ICE                   --> ', hydra_MI2_obj%run_ICE
    write(*,'(A,L1)')  ' -  run IGLOO                 --> ', hydra_MI2_obj%run_IGLOO
    write(*,'(A,L1)')  ' -  relaxation                --> ', hydra_MI2_obj%relaxation

    write(*,'(A)',advance='no') ' -  MOSE-iter                 --> '
    do i = 1, size(hydra_obj_sim_param%MOSE_iter)
      write(*,'(I0,1X)',advance='no') hydra_obj_sim_param%MOSE_iter(i)
    enddo
    write(*,*)

    write(*,'(A)',advance='no') ' -  N-coupling                --> '
    do i = 1, size(hydra_obj_sim_param%N_Coupling)
      write(*,'(I0,1X)',advance='no') hydra_obj_sim_param%N_Coupling(i)
    enddo
    write(*,*)

    if (hydra_MI2_obj%run_ICE) then
      write(*,'(A)',advance='no') ' -  ICE-iter                  --> '
      do i = 1, size(hydra_MI2_obj%ICE_iter)
        write(*,'(I0,1X)',advance='no') hydra_MI2_obj%ICE_iter(i)
      enddo
      write(*,*)
    endif

    if (hydra_MI2_obj%run_IGLOO) then
      write(*,'(A)',advance='no') ' -  IGLOO-iter (0/1/2)        --> '
      do i = 1, size(hydra_MI2_obj%IGLOO_iter)
        write(*,'(I0,1X)',advance='no') hydra_MI2_obj%IGLOO_iter(i)
      enddo
      write(*,*)
      write(*,'(A)',advance='no') ' -  relaxation-coeff          --> '
      do i = 1, size(hydra_MI2_obj%relCoeff)
        write(*,'(F6.3,1X)',advance='no') hydra_MI2_obj%relCoeff(i)
      enddo
      write(*,*)
    endif
    write(*,*)

  end subroutine Report_MI2_Config


  !> A material name that ICE and IGLOO both carry and burn differently, where a phase file does not say the model
  !  (no combustion token on the material line: the solver's section decided): a message, nothing stops.
  subroutine Check_Burning_Defaults()
    use ICE_Global_m,    only: ICE_phase_prefix
    use IGLOO_variables, only: IGLOO_phase_prefix
    implicit none
    character(len=9), parameter :: ICE_WORD(0:2) = ['none     ', 'direct   ', 'Beckstead']
    character(len=9), parameter :: IGL_WORD(0:2) = ['none     ', 'Beckstead', 'direct   ']
    integer :: mi, mg
    logical :: saidICE, saidIGLOO

    if ( .not. (hydra_MI2_obj%run_ICE .and. hydra_MI2_obj%run_IGLOO) ) return
    do mi = 1, size(ICE_obj_condensed)
      do mg = 1, size(IGLOO%material)
        if ( trim(adjustl(ICE_obj_condensed(mi)%name)) /= trim(adjustl(IGLOO%material(mg)%matName)) ) cycle
        if ( IGLOO%material(mg)%ngroups == 0 ) cycle
        if ( trim(ICE_WORD(ICE_obj_condensed(mi)%combSelect)) == trim(IGL_WORD(IGLOO%material(mg)%combSelect)) ) cycle
        saidICE   = Line_Says_Model( 'INPUT/'//trim(ICE_phase_prefix)//'phase.txt',   ICE_obj_condensed(mi)%name )
        saidIGLOO = Line_Says_Model( 'INPUT/'//trim(IGLOO_phase_prefix)//'phase.txt', IGLOO%material(mg)%matName )
        if ( saidICE .and. saidIGLOO ) cycle
        write(*,'(A)') ' [WARNING] MI2: material '//trim(adjustl(ICE_obj_condensed(mi)%name))//' burns as '// &
          trim(ICE_WORD(ICE_obj_condensed(mi)%combSelect))//' in ICE and as '// &
          trim(IGL_WORD(IGLOO%material(mg)%combSelect))//' in IGLOO, and '// &
          trim(merge('neither phase file says    ', 'one phase file does not say', .not. (saidICE .or. saidIGLOO)))// &
          ' its model: there the solver''s section decided ([ICE-Physics] combustion, [IGLOO-Models] combustion). '// &
          'Write combustion in both [GPB-Phase*] sections (none for a material that does not burn) to say it is meant'
      enddo
    enddo

  end subroutine Check_Burning_Defaults


  !> True when the line of the material in a phase file carries a combustion token.
  logical function Line_Says_Model( file, material ) result(says)
    implicit none
    character(len=*), intent(in) :: file, material
    character(len=1024) :: line
    character(len=64)   :: first
    integer :: u, ios
    logical :: exists

    says = .false.
    inquire(file=file, exist=exists)
    if ( .not. exists ) return
    open(newunit=u, file=file, status='old', action='read')
    read(u,'(A)',iostat=ios) line                 ! the type word
    do
      read(u,'(A)',iostat=ios) line
      if ( ios /= 0 ) exit
      if ( len_trim(line) == 0 ) cycle
      read(line,*,iostat=ios) first
      if ( ios /= 0 ) cycle
      if ( trim(first) /= trim(adjustl(material)) ) cycle
      says = index(line, ' combustion=') > 0
      exit
    enddo
    close(u)

  end function Line_Says_Model


  !> Refuse a cycle that advances MOSE and ICE together while either is time-accurate: the two
  !  solvers do not share a time step, so their physical times would part.
  subroutine Check_Time_Accuracy()
    implicit none
    integer :: c

    if ( .not. hydra_MI2_obj%run_ICE ) return
    if ( .not. (MOSE_obj_time_scheme%time_accurate .or. ICE_obj_time_scheme%time_accurate) ) return
    do c = 1, size(hydra_obj_sim_param%MOSE_iter)
      if ( hydra_obj_sim_param%MOSE_iter(c) > 0 .and. hydra_MI2_obj%ICE_iter(c) > 0 ) then
        write(*,'(A,I0,A)') '[ERROR] [MI2] cycle ', c, ' advances MOSE and ICE together with time-accurate on:'
        write(*,'(A)')      '        they do not share a time step, so their physical times would part. Advance'
        write(*,'(A)')      '        them in separate cycles, or set time-accurate = false in [MOSE-Numerics] and [ICE-Numerics].'
        error stop 1
      endif
    enddo

  end subroutine Check_Time_Accuracy


  !> [MI2-Coupling] report: the MOSE block of this coupling step and what the coupler handed the gas
  !  during it (the sources in force, before the solves of this step replace them).
  subroutine Report_Coupling_Step()
    implicit none
    real(R8) :: F(3), Q, M, taumin, taumax

    call Gas_Momentum    ( MOSE_simulation%domain(1), momRep1 )
    call Gas_Total_Energy( MOSE_simulation%domain(1), enRep1 )
    write(*,'(A,I0,A,I0,2(A,ES24.16))') ' [MI2-report] cycle ', n_cycle, ' step ', n_couple, &
                                         ' t0 ', t_block0, ' t1 ', MOSE_simulation%domain(1)%time
    write(*,'(A,6ES24.16)') ' [MI2-report] gas-momentum', momRep0, momRep1
    write(*,'(A,2ES24.16)') ' [MI2-report] gas-energy', enRep0, enRep1
    if ( hydra_MI2_obj%run_ICE .and. ICE_has_sources ) then
      call ICE_To_Gas_Totals( MOSE_simulation%domain(1), ICE_simulation%domain(1), F, Q, M, taumin, taumax )
      write(*,'(A,6ES24.16)') ' [MI2-report] ice-to-gas', F, Q, taumin, taumax
      if ( has_stoichiometry_ICE ) write(*,'(A,ES24.16)') ' [MI2-report] ice-to-gas-mass', M
    endif
    if ( hydra_MI2_obj%run_IGLOO .and. IGLOO_has_sources .and. .not. IGLOO_src_eulerian ) then
      call IGLOO_To_Gas_Totals( MOSE_simulation%domain(1), IGLOO, F, Q, M )
      write(*,'(A,5ES24.16)') ' [MI2-report] igloo-to-gas', F, Q, M
    endif

  end subroutine Report_Coupling_Step


  !> What actually got built, after setup.
  subroutine Report_Setup()
    implicit none

    write(*,*)
    write(*,'(A)') ' HYDRA-MI2 solvers'
    write(*,'(A,I0)') ' -  MOSE  blocks              --> ', MOSE_simulation%domain(1)%nb
    if (hydra_MI2_obj%run_ICE) then
      write(*,'(A,I0)') ' -  ICE   blocks              --> ', ICE_simulation%domain(1)%nb
    endif
    if (hydra_MI2_obj%run_IGLOO) then
      write(*,'(A,I0)') ' -  IGLOO geometry blocks     --> ', size(IGLOO%geoblock)
      write(*,'(A,I0)') ' -  IGLOO materials           --> ', size(IGLOO%material)
    endif
    write(*,*)

  end subroutine Report_Setup

end program hydra_MI2_program
