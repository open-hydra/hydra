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
                                     has_stoichiometry, has_stoichiometry_ICE, Build_Datum_Offsets, T_datum_ref, &
                                     Smoke_Record
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
                                     mdotHanded, nHanded, goneShare, Gone_Delivered, Gone_Stored, Gone_Clock, &
                                     Gone_To_Rate, Source_Clock, Source_Mean, &
                                     smokeIG, smokeIGOn, IGLOO_Smoke_Accumulate, IGLOO_Smoke_Deliver, nSmokeNoState, &
                                     handRates, ICE_Exchange_File, Write_ICE_Exchange, Read_ICE_Exchange
  use hydraMI2_Build_Stamp,    only: Print_Build_Stamp
  use Lib_ORION_data,          only: orion_data

  use MOSE_Advanced_Types_m, only: MOSE_simulation_type, MOSE_domain_type
  use MOSE_Config_Types_m,   only: MOSE_obj_sim_param => obj_sim_param, MOSE_obj_io => obj_io, &
                                   MOSE_obj_time_scheme => obj_time_scheme, MOSE_obj_multigrid => obj_multigrid
  use MOSE_Procedures_m,     only: MOSE_type

  use ICE_Advanced_Types_m,  only: ICE_simulation_type
  use ICE_Config_Types_m,    only: ICE_obj_sim_param => obj_sim_param, ICE_obj_time_scheme => obj_time_scheme, &
                                   ICE_obj_condensed => obj_condensed, ICE_obj_io => obj_io, &
                                   ICE_obj_multigrid => obj_multigrid
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
  real(R8) :: t_ice0   = 0._R8     ! ICE's time at the start of its block
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
  !> The solver made the iteration the schedule gives as its last of the run (its postprocess is the end-of-run output)
  logical :: MOSE_ended = .false., ICE_ended = .false.
  logical :: lastIter
  integer :: itermaxKept
  !> The iteration of each solver's field on disk (-1: none of this run or of the run it resumes) and its output count:
  !  a coupling step that leaves both fields at their solver's iteration writes ICE's hand-off with them
  integer :: MOSE_fieldIter = -1, ICE_fieldIter = -1, MOSE_outNum = 0, ICE_outNum = 0
  !> Each solver's iteration after its set-up: a restarted solver that makes no iteration writes no end-of-run output
  integer :: MOSE_iter0 = 0, ICE_iter0 = 0
  !> A solver that restarts and resumes the counters of the run its field comes from (not on a change of kind, steady
  !  <-> time-accurate, nor from a field without them): this run continues that one
  logical :: MOSE_continues = .false., ICE_continues = .false.
  logical :: eulCoupling       = .false.   ! true when this cycle uses the Eulerian coupling mode
  !> Coupling mode in force when the sources MOSE is about to consume were produced.
  logical :: IGLOO_src_eulerian = .false.
  integer :: nAuditSkipped      = 0        ! coupling steps whose sources were Eulerian, so not audited
  integer :: nRunMOSE           = 0        ! MOSE iterations of the block that ran: the share of ICE's removals delivered
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
  ! MOSE's grid levels: only the multigrid cycle is coupled. MOSE then steps its finest grid, the one the coupler hands
  ! the condensed solvers, from the first iteration; a start on the coarse grids is refused
  if ( MOSE_obj_multigrid%MGL > 1 ) then
    if ( trim(MOSE_obj_multigrid%cycle) == 'none' ) then
      write(*,'(A,I0,A)') ' [ERROR] MI2: [MOSE-Multigrid] levels = ', MOSE_obj_multigrid%MGL, ' without a multigrid '// &
        'cycle: MOSE would start on its coarse grids, while the coupler hands IGLOO and ICE the gas of MOSE''s finest '// &
        'grid and adds their sources cell by cell of that grid. Set cycle = v or w (steady runs), or levels = 1'
      error stop 1
    endif
    if ( MOSE_obj_multigrid%iter_threshold(1) <= 0 ) then
      write(*,'(A)') ' [ERROR] MI2: [MOSE-Multigrid] level1-iter is not set: with levels > 1 MOSE takes the iteration '// &
        'limit of its finest grid from it, and with its default (0) the gas makes one iteration per coupling step. '// &
        'Set level1-iter (the value iter-threshold would have)'
      error stop 1
    endif
    MOSE_obj_multigrid%MG_level = 1
    write(*,'(A,I0,A)') ' MI2: MOSE''s multigrid cycle on ', MOSE_obj_multigrid%MGL, ' levels: MOSE steps its finest '// &
      'grid from the first iteration; a coarse level takes IGLOO''s and ICE''s sources summed over its sub-cells, in '// &
      'the residual MOSE restricts to it'
  endif

  ! ICE and IGLOO ingest the gas field as ORION data; each is set up only when enabled
  if (hydra_MI2_obj%run_ICE)   call ICE%setup( ICE_simulation, MOSE_simulation%IOfield(1) )
  ! ICE's grid levels: the coupler hands the gas to ICE's finest level and takes ICE's source from it
  if ( hydra_MI2_obj%run_ICE ) then
    if ( ICE_obj_multigrid%MGL > 1 ) then
      write(*,'(A,I0,A)') ' [ERROR] MI2: [ICE-Multigrid] levels = ', ICE_obj_multigrid%MGL, ': the coupler hands the '// &
        'gas to ICE''s finest grid and takes ICE''s source from it, while ICE steps a coarse level on that level''s own '// &
        'grid, which the gas never reaches. Use levels = 1'
      error stop 1
    endif
  endif
  MOSE_continues = .not. MOSE_obj_sim_param%newrun .and. MOSE_obj_io%restart%found .and. &
                   MOSE_obj_sim_param%iter_general == MOSE_obj_io%restart%cycle
  if ( hydra_MI2_obj%run_ICE ) ICE_continues = .not. ICE_obj_sim_param%newrun .and. ICE_obj_io%restart%found .and. &
                                               ICE_obj_sim_param%iter_general == ICE_obj_io%restart%cycle
  ! a restart with IGLOO is a new IGLOO start: stopped when it would not continue a time-accurate cycle
  if ( hydra_MI2_obj%run_IGLOO ) call Check_IGLOO_Restart()
  ! an ICE restart takes back what ICE's last coupling step handed the gas and left for its next, not in its field file
  if ( hydra_MI2_obj%run_ICE ) then
    if ( .not. ICE_obj_sim_param%newrun ) then
      block
        logical :: used
        character(len=:), allocatable :: why
        ! the hand-off is of a coupling step of the run the fields come from only when both resume its counters
        used = .false.
        if ( .not. ICE_continues .or. .not. (MOSE_obj_sim_param%newrun .or. MOSE_continues) ) then
          why = trim(merge('ICE ', 'MOSE', .not. ICE_continues))//' restarts without the counters of the run it '// &
            'continues (a run of the other kind, steady or time-accurate, or a field without them)'
        else
          call Read_ICE_Exchange( ICE_simulation, MOSE_obj_sim_param%iter_general, .not. MOSE_obj_sim_param%newrun, &
                                  used, ICE_has_sources, why )
          if ( used ) write(*,'(A,I0,A,I0,A)') ' [MI2] ICE restarts with what its last coupling step handed the '// &
            'gas, from '//ICE_Exchange_File()//' (ICE iteration ', ICE_obj_sim_param%iter_general, ', MOSE iteration ', &
            MOSE_obj_sim_param%iter_general, ')'
        endif
        ! without the hand-off the first MOSE block runs without ICE's source: stopped when that ends a continuation
        if ( .not. used ) call ICE_Handoff_Unused( why )
      end block
    endif
  endif
  ! the fields a restart read are on disk at their solver's iteration (with the counters of the run they come from)
  if ( .not. MOSE_obj_sim_param%newrun .and. MOSE_obj_io%restart%found .and. &
       MOSE_obj_sim_param%iter_general == MOSE_obj_io%restart%cycle ) MOSE_fieldIter = MOSE_obj_sim_param%iter_general
  if ( hydra_MI2_obj%run_ICE ) then
    if ( .not. ICE_obj_sim_param%newrun .and. ICE_obj_io%restart%found .and. &
         ICE_obj_sim_param%iter_general == ICE_obj_io%restart%cycle ) ICE_fieldIter = ICE_obj_sim_param%iter_general
  endif
  MOSE_iter0  = MOSE_obj_sim_param%iter_general
  MOSE_outNum = MOSE_obj_io%output_number
  if ( hydra_MI2_obj%run_ICE ) then
    ICE_iter0  = ICE_obj_sim_param%iter_general
    ICE_outNum = ICE_obj_io%output_number
  endif
  ! what the coupler hands between ICE and the gas or IGLOO (the removals, IGLOO's smoke): rates when MOSE or ICE steps
  ! in pseudo-time, so that a converged steady field holds per unit time; amounts when both are time-accurate
  handRates = hydra_MI2_obj%run_ICE .and. &
              .not. (MOSE_obj_time_scheme%time_accurate .and. ICE_obj_time_scheme%time_accurate)
  ! material -> gas species of ICE's phase, from the file ATLAS writes: a material that converts needs its record
  if (hydra_MI2_obj%run_ICE) then
    block
      use ICE_Global_m,       only: ICE_phase_prefix
      use ICE_Lib_Properties, only: mat_e, mat_cp_const, mat_lv
      use ICE_Config_Types_m, only: COMB_BECKSTEAD
      character(len=64), allocatable :: matnames(:)
      logical,           allocatable :: evap(:), burn(:), beck(:)
      real(R8),          allocatable :: eICEref(:), Lvmat(:), qmat(:), Thot(:), eICEhot(:), beta(:)
      logical,           allocatable :: Lvconst(:)
      integer :: im, nm
      nm = size(ICE_obj_condensed)
      allocate(matnames(nm), evap(nm), burn(nm), eICEref(nm), Lvmat(nm), qmat(nm), Thot(nm), eICEhot(nm), Lvconst(nm))
      allocate(beck(nm), beta(nm))
      do im = 1, nm
        matnames(im) = ICE_obj_condensed(im)%name
        evap(im)     = ICE_obj_condensed(im)%evapSelect > 0
        burn(im)     = ICE_obj_condensed(im)%combSelect > 0
        beck(im)     = ICE_obj_condensed(im)%combSelect == COMB_BECKSTEAD
        beta(im)     = merge(ICE_obj_condensed(im)%betaPart, 0._R8, beck(im))
        ! where the material converts: its boiling temperature, or the one it burns from (Beckstead: T-ign)
        Thot(im) = ICE_obj_condensed(im)%Tboil
        if ( beck(im) ) Thot(im) = ICE_obj_condensed(im)%Tign
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
      call Setup_ICE_Smoke( matnames, nm )
      call Report_Stoichiometry( MAP_ICE, matnames, nm )
      call Build_Datum_Offsets( MAP_ICE, matnames, nm, beck, beta, eICEref, Lvmat, qmat, Thot, eICEhot, Lvconst )
      deallocate(matnames, evap, burn, beck, beta, eICEref, Lvmat, qmat, Thot, eICEhot, Lvconst)
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
      ! the smoke of IGLOO's burning materials is born in ICE families
      call Setup_IGLOO_Smoke( matnames, size(IGLOO%material) )
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
        use IGLOO_Lib_Evaporation, only: IGLOO_COMB_BECKSTEAD => COMB_BECKSTEAD
        real(R8), allocatable :: hIGLOOref(:), Lvmat(:), qmat(:), Thot(:), hIGLOOhot(:), beta(:)
        logical,  allocatable :: beck(:)
        logical,  allocatable :: Lvconst(:)
        integer :: it
        allocate(hIGLOOref(size(IGLOO%material)), Lvmat(size(IGLOO%material)), qmat(size(IGLOO%material)), &
                 Thot(size(IGLOO%material)), hIGLOOhot(size(IGLOO%material)), Lvconst(size(IGLOO%material)))
        allocate(beck(size(IGLOO%material)), beta(size(IGLOO%material)))
        do im = 1, size(IGLOO%material)
          Lvmat(im)   = IGLOO%material(im)%Lv
          Lvconst(im) = .not. IGLOO%material(im)%LvVariable
          qmat(im)    = IGLOO%material(im)%qComb
          ! where the material converts: its boiling temperature, or the one it burns from (Beckstead: T-ign)
          Thot(im) = IGLOO%material(im)%Tboil
          beck(im) = IGLOO%material(im)%combSelect == IGLOO_COMB_BECKSTEAD
          if ( beck(im) ) Thot(im) = IGLOO%material(im)%Tign
          beta(im) = merge(IGLOO%material(im)%betaPart, 0._R8, beck(im))
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
        call Build_Datum_Offsets( MAP_IGLOO, matnames, size(IGLOO%material), beck, beta, hIGLOOref, Lvmat, qmat, &
                                  Thot, hIGLOOhot, Lvconst )
        deallocate(hIGLOOref, Lvmat, qmat, Thot, hIGLOOhot, Lvconst)
        deallocate(beck, beta)
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

  call Report_Chimera()

  call Check_Time_Accuracy()

  call Check_Burning_Defaults()

  ! a new ICE run, past every set-up refusal: a hand-off left in OUTPUT/ by an earlier run is not of the fields this
  ! run will write
  if ( hydra_MI2_obj%run_ICE ) then
    if ( ICE_obj_sim_param%newrun ) then
      block
        logical :: there
        integer :: u
        inquire( file=ICE_Exchange_File(), exist=there )
        if ( there ) then
          open( newunit=u, file=ICE_Exchange_File(), status='old' )
          close( u, status='delete' )
          write(*,'(A)') ' [MI2] ICE starts anew: '//ICE_Exchange_File()//' of an earlier run removed'
        endif
      end block
    endif
  endif

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

      ! the removals ICE stores reach the gas over this block: a steady run, as the rate of the store; a time-accurate
      ! one, each MOSE iteration delivers 1/MOSE-iter of them
      nRunMOSE = 0
      if ( handRates ) then
        block
          real(R8) :: Mg, Pg(3), Eg
          call Gone_To_Rate( MOSE_simulation%domain(1), ICE_simulation%domain(1), Mg, Pg, Eg )
          if ( hydra_MI2_obj%report .and. has_stoichiometry_ICE ) &
            write(*,'(A,5ES24.16)') ' [MI2-report] ice-removed-rate', Mg, Pg, Eg
        end block
      elseif ( hydra_MI2_obj%run_ICE .and. ICE_has_sources .and. hydra_obj_sim_param%MOSE_iter(n_cycle) > 0 ) then
        goneShare = 1._R8 / real(hydra_obj_sim_param%MOSE_iter(n_cycle), R8)
      endif
      if ( MOSE_ended .and. hydra_obj_sim_param%MOSE_iter(n_cycle) > 0 ) call Resumed( 'MOSE', MOSE_ended )
      do i_MOSE = 1, hydra_obj_sim_param%MOSE_iter(n_cycle)
        if ( MOSE_obj_sim_param%TODO <= 2 ) then
          ! the run's last MOSE iteration ends on MOSE's own iteration limit (TODO = 3), as a run of MOSE alone: its
          ! postprocess is the end-of-run output
          lastIter = Last_Of_Run( hydra_obj_sim_param%MOSE_iter, i_MOSE )
          if ( lastIter ) then
            itermaxKept = MOSE_simulation%domain(1)%itermax
            MOSE_simulation%domain(1)%itermax = MOSE_obj_sim_param%iter_from_call + 1
          endif
          call MOSE%solve( MOSE_simulation, MOSE_external_routine )
          if ( lastIter ) then
            MOSE_simulation%domain(1)%itermax = itermaxKept
            MOSE_ended = .true.
          endif
          nRunMOSE = nRunMOSE + 1
          if ( smokeIGOn .and. hydra_MI2_obj%run_IGLOO .and. IGLOO_has_sources .and. .not. IGLOO_src_eulerian ) &
            call IGLOO_Smoke_Accumulate( MOSE_simulation%domain(1), IGLOO )
          if ( MOSE_obj_sim_param%TODO <= 2 ) call MOSE%postprocess( MOSE_simulation )
          if ( MOSE_obj_io%output_number /= MOSE_outNum ) then
            MOSE_outNum    = MOSE_obj_io%output_number
            MOSE_fieldIter = MOSE_obj_sim_param%iter_general
          endif
        endif
      enddo
      if ( goneShare > 0._R8 ) then
        block
          real(R8) :: Mg, Pg(3), Eg
          call Gone_Delivered( MOSE_simulation%domain(1), ICE_simulation%domain(1), real(nRunMOSE, R8)*goneShare, &
                               Mg, Pg, Eg )
          if ( hydra_MI2_obj%report ) write(*,'(A,5ES24.16)') ' [MI2-report] ice-removed-delivered', Mg, Pg, Eg
        end block
        goneShare = 0._R8
      endif
      ! the smoke IGLOO's burning made over this block: out of the gas's books, and to ICE in the hand-off's form
      if ( smokeIGOn ) then
        block
          real(R8) :: Ms, Ps(3), Es, Rs(6)
          call IGLOO_Smoke_Deliver( MOSE_simulation%domain(1), ICE_simulation%domain(1), IGLOO, &
                                    IGLOO_has_sources .and. .not. IGLOO_src_eulerian, Ms, Ps, Es, Rs )
          if ( hydra_MI2_obj%report ) then
            if ( handRates ) then
              write(*,'(A,5ES24.16)') ' [MI2-report] igloo-smoke-withheld', Ms, Ps, Es
              write(*,'(A,6ES24.16)') ' [MI2-report] igloo-smoke-rate', Rs
            else
              write(*,'(A,5ES24.16)') ' [MI2-report] igloo-smoke-delivered', Ms, Ps, Es
            endif
          endif
          if ( nSmokeNoState > 0 ) write(*,'(A,I0,A)') ' [WARNING] MI2: in ', nSmokeNoState, ' cells of this block''s '// &
            'iterations a smoke was born with no parcel''s state in IGLOO''s columns: at rest, at T_datum_ref'
          nSmokeNoState = 0
        end block
      endif

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

        t_ice0 = ICE_simulation%domain(1)%time
        if ( ICE_ended .and. hydra_MI2_obj%ICE_iter(n_cycle) > 0 ) call Resumed( 'ICE', ICE_ended )
        do i_ICE = 1, hydra_MI2_obj%ICE_iter(n_cycle)
          if ( ICE_obj_sim_param%TODO <= 2 ) then
            ! the run's last ICE iteration ends on ICE's own iteration limit, as MOSE's
            lastIter = Last_Of_Run( hydra_MI2_obj%ICE_iter, i_ICE )
            if ( lastIter ) then
              itermaxKept = ICE_simulation%domain(1)%itermax
              ICE_simulation%domain(1)%itermax = ICE_obj_sim_param%iter_from_call + 1
            endif
            call ICE%solve( ICE_simulation, Dummy_Function )
            if ( lastIter ) then
              ICE_simulation%domain(1)%itermax = itermaxKept
              ICE_ended = .true.
            endif
            ! a steady run: the time the removals of this iteration were made over
            if ( handRates ) call Gone_Clock( ICE_simulation%domain(1) )
            ! a time-accurate ICE: the gas takes the block's rates, the steps' weighted by their dt
            if ( ICE_obj_time_scheme%time_accurate ) call Source_Clock( ICE_simulation%domain(1) )
            if ( ICE_obj_sim_param%TODO <= 2 ) call ICE%postprocess( ICE_simulation )
            if ( ICE_obj_io%output_number /= ICE_outNum ) then
              ICE_outNum    = ICE_obj_io%output_number
              ICE_fieldIter = ICE_obj_sim_param%iter_general
            endif
            ! after a solve with twcoupled set, cond_phase%source holds the rates of ICE's step
            ICE_has_sources = .true.
            ICE_src_cycle   = n_cycle
            ICE_src_step    = n_couple
          endif
        enddo

        if ( ICE_obj_time_scheme%time_accurate ) call Source_Mean( ICE_simulation%domain(1) )
        ! the ICE time the smoke's rates were integrated over
        if ( handRates .and. hydra_MI2_obj%report .and. hydra_MI2_obj%ICE_iter(n_cycle) /= 0 ) &
          write(*,'(A,2ES24.16)') ' [MI2-report] ice-time', t_ice0, ICE_simulation%domain(1)%time
        ! refuse a mass source of a material without a record
        if ( ICE_has_sources ) call Check_ICE_Mass_Channel( ICE_simulation%domain(1) )
        ! refuse a non-finite source before MOSE applies it
        if ( ICE_has_sources ) call Check_ICE_Source_Finite( MOSE_simulation%domain(1), ICE_simulation%domain(1) )
        ! both fields on disk are of this step's end: ICE's hand-off with them, a run stopped later restarts from here.
        ! With sol-overwrite = false a restart reads the unnumbered fields of the last end of a run, which the hand-off
        ! of that end matches: no hand-off of the numbered outputs in between
        if ( MOSE_fieldIter == MOSE_obj_sim_param%iter_general .and. ICE_fieldIter == ICE_obj_sim_param%iter_general &
             .and. MOSE_obj_io%sol_overwrite .and. ICE_obj_io%sol_overwrite ) &
          call Write_ICE_Exchange( ICE_simulation, ICE_has_sources, MOSE_obj_sim_param%iter_general, quiet=.true. )
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
  ! each solver's end-of-run output (TODO = 3): its last iteration's rows and the solution due there, written once, then
  ! its final solution
  write(*,*)
  write(*,'(A)') ' HYDRA-MI2: FINAL POSTPROCESS'
  ! a restarted solver that made no iteration: its end-of-run output is the one it restarted from, not written again
  MOSE_obj_sim_param%TODO = 3
  if ( .not. MOSE_obj_sim_param%newrun .and. MOSE_obj_sim_param%iter_general == MOSE_iter0 ) then
    write(*,'(A)') ' [MI2] MOSE made no iteration in this run: its end-of-run output is the one it restarted from'
  else
    call MOSE%postprocess( MOSE_simulation )
  endif
  if (hydra_MI2_obj%run_ICE) then
    ICE_obj_sim_param%TODO = 3
    if ( .not. ICE_obj_sim_param%newrun .and. ICE_obj_sim_param%iter_general == ICE_iter0 ) then
      write(*,'(A)') ' [MI2] ICE made no iteration in this run: its end-of-run output is the one it restarted from'
    else
      call ICE%postprocess( ICE_simulation )
    endif
    call Write_ICE_Exchange( ICE_simulation, ICE_has_sources, MOSE_obj_sim_param%iter_general )
  endif
  ! IGLOO's fields are those of its last solve: none when no cycle advanced it
  if (hydra_MI2_obj%run_IGLOO .and. .not. IGLOO_first_sweep) call IGLOO%writeout()

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

    ! a coarse level of MOSE's multigrid cycle: the sources reach it in the residual MOSE restricts from the finest grid
    if ( MOSE_domain%mg_level /= 1 ) return

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


  !> Whether iteration i of the current coupling step is the solver's last of the run: the step ends its cycle, i is
  !  the cycle's X-iter and no later cycle advances the solver (the schedule as last read).
  logical function Last_Of_Run ( Xiter, i )
    implicit none
    integer, intent(in) :: Xiter(:), i

    Last_Of_Run = n_couple >= hydra_obj_sim_param%N_Coupling(n_cycle) .and. i == Xiter(n_cycle) .and. &
                  all( Xiter(n_cycle+1:N_cycles) <= 0 )

  end function Last_Of_Run


  !> A re-read of input.ini that advances a solver after the iteration the schedule had made its last: the run goes on,
  !  that iteration's rows and due solution were not written.
  subroutine Resumed ( solver, ended )
    implicit none
    character(len=*), intent(in)    :: solver
    logical,          intent(inout) :: ended

    write(*,'(A)') ' [WARNING] MI2: input.ini, read again, advances '//solver//' after the iteration the schedule had '// &
      'made its last: the rows and the solution due at that iteration were not written'
    ended = .false.

  end subroutine Resumed


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
    if ( iceFrozen .and. ICE_src_cycle > 0 ) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      ' [MI2] ICE source FROZEN this cycle (ICE-iter = 0): the gas keeps the sources produced in cycle ', &
      ICE_src_cycle, ', coupling step ', ICE_src_step, '  [cycle ', n_cycle, ', step ', n_couple, &
      ', MOSE iter ', iterNow, ']'
    if ( iceFrozen .and. ICE_src_cycle == 0 ) write(*,'(A,I0,A,I0,A,I0,A)') &
      ' [MI2] ICE source FROZEN this cycle (ICE-iter = 0): the gas keeps the sources read at the restart  [cycle ', &
      n_cycle, ', step ', n_couple, ', MOSE iter ', iterNow, ']'

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
    ! the two solvers number their models alike (COMB_BECKSTEAD = 1, COMB_DIRECT = 2)
    character(len=9), parameter :: ICE_WORD(0:2) = ['none     ', 'Beckstead', 'direct   ']
    character(len=9), parameter :: IGL_WORD(0:2) = ICE_WORD
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


  !> newrun = false with IGLOO: IGLOO keeps neither its parcels nor the sources it hands the gas from one run to the
  !  next, so it starts empty on the restored fields. Stop when a cycle in which IGLOO iterates continues a
  !  time-accurate run (a cold IGLOO start does not); otherwise one warning; nothing when IGLOO never iterates.
  subroutine Check_IGLOO_Restart()
    implicit none
    character(len=:), allocatable :: why
    logical :: moseRe, iceRe
    integer :: c

    moseRe = .not. MOSE_obj_sim_param%newrun
    iceRe  = .false.
    if ( hydra_MI2_obj%run_ICE ) iceRe = .not. ICE_obj_sim_param%newrun
    if ( .not. (moseRe .or. iceRe) .or. .not. any( hydra_MI2_obj%IGLOO_iter /= 0 ) ) return

    do c = 1, size( hydra_MI2_obj%IGLOO_iter )
      if ( hydra_MI2_obj%IGLOO_iter(c) == 0 ) cycle
      why = Cycle_Clock( c, 'IGLOO' )
      if ( len(why) == 0 ) cycle
      write(*,'(A,I0,A)') ' [ERROR] MI2: newrun = false with IGLOO, and IGLOO iterates in cycle ', c, ', a time-'// &
        'accurate cycle: '//why//'. IGLOO keeps no parcels and no sources from one run to the next, so it would start '// &
        'empty on the restored fields, and a cold IGLOO start does not continue a time-accurate run. To go on from '// &
        'these fields, start a new run from them: '//Way_Out()
      error stop 1
    enddo

    write(*,*)
    write(*,'(A)') ' [WARNING] MI2: newrun = false with IGLOO: IGLOO keeps no parcels and no sources from one run '// &
      'to the next, so it starts empty on the restored fields: a new IGLOO start. No time-accurate run is continued here '// &
      '(a steady schedule, or counters not resumed): a steady run converges again, a time-accurate one is a new run '// &
      'from the restored fields'
    write(*,*)

  end subroutine Check_IGLOO_Restart


  !> ICE restarts and its hand-off cannot be used (why): the first MOSE block runs without ICE's source. Stop when a
  !  cycle in which ICE iterates continues a time-accurate run, or when ICE iterates in no cycle (the hand-off would be
  !  the only source of ICE the gas gets: another end state); one warning otherwise.
  subroutine ICE_Handoff_Unused( why )
    implicit none
    character(len=*), intent(in) :: why
    character(len=*), parameter :: E = ' [ERROR] MI2: ICE restarts ([ICE-Parameters] newrun = false) and its hand-off '// &
                                       'cannot be used: '
    character(len=:), allocatable :: clock
    integer :: c

    if ( .not. any( hydra_MI2_obj%ICE_iter > 0 ) ) then
      write(*,'(A)') E//why//'. ICE iterates in no cycle of this run, so the hand-off would be the only source of ICE '// &
        'the gas gets, and the run would end on another state. Let ICE iterate in a cycle, or start a new run from '// &
        'these fields: '//Way_Out()
      error stop 1
    endif
    do c = 1, size( hydra_MI2_obj%ICE_iter )
      if ( hydra_MI2_obj%ICE_iter(c) == 0 ) cycle
      clock = Cycle_Clock( c, 'ICE' )
      if ( len(clock) == 0 ) cycle
      write(*,'(A,I0,A)') E//why//'. ICE iterates in cycle ', c, ', a time-accurate cycle: '//clock//'. Without the '// &
        'hand-off the first MOSE block would run without ICE''s source, and the run would not continue the one it '// &
        'restarts. To go on from these fields, start a new run from them: '//Way_Out()
      error stop 1
    enddo
    write(*,*)
    write(*,'(A)') ' [WARNING] MI2: ICE''s hand-off is not used: '//why//'. The first MOSE block runs without ICE''s '// &
      'source. No time-accurate run is continued here (a steady schedule, or counters not resumed): a steady run '// &
      'converges again, a time-accurate one is a new run from the restored fields'
    write(*,*)

  end subroutine ICE_Handoff_Unused


  !> Why cycle c, in which IGLOO or ICE (who) iterates, continues a time-accurate run ('' if it does not): the gas it
  !  hands its sources to steps in physical time and MOSE resumed its counters (the sources stay in force until the
  !  next solve, whichever cycle MOSE iterates in), or ICE iterates in it in physical time and resumed its counters
  function Cycle_Clock( c, who ) result( why )
    integer,          intent(in) :: c
    character(len=*), intent(in) :: who
    character(len=:), allocatable :: why
    why = ''
    if ( MOSE_continues .and. MOSE_obj_time_scheme%time_accurate .and. any( hydra_obj_sim_param%MOSE_iter > 0 ) ) then
      why = 'the gas '//who//' hands its sources to steps in physical time ([MOSE-Numerics] time-accurate = true) '// &
            'and continues the restored run (MOSE resumes its counters)'
    elseif ( ICE_continues ) then
      if ( hydra_MI2_obj%ICE_iter(c) > 0 .and. ICE_obj_time_scheme%time_accurate ) &
        why = 'ICE iterates in it in physical time ([ICE-Numerics] time-accurate = true) and continues the restored '// &
              'run (ICE resumes its counters)'
    endif
  end function Cycle_Clock


  !> The way out of a restart that cannot continue: a new run from the restored fields, each copied as the initial
  !  condition (keeping the one it replaces), in the format it was written in, newrun = true in the sections of the
  !  solvers that restart. Items separated by '; '
  function Way_Out() result( s )
    use MOSE_Global_m,   only: MOSE_phase_prefix
    use ICE_Global_m,    only: ICE_phase_prefix
    use ICE_IO_Solution, only: ICE_condinit => condinit
    character(len=:), allocatable :: s, keys
    logical :: moseRe, iceRe
    moseRe = .not. MOSE_obj_sim_param%newrun
    iceRe  = .false.
    if ( hydra_MI2_obj%run_ICE ) iceRe = .not. ICE_obj_sim_param%newrun
    s    = ''
    keys = ''
    if ( moseRe ) then
      s    = Field_As_IC( trim(MOSE_obj_io%nameinit), MOSE_phase_prefix, MOSE_obj_io%ini_format, &
                          MOSE_obj_io%sol_format, '[MOSE-IO]' )
      keys = '[MOSE-Parameters]'
    endif
    if ( iceRe ) then
      s    = s//Field_As_IC( trim(ICE_condinit), ICE_phase_prefix, ICE_obj_io%ini_format, ICE_obj_io%sol_format, &
                             '[ICE-IO]' )
      keys = keys//trim(merge(' and ', '     ', moseRe))//' [ICE-Parameters]'
    endif
    s = s//'set newrun = true in '//trim(adjustl(keys))
  end function Way_Out


  !> The items of Way_Out for one restored field: keep a copy of the initial condition it replaces, copy it there (a
  !  .vtm with the block files it lists, OUTPUT/vtk/<name>Block<b>.vts, which it names relative to its directory), and
  !  the ic-format of the field's format when the two differ
  function Field_As_IC( field, prefix, iniFmt, solFmt, ioSection ) result( s )
    character(len=*), intent(in)  :: field, prefix, iniFmt, solFmt, ioSection
    character(len=:), allocatable :: s, ext, ic, base
    logical :: there
    ext  = field(index(field, '.', back=.true.):)
    base = field(index(field, '/', back=.true.)+1:index(field, '.', back=.true.)-1)
    ic   = 'INPUT/'//trim(prefix)//'ic'//ext
    s    = ''
    inquire( file=ic, exist=there )
    if ( there ) s = 'keep a copy of '//ic//' (the initial condition it replaces); '
    s = s//'copy '//field//' to '//ic
    if ( ext == '.vtm' ) s = s//' and OUTPUT/vtk/'//base//'Block*.vts to INPUT/vtk/'
    s = s//'; '
    if ( trim(adjustl(iniFmt)) /= trim(adjustl(solFmt)) ) &
      s = s//'set ic-format = '//trim(adjustl(solFmt))//' in '//ioSection//'; '
  end function Field_As_IC


  !> The smoke of ICE's burning materials (the tail of their combustion records): the smoke material is a material
  !  of ICE's phase with exactly one family, of the MK closure, that neither evaporates nor burns; the burning material
  !  gives d-smoke. Sets ICE's smokeFam and nuSmoke and the smoke's datum offset on ICE's own energy of it.
  subroutine Setup_ICE_Smoke( matnames, nm )
    use ICE_Global_m,           only: ICE_phase_prefix
    use hydraMI2_Stoichiometry, only: Smoke_Record, Set_Smoke_Datum, MAP_ICE
    implicit none
    character(len=*), intent(in) :: matnames(:)
    integer,          intent(in) :: nm
    character(len=64) :: sname, sphase
    character(len=*), parameter :: E = ' [ERROR] MI2: '
    real(R8) :: snu, sh, dh
    integer  :: m, gs

    do m = 1, nm
      if ( .not. Smoke_Record( MAP_ICE, m, sname, snu, sh, sphase ) ) cycle
      call Smoke_Family( trim(matnames(m)), 'INPUT/'//trim(ICE_phase_prefix)//'products.txt', sname, sphase, sh, gs, dh )
      if ( .not. ICE_obj_condensed(m)%dSmoke > 0._R8 ) then
        write(*,'(A)') E//trim(matnames(m))//' burns into the smoke '//trim(sname)//' and its row of INPUT/'// &
          trim(ICE_phase_prefix)//'models.txt gives no d-smoke, the diameter each smoke drop is born with'
        error stop 1
      endif
      ICE_obj_condensed(m)%smokeFam = gs
      ICE_obj_condensed(m)%nuSmoke  = snu
      call Set_Smoke_Datum( m, dh )
      write(*,'(A,I0,A,ES15.7,A,ES12.4,A,ES15.7,A)') ' [MI2] ICE material '//trim(matnames(m))//' burns into the '// &
        'smoke '//trim(sname)//' (family ', gs, '): ', snu, ' kg per kg burnt, d-smoke ', ICE_obj_condensed(m)%dSmoke, &
        ' m, datum offset ', dh, ' J/kg'
    enddo

  end subroutine Setup_ICE_Smoke


  !> The ICE family a smoke is born in, for the ICE path and the IGLOO path alike: the smoke material is a material of
  !  ICE's phase (the record's phase) with exactly one family, of the MK closure, that neither evaporates nor burns.
  !  Returns the family and the smoke's datum offset: its formation enthalpy sh less ICE's energy of it at T_datum_ref.
  subroutine Smoke_Family( owner, file, sname, sphase, sh, gs, dh )
    use ICE_Global_m,           only: ICE_phase_prefix, ICE_mat_of => mat_of, ICE_ngroups => ngroups
    use ICE_Lib_Properties,     only: mat_e, mat_cp_const
    use hydraMI2_Stoichiometry, only: T_datum_ref
    implicit none
    character(len=*), intent(in)  :: owner, file, sname, sphase
    real(R8),         intent(in)  :: sh
    integer,          intent(out) :: gs
    real(R8),         intent(out) :: dh
    character(len=*), parameter :: E = ' [ERROR] MI2: '
    integer  :: g, ms, nfam
    real(R8) :: es

    ms = 0
    do g = 1, size(ICE_obj_condensed)
      if ( trim(ICE_obj_condensed(g)%name) == trim(sname) ) ms = g
    enddo
    if ( ms == 0 ) then
      write(*,'(A)') E//trim(file)//': '//owner//' burns into the '// &
        'smoke '//trim(sname)//', which is not a material of ICE''s phase: declare it in INPUT/'// &
        trim(ICE_phase_prefix)//'phase.txt with one population and its family with closure = MK'
      error stop 1
    endif
    if ( trim(sphase) /= trim(ICE_phase_prefix) ) then
      write(*,'(A)') E//'the smoke of '//owner//' is born in the phase '//trim(sphase)// &
        ', and ICE runs the phase '//trim(ICE_phase_prefix)//': run ATLAS GPB again'
      error stop 1
    endif
    nfam = 0; gs = 0
    do g = 1, ICE_ngroups
      if ( ICE_mat_of(g) == ms ) then
        nfam = nfam + 1; gs = g
      endif
    enddo
    if ( nfam /= 1 ) then
      write(*,'(A,I0,A)') E//'the smoke '//trim(sname)//' of '//owner//' has ', nfam, &
        ' families in ICE: it is born in one, give its material one population'
      error stop 1
    endif
    if ( trim(ICE_obj_time_scheme%model(gs)) /= 'MK' ) then
      write(*,'(A,I0,A)') E//'the smoke '//trim(sname)//' is born in ICE family ', gs, ', whose closure is '// &
        trim(ICE_obj_time_scheme%model(gs))//': the smoke takes mass and number at a fixed diameter, which the '// &
        'MK closure carries; set closure = MK'
      error stop 1
    endif
    if ( ICE_obj_condensed(ms)%evapSelect > 0 .or. ICE_obj_condensed(ms)%combSelect > 0 ) then
      write(*,'(A)') E//'the smoke '//trim(sname)//' evaporates or burns in ICE: a smoke does neither'
      error stop 1
    endif
    if ( ICE_obj_condensed(ms)%cs_varies ) then
      es = mat_e( ICE_obj_condensed(ms), T_datum_ref )
    else
      es = mat_cp_const( ICE_obj_condensed(ms) ) * T_datum_ref
    endif
    dh = sh - es

  end subroutine Smoke_Family


  !> The smoke of IGLOO's burning materials (the tail of a combustion record, and of the release record of a cap oxide
  !  handed at burnout): ICE runs and holds its family (Smoke_Family), the coupling is not relaxed, the material gives
  !  d-smoke. Fills smokeIG; one report line per smoke.
  subroutine Setup_IGLOO_Smoke( matnames, nm )
    use IGLOO_variables,        only: IGLOO_phase_prefix
    use hydraMI2_Stoichiometry, only: Smoke_Record, Release_Smoke_Record, MAP_IGLOO
    implicit none
    character(len=*), intent(in) :: matnames(:)
    integer,          intent(in) :: nm
    character(len=*), parameter :: E = ' [ERROR] MI2: '
    character(len=64) :: sname, sphase
    character(len=256) :: file
    real(R8) :: snu, sh, dh
    integer  :: m, p, gs, n
    logical  :: has

    allocate(smokeIG(max(nm,1)))
    file = 'INPUT/'//trim(IGLOO_phase_prefix)//'products.txt'
    do m = 1, nm
      do n = 1, 2
        if ( n == 1 ) then
          has = Smoke_Record( MAP_IGLOO, m, sname, snu, sh, sphase )
        else
          p = IGLOO%material(m)%prodID
          has = .false.
          if ( IGLOO%material(m)%combSelect > 0 .and. p > 0 .and. .not. IGLOO%material(m)%resOn ) &
            has = Release_Smoke_Record( MAP_IGLOO, p, sname, snu, sh, sphase )
        endif
        if ( .not. has ) cycle
        if ( .not. hydra_MI2_obj%run_ICE ) then
          write(*,'(A)') E//trim(file)//': the IGLOO material '//trim(matnames(m))//' burns into the smoke '// &
            trim(sname)//', and the smoke has no family: ICE is not run. Run ICE with '//trim(sname)// &
            ' in its phase, or burn the material into gas species only'
          error stop 1
        endif
        if ( hydra_MI2_obj%relaxation ) then
          write(*,'(A)') E//'the IGLOO material '//trim(matnames(m))//' burns into the smoke '//trim(sname)// &
            ': relaxation = true blends two coupling steps of the gas''s share, and the smoke''s is not blended; '// &
            'set relaxation = false'
          error stop 1
        endif
        if ( .not. IGLOO%material(m)%dSmoke > 0._R8 ) then
          write(*,'(A)') E//'the IGLOO material '//trim(matnames(m))//' burns into the smoke '//trim(sname)// &
            ' and its row of INPUT/'//trim(IGLOO_phase_prefix)//'models.txt gives no d-smoke, the diameter each '// &
            'smoke drop is born with'
          error stop 1
        endif
        call Smoke_Family( 'the IGLOO material '//trim(matnames(m)), file, sname, sphase, sh, gs, dh )
        smokeIG(m)%fam(n)  = gs
        smokeIG(m)%nu(n)   = snu
        smokeIG(m)%dh(n)   = dh
        smokeIG(m)%dSmoke  = IGLOO%material(m)%dSmoke
        smokeIGOn = .true.
        write(*,'(A,I0,A,ES15.7,A,ES12.4,A,ES15.7,A)') ' [MI2] IGLOO material '//trim(matnames(m))//' '// &
          trim(merge('burns into' // repeat(' ', 17), 'hands its cap at burnout to', n == 1))//' the smoke '//trim(sname)// &
          ' (ICE family ', gs, &
          '): ', snu, ' kg per kg, d-smoke ', smokeIG(m)%dSmoke, ' m, datum offset ', dh, ' J/kg'
      enddo
    enddo
    if ( smokeIGOn ) then
      if ( handRates ) then
        write(*,'(A)') ' [MI2] IGLOO''s smoke reaches ICE as rates: a source of each smoke family over ICE''s own steps'
      else
        write(*,'(A)') ' [MI2] IGLOO''s smoke reaches ICE as amounts: what each MOSE block made, after the block'
      endif
    endif

  end subroutine Setup_IGLOO_Smoke


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
      ! what ICE's removals still hold for the gas, on the gas's scale (mass through the rows, momentum, energy)
      if ( has_stoichiometry_ICE ) then
        call Gone_Stored( MOSE_simulation%domain(1), ICE_simulation%domain(1), M, F, Q )
        write(*,'(A,5ES24.16)') ' [MI2-report] ice-removed-stored', M, F, Q
      endif
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


  !> An overset mesh (gas faces of kind 102): the sums over the domain the coupler prints take every cell of every
  !  block, so they count an overlap once per block that covers it. One warning; nothing on a mesh without such a face.
  subroutine Report_Chimera()
    implicit none
    integer :: i, n, nb
    logical, allocatable :: has(:)

    associate ( dom => MOSE_simulation%domain(1) )
      if ( .not. allocated(dom%bc) ) return
      allocate(has(dom%nb)); has = .false.
      n = 0
      do i = 1, size(dom%bc)
        if ( dom%bc(i)%type /= 102 ) cycle
        n = n + 1
        has(dom%bc(i)%b) = .true.
      enddo
      nb = count(has)
    end associate
    if ( n == 0 ) return

    write(*,*)
    write(*,'(A,I0,A,I0,A)') ' [WARNING] MI2: an overset mesh: ', n, ' chimera (102) boundary cells on ', nb, &
      ' block(s). Every sum over the domain the coupler prints takes every cell of every block, so it counts an '// &
      'overlap once per block that covers it: the gas mass, momentum and energy, what ICE and IGLOO hand the gas, '// &
      'the coupling report and the mass and energy balances at the end of the run. The coupler makes no correction '// &
      'for the overlap'
    write(*,*)

  end subroutine Report_Chimera

end program hydra_MI2_program
