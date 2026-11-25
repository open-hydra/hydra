!> hydra-MI2's input: input.ini read (Read_Inifile_MI2), re-read every [HYDRA-IO] ini-diter coupling steps
!  (Read_Inifile_MI2_Runtime), and what MI2 refuses or warns about in the solvers' own sections once their setup has
!  read them (Check_MOSE_Grid_Levels, Check_ICE_Grid_Levels, Check_Coupled_Inputs).
module hydraMI2_Read_Ini
  use iso_fortran_env, only: R8 => real64

  implicit none
  private
  public :: Read_Inifile_MI2, Read_Inifile_MI2_Runtime
  public :: Check_MOSE_Grid_Levels, Check_ICE_Grid_Levels, Check_Coupled_Inputs, Check_MOSE_Level1_Iter

contains

  !> MI2's input reader: the common registrations plus MI2's own, populating hydra_obj_* and
  !  hydra_MI2_obj from input.ini. Echoes the coupling schedule as resolved and refuses what MI2
  !  cannot run.
  subroutine Read_Inifile_MI2()
    use Finer,                only: file_ini
    use hydra_Read_Sim_Param, only: Register_Sim_Param
    use hydra_Read_IO,        only: Register_IO_Fields
    use hydra_Read_Numerics,  only: Register_Numerics
    use hydra_Backend_INI,    only: Load_Ini, Scan_Ini
    use hydra_Config_Types_m, only: hydra_obj_multigrid
    use hydraMI2_Read_Params, only: Register_MI2_Params
    implicit none
    type(file_ini) :: fini
    integer :: nmgl

    call fini%load(filename='input.ini')

    ! Size the common arrays and read the multigrid depth
    call Scan_Ini(fini, nmgl)
    ! Size the MI2 arrays before registration
    call Scan_Ini_MI2(fini)

    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Numerics(nmgl)
    call Register_MI2_Params()
    call Release_Unused_Keys()

    call Load_Ini(fini)

    call Normalize_MI2_Params()

    call Report_Schedule()

    ! Multigrid is not supported: IGLOO has no levels
    if (hydra_obj_multigrid%MGL > 1) then
      write(*,'(A,I0,A)') '[ERROR] hydra-MI2 does not support multigrid (HYDRA-Multigrid levels = ', &
                          hydra_obj_multigrid%MGL, '). IGLOO has no multigrid levels. Use levels = 1.'
      error stop 1
    endif

  end subroutine Read_Inifile_MI2


  !> Hot re-read of input.ini every [HYDRA-IO] ini-diter coupling steps (0: never), so a long run
  !  can be steered without restarting: refreshes the values of the registered keys (sizes are not
  !  re-scanned).
  subroutine Read_Inifile_MI2_Runtime(n_couple)
    use Finer,                   only: file_ini
    use hydra_Backend_INI,       only: Load_Ini
    use hydra_Config_Types_m,    only: hydra_obj_io
    use hydra_Input_Registry,    only: Validate_Registry
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    integer, intent(in) :: n_couple      !> coupling step of the cycle
    type(file_ini) :: fini
    character(len=1024) :: out

    if (hydra_obj_io%ini_diter <= 0) return
    if (mod(n_couple, hydra_obj_io%ini_diter) /= 0) return

    call fini%load(filename='input.ini')
    call Load_Ini(fini)
    out = Validate_Registry()
    if (trim(out) /= '') then
      write(*,'(A)') trim(out)
      error stop 1
    endif

    call Normalize_MI2_Params()
    ! relaxation read again under a start on MOSE's coarse grids
    if (hydra_MI2_obj%coarse_start) call Check_Coarse_Start()

  end subroutine Read_Inifile_MI2_Runtime


  !> The shared keys of solvers MI2 does not run (FUSS, the multigrid levels) are not required.
  subroutine Release_Unused_Keys()
    use hydra_Input_Registry, only: reg
    implicit none
    integer :: i

    do i = 1, reg%size
      if (reg%params(i)%name == 'FUSS-iter' .or. reg%params(i)%section == 'HYDRA-Multigrid') &
        reg%params(i)%required = .false.
    enddo

  end subroutine Release_Unused_Keys


  !> Allocate the MI2 arrays at the length the file carries and record whether each key is present.
  subroutine Scan_Ini_MI2(fini)
    use Finer, only: file_ini
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydra_Backend_INI,       only: Count_Values_Checked
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: n

    associate(o => hydra_MI2_obj)

    ! --- ICE-iter, IGLOO-iter: an absent key disables its solver and leaves a one-entry list ---
    o%run_ICE = fini%has_option(option_name='ICE-iter')
    n = 1
    if (o%run_ICE) n = max(Count_Values_Checked(fini,'HYDRA-Parameters','ICE-iter'), 1)
    if (allocated(o%ICE_iter)) deallocate(o%ICE_iter)
    allocate(o%ICE_iter(1:n))

    o%run_IGLOO = fini%has_option(option_name='IGLOO-iter')
    n = 1
    if (o%run_IGLOO) n = max(Count_Values_Checked(fini,'HYDRA-Parameters','IGLOO-iter'), 1)
    if (allocated(o%IGLOO_iter)) deallocate(o%IGLOO_iter)
    allocate(o%IGLOO_iter(1:n))

    ! --- relaxation-coeff: one value, or one per cycle ---
    o%n_relax_given = 0
    if (fini%has_option(option_name='relaxation-coeff')) &
      o%n_relax_given = Count_Values_Checked(fini,'HYDRA-Parameters','relaxation-coeff')
    if (allocated(o%relCoeff_raw)) deallocate(o%relCoeff_raw)
    allocate(o%relCoeff_raw(1:max(o%n_relax_given,1)))

    end associate

  end subroutine Scan_Ini_MI2


  !> Cross-key normalisation and validation after the read: IGLOO mode clamp, relaxation
  !  coefficients per cycle, cycle-length agreement with MOSE-iter.
  subroutine Normalize_MI2_Params()
    use hydra_Config_Types_m,    only: hydra_obj_sim_param
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    integer :: i

    associate(o => hydra_MI2_obj)

    ! --- IGLOO mode clamp: anything <1 is OFF, anything >1 is Eulerian coupling ---
    if (o%run_IGLOO) then
      do i = 1, size(o%IGLOO_iter)
        if     (o%IGLOO_iter(i) < 1) then; o%IGLOO_iter(i) = 0
        elseif (o%IGLOO_iter(i) > 1) then; o%IGLOO_iter(i) = 2
        endif
      enddo
    else
      o%IGLOO_iter = 0
    endif

    if (.not. o%run_ICE) o%ICE_iter = 0

    ! --- relaxation coefficients, one per cycle ---
    if (allocated(o%relCoeff)) deallocate(o%relCoeff)
    allocate(o%relCoeff(1:size(o%IGLOO_iter)))

    if (o%relaxation) then
      if (o%n_relax_given <= 1) then
        ! absent or a single scalar -> broadcast
        o%relCoeff = o%relCoeff_raw(1)
      elseif (o%n_relax_given == size(o%IGLOO_iter)) then
        o%relCoeff = o%relCoeff_raw(1:size(o%IGLOO_iter))
      else
        write(*,'(A)') '[ERROR] HYDRA-Parameters: relaxation-coeff must have either 1 value or '// &
                       'one per IGLOO-iter cycle.'
        error stop 1
      endif
    else
      ! relaxation off: coefficient 1
      o%relCoeff = 1._R8
    endif

    ! --- cycle-length agreement: every list has one entry per cycle, as MOSE-iter ---
    call Same_Length_As_MOSE_Iter( 'ICE-iter',   size(o%ICE_iter),   o%run_ICE )
    call Same_Length_As_MOSE_Iter( 'IGLOO-iter', size(o%IGLOO_iter), o%run_IGLOO )
    call Same_Length_As_MOSE_Iter( 'N-coupling', size(hydra_obj_sim_param%N_Coupling), .true. )

    end associate

  contains

    !> Refuse a per-cycle list whose length is not that of MOSE-iter.
    subroutine Same_Length_As_MOSE_Iter(key, n, active)
      character(len=*), intent(in) :: key
      integer,          intent(in) :: n
      logical,          intent(in) :: active

      if (.not. active .or. n == size(hydra_obj_sim_param%MOSE_iter)) return
      write(*,'(A,I0,A,I0)') '[ERROR] HYDRA-Parameters: '//key//' has ', n, &
                             ' values but MOSE-iter has ', size(hydra_obj_sim_param%MOSE_iter)
      error stop 1

    end subroutine Same_Length_As_MOSE_Iter

  end subroutine Normalize_MI2_Params


  !> Echo the coupling schedule exactly as MI2 resolved it, after clamping and defaulting.
  subroutine Report_Schedule()
    use hydra_Config_Types_m,    only: hydra_obj_sim_param
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    character(len=*), parameter :: INTS = '(A,*(I0," "))'

    write(*,*)
    write(*,'(A)') ' HYDRA-MI2 coupling schedule'
    write(*,'(A,I0)')  ' -  cycles (from MOSE-iter)   --> ', size(hydra_obj_sim_param%MOSE_iter)
    write(*,'(A,L1)')  ' -  run ICE                   --> ', hydra_MI2_obj%run_ICE
    write(*,'(A,L1)')  ' -  run IGLOO                 --> ', hydra_MI2_obj%run_IGLOO
    write(*,'(A,L1)')  ' -  relaxation                --> ', hydra_MI2_obj%relaxation
    write(*,INTS)      ' -  MOSE-iter                 --> ', hydra_obj_sim_param%MOSE_iter
    write(*,INTS)      ' -  N-coupling                --> ', hydra_obj_sim_param%N_Coupling
    if (hydra_MI2_obj%run_ICE) then
      write(*,INTS)    ' -  ICE-iter                  --> ', hydra_MI2_obj%ICE_iter
    endif
    if (hydra_MI2_obj%run_IGLOO) then
      write(*,INTS)    ' -  IGLOO-iter (0/1/2)        --> ', hydra_MI2_obj%IGLOO_iter
      write(*,'(A,*(F6.3," "))') ' -  relaxation-coeff          --> ', hydra_MI2_obj%relCoeff
    endif
    write(*,*)

  end subroutine Report_Schedule


  ! ---------------------------------------------------------------------------------------------
  ! The solvers' own sections, after their setup
  ! ---------------------------------------------------------------------------------------------

  !> MOSE's grid levels, after MOSE's setup. A start on its coarse grids (a level<m>-iter > 0, m > 1) runs as MOSE
  !  alone runs it: the coupler hands the condensed solvers the finest grid filled from the level MOSE steps and adds
  !  their sources to that level (coarse_start). Otherwise MOSE steps its finest grid, the one the coupler hands the
  !  condensed solvers, from the first iteration (its multigrid cycle alone).
  subroutine Check_MOSE_Grid_Levels ()
    use MOSE_Config_Types_m,     only: MOSE_multigrid => obj_multigrid
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none

    if ( MOSE_multigrid%MGL <= 1 ) return
    call Check_MOSE_Level1_Iter()
    hydra_MI2_obj%coarse_start = any( MOSE_multigrid%iter_threshold(2:MOSE_multigrid%MGL) > 0 )
    if ( hydra_MI2_obj%coarse_start ) then
      call Check_Coarse_Start()
      write(*,'(A,I0,A)') ' MI2: MOSE starts on its coarse grids (levels = ', MOSE_multigrid%MGL, '): the coupler '// &
        'hands IGLOO and ICE the gas of the finest grid, filled from the level MOSE steps, and adds their sources to '// &
        'that level summed over its sub-cells'
      return
    endif
    MOSE_multigrid%MG_level = 1
    if ( trim(MOSE_multigrid%cycle) /= 'none' ) write(*,'(A,I0,A)') ' MI2: MOSE''s multigrid cycle on ', &
      MOSE_multigrid%MGL, ' levels: MOSE steps its finest grid from the first iteration; a coarse level takes the '// &
      'condensed solvers'' sources in the residual MOSE restricts to it'

  end subroutine Check_MOSE_Grid_Levels


  !> MOSE's grid levels without level1-iter, refused: at MOSE's set-up, and after MOSE reads input.ini again.
  subroutine Check_MOSE_Level1_Iter ()
    use MOSE_Config_Types_m, only: MOSE_multigrid => obj_multigrid
    implicit none

    if ( MOSE_multigrid%MGL <= 1 ) return
    if ( MOSE_multigrid%iter_threshold(1) <= 0 ) then
      write(*,'(A)') ' [ERROR] MI2: [MOSE-Multigrid] level1-iter is not set: with levels > 1 it is the iteration '// &
        'limit of MOSE''s finest grid, and with its default (0) the gas makes one iteration per coupling step'
      error stop 1
    endif

  end subroutine Check_MOSE_Level1_Iter


  !> A start on MOSE's coarse grids where the coupler does not follow it: a time-accurate run (MOSE does not carry the
  !  physical time across levels), and IGLOO's relaxation (its store is cell by cell of the level MOSE steps). At
  !  MOSE's set-up, and after input.ini is read again.
  subroutine Check_Coarse_Start ()
    use MOSE_Config_Types_m,     only: MOSE_time_scheme => obj_time_scheme, MOSE_sim => obj_sim_param
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    character(len=*), parameter :: E = ' [ERROR] MI2: [MOSE-Multigrid] level<m>-iter > 0 (m > 1) starts MOSE on its '// &
                                       'coarse grids, '

    if ( MOSE_time_scheme%time_accurate .or. MOSE_sim%HYDRA_time_accurate ) then
      write(*,'(A)') E//'and the run is time-accurate ([MOSE-Numerics] time-accurate = true): MOSE does not carry the '// &
        'physical time from one level to the next. Set level<m>-iter = 0 for every m > 1 (levels = 1), or run steady'
      error stop 1
    endif
    if ( hydra_MI2_obj%run_IGLOO .and. hydra_MI2_obj%relaxation ) then
      write(*,'(A)') E//'and relaxation = true blends each coupling step''s IGLOO sources with the last step''s, cell by '// &
        'cell of the level MOSE steps, which changes at each switch of level. Set relaxation = false, or '// &
        'level<m>-iter = 0 for every m > 1'
      error stop 1
    endif

  end subroutine Check_Coarse_Start


  !> ICE's grid levels, after ICE's setup: the coupler hands the gas to ICE's finest level and takes ICE's source from
  !  it, so ICE runs on one level.
  subroutine Check_ICE_Grid_Levels ()
    use ICE_Config_Types_m, only: ICE_multigrid => obj_multigrid
    implicit none

    if ( ICE_multigrid%MGL <= 1 ) return
    write(*,'(A,I0,A)') ' [ERROR] MI2: [ICE-Multigrid] levels = ', ICE_multigrid%MGL, ': the coupler exchanges with '// &
      'ICE''s finest grid only, which a coarse level does not step. Use levels = 1'
    error stop 1

  end subroutine Check_ICE_Grid_Levels


  !> After every solver's setup: a cycle that advances MOSE and ICE together while either is time-accurate (they do
  !  not share a time step) and the Eulerian coupling of an IGLOO material that converts are refused; a material both
  !  condensed solvers carry and burn differently by their sections' defaults is warned about.
  subroutine Check_Coupled_Inputs ( IGLOO )
    use hydra_Config_Types_m,    only: hydra_obj_sim_param
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydraMI2_Chemistry,      only: chemIGLOO
    use MOSE_Config_Types_m,     only: MOSE_time_scheme => obj_time_scheme
    use ICE_Config_Types_m,      only: ICE_time_scheme  => obj_time_scheme
    use IGLOO_module,            only: obj_IGLOO
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    integer :: c

    if ( hydra_MI2_obj%run_ICE .and. (MOSE_time_scheme%time_accurate .or. ICE_time_scheme%time_accurate) ) then
      do c = 1, size(hydra_obj_sim_param%MOSE_iter)
        if ( hydra_obj_sim_param%MOSE_iter(c) > 0 .and. hydra_MI2_obj%ICE_iter(c) > 0 ) then
          write(*,'(A,I0,A)') '[ERROR] [MI2] cycle ', c, ' advances MOSE and ICE together with time-accurate on: they '// &
            'do not share a time step, so their physical times would part. Advance them in separate cycles, or set '// &
            'time-accurate = false in [MOSE-Numerics] and [ICE-Numerics].'
          error stop 1
        endif
      enddo
    endif

    if ( chemIGLOO%active .and. any(hydra_MI2_obj%IGLOO_iter == 2) ) then
      write(*,'(A)') ' [ERROR] MI2: an IGLOO material evaporates or burns and a cycle has IGLOO-iter = 2: the '// &
        'Eulerian coupling hands the gas momentum and energy only. Use IGLOO-iter = 1 in every cycle that runs IGLOO'
      error stop 1
    endif

    if ( hydra_MI2_obj%run_ICE .and. hydra_MI2_obj%run_IGLOO ) call Check_Burning_Defaults( IGLOO )

  end subroutine Check_Coupled_Inputs


  !> A material name that ICE and IGLOO both carry and burn differently, where a phase file does not say the model
  !  (no combustion token on the material line: the solver's section decided): a warning.
  subroutine Check_Burning_Defaults ( IGLOO )
    use ICE_Global_m,       only: ICE_phase_prefix
    use ICE_Config_Types_m, only: ICE_condensed => obj_condensed
    use IGLOO_module,       only: obj_IGLOO
    use IGLOO_variables,    only: IGLOO_phase_prefix
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    ! the two solvers number their models alike (COMB_BECKSTEAD = 1, COMB_DIRECT = 2)
    character(len=9), parameter :: WORD(0:2) = ['none     ', 'Beckstead', 'direct   ']
    integer :: mi, mg
    logical :: saidICE, saidIGLOO

    do mi = 1, size(ICE_condensed)
      do mg = 1, size(IGLOO%material)
        if ( trim(adjustl(ICE_condensed(mi)%name)) /= trim(adjustl(IGLOO%material(mg)%matName)) ) cycle
        if ( IGLOO%material(mg)%ngroups == 0 ) cycle
        if ( ICE_condensed(mi)%combSelect == IGLOO%material(mg)%combSelect ) cycle
        saidICE   = Line_Says_Model( 'INPUT/'//trim(ICE_phase_prefix)//'phase.txt',   ICE_condensed(mi)%name )
        saidIGLOO = Line_Says_Model( 'INPUT/'//trim(IGLOO_phase_prefix)//'phase.txt', IGLOO%material(mg)%matName )
        if ( saidICE .and. saidIGLOO ) cycle
        write(*,'(A)') ' [WARNING] MI2: material '//trim(adjustl(ICE_condensed(mi)%name))//' burns as '// &
          trim(WORD(ICE_condensed(mi)%combSelect))//' in ICE and as '//trim(WORD(IGLOO%material(mg)%combSelect))// &
          ' in IGLOO, and '//trim(merge('neither phase file says    ', 'one phase file does not say', &
          .not. (saidICE .or. saidIGLOO)))//' its model: the solver''s section decided. Write combustion in both '// &
          '[GPB-Phase*] sections (none for a material that does not burn) to say it is meant'
      enddo
    enddo

  contains

    !> True when the line of the material in a phase file carries a combustion token.
    logical function Line_Says_Model ( file, material ) result(says)
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
        read(line,*,iostat=ios) first
        if ( ios /= 0 ) cycle
        if ( trim(first) /= trim(adjustl(material)) ) cycle
        says = index(line, ' combustion=') > 0
        exit
      enddo
      close(u)
    end function Line_Says_Model

  end subroutine Check_Burning_Defaults

end module hydraMI2_Read_Ini
