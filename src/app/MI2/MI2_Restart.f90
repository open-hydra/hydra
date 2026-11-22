!> A run of hydra-MI2 that restarts (newrun = false in a solver's section), and what ICE leaves for one.
!
!  The rule: a restart either continues the run its fields come from exactly, or it is stopped when it would not,
!  with the way out (a new run from the restored fields). A solver continues a run when it restarts and resumes
!  that run's counters (not across a change of kind, steady <-> time-accurate, nor from a field without them).
!    IGLOO keeps no parcels and no sources from one run to the next: a restart is a new IGLOO start, stopped when a
!          cycle in which IGLOO iterates continues a time-accurate run (Check_IGLOO_Restart).
!    ICE   leaves, with its fields, the hand-off file (ICE_Exchange_File): what its last coupling step handed the gas
!          and left for its next, which its field file cannot give back. A restart reads it back, or is stopped
!          when a time-accurate run cannot be continued without it (ICE_Handoff_Unused).
!  The hand-off is written at the end of every coupling step that leaves both solvers' fields on disk, and at the
!  end of a run with ICE.
module hydraMI2_Restart
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Setup_Restart, Remove_Stale_Handoff
  public :: Note_Outputs, Save_Handoff_If_Fields_Current, Write_ICE_Exchange
  public :: Made_No_Iteration

  interface
    !> libc rename(2): the hand-off's temporary file onto the hand-off, in one step
    function c_rename ( from, to ) bind(C, name='rename') result( rc )
      use, intrinsic :: iso_c_binding, only: c_char, c_int
      character(kind=c_char), intent(in) :: from(*), to(*)
      integer(c_int) :: rc
    end function c_rename
  end interface

  !> A solver that restarts and resumes the counters of the run its field comes from: this run continues that one.
  logical :: MOSE_continues = .false., ICE_continues = .false.
  !> The iteration of each solver's field on disk (-1: none of this run or of the run it resumes) and its output
  !  count: a coupling step that leaves both fields at their solver's iteration writes ICE's hand-off with them.
  integer :: MOSE_fieldIter = -1, ICE_fieldIter = -1, MOSE_outNum = 0, ICE_outNum = 0
  !> Each solver's iteration after its set-up: a restarted solver that makes no iteration writes no end-of-run output.
  integer :: MOSE_iter0 = 0, ICE_iter0 = 0

  !> A variable of the hand-off file: its name, what it holds (XS source, XH heat, XB birth, XG the removals' store,
  !  XU the relaxation time, XR a smoke family's srate, XA and XT goneRate's a and t), family, component.
  type :: xvar_t
    character(len=32) :: name = ''
    integer :: what = 0, g = 0, n = 0
  end type xvar_t
  integer, parameter :: XS = 1, XH = 2, XB = 3, XG = 4, XR = 5, XA = 6, XT = 7, XU = 8
  character(len=24), parameter :: XAUX(2) = [character(len=24) :: 'MI2_ICE_HAS_SOURCES', 'MI2_MOSE_ITERATION']

  !> What a restart that is a new start for a solver says of the run.
  character(len=*), parameter :: NOT_CONTINUED = 'No time-accurate run is continued here (a steady schedule, or '// &
    'counters not resumed): a steady run converges again, a time-accurate one is a new run from the restored fields'

contains

  ! ---------------------------------------------------------------------------------------------
  ! Setup
  ! ---------------------------------------------------------------------------------------------

  !> After the solvers' setup: whether each solver continues the run it restarts from, IGLOO's restart rule, and an
  !  ICE restart's hand-off read back (ICE_has_sources: whether it held sources MOSE may read); then the iteration
  !  and output count each solver starts from.
  subroutine Setup_Restart ( ICE_simulation, ICE_has_sources )
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use MOSE_Config_Types_m,     only: MOSE_sim => obj_sim_param, MOSE_io => obj_io
    use ICE_Config_Types_m,      only: ICE_sim => obj_sim_param, ICE_io => obj_io
    use ICE_Advanced_Types_m,    only: ICE_simulation_type
    implicit none
    type(ICE_simulation_type), intent(inout) :: ICE_simulation
    logical,                   intent(inout) :: ICE_has_sources

    MOSE_continues = .not. MOSE_sim%newrun .and. MOSE_io%restart%found .and. &
                     MOSE_sim%iter_general == MOSE_io%restart%cycle
    if ( hydra_MI2_obj%run_ICE ) ICE_continues = .not. ICE_sim%newrun .and. ICE_io%restart%found .and. &
                                                 ICE_sim%iter_general == ICE_io%restart%cycle

    ! a restart with IGLOO is a new IGLOO start: stopped when it would not continue a time-accurate cycle
    if ( hydra_MI2_obj%run_IGLOO ) call Check_IGLOO_Restart()

    ! an ICE restart takes back what ICE's last coupling step handed the gas and left for its next
    if ( hydra_MI2_obj%run_ICE ) then
      if ( .not. ICE_sim%newrun ) call Restore_ICE_Handoff( ICE_simulation, ICE_has_sources )
    endif

    ! the fields a restart read are on disk at their solver's iteration (with the counters of the run they come from)
    if ( MOSE_continues ) MOSE_fieldIter = MOSE_sim%iter_general
    if ( ICE_continues )  ICE_fieldIter  = ICE_sim%iter_general
    MOSE_iter0  = MOSE_sim%iter_general
    MOSE_outNum = MOSE_io%output_number
    if ( hydra_MI2_obj%run_ICE ) then
      ICE_iter0  = ICE_sim%iter_general
      ICE_outNum = ICE_io%output_number
    endif

  end subroutine Setup_Restart


  !> An ICE restart: the hand-off is of a coupling step of the run the fields come from only when both solvers resume
  !  its counters (MOSE may also start anew). Read back when it is; without it the first MOSE block runs without
  !  ICE's source, which ICE_Handoff_Unused allows or stops.
  subroutine Restore_ICE_Handoff ( ICE_simulation, ICE_has_sources )
    use MOSE_Config_Types_m,  only: MOSE_sim => obj_sim_param
    use ICE_Config_Types_m,   only: ICE_sim => obj_sim_param
    use ICE_Advanced_Types_m, only: ICE_simulation_type
    implicit none
    type(ICE_simulation_type), intent(inout) :: ICE_simulation
    logical,                   intent(inout) :: ICE_has_sources
    logical :: used
    character(len=:), allocatable :: why

    used = .false.
    if ( .not. ICE_continues .or. .not. (MOSE_sim%newrun .or. MOSE_continues) ) then
      why = trim(merge('ICE ', 'MOSE', .not. ICE_continues))//' restarts without the counters of the run it '// &
        'continues (a run of the other kind, steady or time-accurate, or a field without them)'
    else
      call Read_ICE_Exchange( ICE_simulation, MOSE_sim%iter_general, .not. MOSE_sim%newrun, used, ICE_has_sources, why )
      if ( used ) write(*,'(A,I0,A,I0,A)') ' [MI2] ICE restarts with what its last coupling step handed the '// &
        'gas, from '//ICE_Exchange_File()//' (ICE iteration ', ICE_sim%iter_general, ', MOSE iteration ', &
        MOSE_sim%iter_general, ')'
    endif
    if ( .not. used ) call ICE_Handoff_Unused( why )

  end subroutine Restore_ICE_Handoff


  !> A new ICE run, past every set-up refusal: a hand-off left in OUTPUT/ by an earlier run is not of the fields this
  !  run will write, so it is removed.
  subroutine Remove_Stale_Handoff ()
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use ICE_Config_Types_m,      only: ICE_sim => obj_sim_param
    implicit none
    logical :: there
    integer :: u

    if ( .not. hydra_MI2_obj%run_ICE ) return
    if ( .not. ICE_sim%newrun ) return
    inquire( file=ICE_Exchange_File(), exist=there )
    if ( .not. there ) return
    open( newunit=u, file=ICE_Exchange_File(), status='old' )
    close( u, status='delete' )
    write(*,'(A)') ' [MI2] ICE starts anew: '//ICE_Exchange_File()//' of an earlier run removed'

  end subroutine Remove_Stale_Handoff


  ! ---------------------------------------------------------------------------------------------
  ! The rules of a restart that cannot continue exactly
  ! ---------------------------------------------------------------------------------------------

  !> newrun = false with IGLOO: IGLOO keeps neither its parcels nor the sources it hands the gas from one run to the
  !  next, so it starts empty on the restored fields. Stop when a cycle in which IGLOO iterates continues a
  !  time-accurate run (a cold IGLOO start does not); otherwise one warning; nothing when IGLOO never iterates.
  subroutine Check_IGLOO_Restart ()
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use MOSE_Config_Types_m,     only: MOSE_sim => obj_sim_param
    use ICE_Config_Types_m,      only: ICE_sim => obj_sim_param
    implicit none
    character(len=:), allocatable :: why
    logical :: iceRe
    integer :: c

    iceRe = .false.
    if ( hydra_MI2_obj%run_ICE ) iceRe = .not. ICE_sim%newrun
    if ( (MOSE_sim%newrun .and. .not. iceRe) .or. all( hydra_MI2_obj%IGLOO_iter == 0 ) ) return

    call Timed_Cycle( hydra_MI2_obj%IGLOO_iter, 'IGLOO', c, why )
    if ( c > 0 ) then
      write(*,'(A,I0,A)') ' [ERROR] MI2: newrun = false with IGLOO, and IGLOO iterates in cycle ', c, ', a time-'// &
        'accurate cycle: '//why//'. IGLOO keeps no parcels and no sources from one run to the next, so it would start '// &
        'empty on the restored fields, and a cold IGLOO start does not continue a time-accurate run. To go on from '// &
        'these fields, start a new run from them: '//Way_Out()
      error stop 1
    endif

    write(*,*)
    write(*,'(A)') ' [WARNING] MI2: newrun = false with IGLOO: IGLOO keeps no parcels and no sources from one run '// &
      'to the next, so it starts empty on the restored fields: a new IGLOO start. '//NOT_CONTINUED
    write(*,*)

  end subroutine Check_IGLOO_Restart


  !> ICE restarts and its hand-off cannot be used (why): the first MOSE block runs without ICE's source. Stop when a
  !  cycle in which ICE iterates continues a time-accurate run, or when ICE iterates in no cycle (the hand-off would be
  !  the only source of ICE the gas gets: another end state); one warning otherwise.
  subroutine ICE_Handoff_Unused ( why )
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
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
    call Timed_Cycle( hydra_MI2_obj%ICE_iter, 'ICE', c, clock )
    if ( c > 0 ) then
      write(*,'(A,I0,A)') E//why//'. ICE iterates in cycle ', c, ', a time-accurate cycle: '//clock//'. Without the '// &
        'hand-off the first MOSE block would run without ICE''s source, and the run would not continue the one it '// &
        'restarts. To go on from these fields, start a new run from them: '//Way_Out()
      error stop 1
    endif
    write(*,*)
    write(*,'(A)') ' [WARNING] MI2: ICE''s hand-off is not used: '//why//'. The first MOSE block runs without ICE''s '// &
      'source. '//NOT_CONTINUED
    write(*,*)

  end subroutine ICE_Handoff_Unused


  !> The first cycle c in which IGLOO or ICE (who, iterating as iters says) iterates and continues a time-accurate
  !  run, and why it does (Cycle_Clock); c = 0 when there is none.
  subroutine Timed_Cycle ( iters, who, c, why )
    implicit none
    integer,                       intent(in)  :: iters(:)
    character(len=*),              intent(in)  :: who
    integer,                       intent(out) :: c
    character(len=:), allocatable, intent(out) :: why

    do c = 1, size(iters)
      if ( iters(c) == 0 ) cycle
      why = Cycle_Clock( c, who )
      if ( len(why) > 0 ) return
    enddo
    c = 0

  end subroutine Timed_Cycle


  !> Why cycle c, in which IGLOO or ICE (who) iterates, continues a time-accurate run ('' if it does not): the gas it
  !  hands its sources to steps in physical time and MOSE resumed its counters (the sources stay in force until the
  !  next solve, whichever cycle MOSE iterates in), or ICE iterates in it in physical time and resumed its counters.
  function Cycle_Clock ( c, who ) result( why )
    use hydra_Config_Types_m,    only: hydra_obj_sim_param
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use MOSE_Config_Types_m,     only: MOSE_time_scheme => obj_time_scheme
    use ICE_Config_Types_m,      only: ICE_time_scheme => obj_time_scheme
    implicit none
    integer,          intent(in) :: c
    character(len=*), intent(in) :: who
    character(len=:), allocatable :: why

    why = ''
    if ( MOSE_continues .and. MOSE_time_scheme%time_accurate .and. any( hydra_obj_sim_param%MOSE_iter > 0 ) ) then
      why = 'the gas '//who//' hands its sources to steps in physical time ([MOSE-Numerics] time-accurate = true) '// &
            'and continues the restored run (MOSE resumes its counters)'
    elseif ( ICE_continues ) then
      if ( hydra_MI2_obj%ICE_iter(c) > 0 .and. ICE_time_scheme%time_accurate ) &
        why = 'ICE iterates in it in physical time ([ICE-Numerics] time-accurate = true) and continues the restored '// &
              'run (ICE resumes its counters)'
    endif

  end function Cycle_Clock


  !> The way out of a restart that cannot continue: a new run from the restored fields, each copied as the initial
  !  condition (keeping the one it replaces), in the format it was written in, newrun = true in the sections of the
  !  solvers that restart. Items separated by '; '.
  function Way_Out () result( s )
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use MOSE_Global_m,           only: MOSE_phase_prefix
    use MOSE_Config_Types_m,     only: MOSE_sim => obj_sim_param, MOSE_io => obj_io
    use ICE_Global_m,            only: ICE_phase_prefix
    use ICE_Config_Types_m,      only: ICE_sim => obj_sim_param, ICE_io => obj_io
    use ICE_IO_Solution,         only: ICE_condinit => condinit
    implicit none
    character(len=:), allocatable :: s, keys
    logical :: moseRe, iceRe

    moseRe = .not. MOSE_sim%newrun
    iceRe  = .false.
    if ( hydra_MI2_obj%run_ICE ) iceRe = .not. ICE_sim%newrun
    s    = ''
    keys = ''
    if ( moseRe ) then
      s    = Field_As_IC( trim(MOSE_io%nameinit), MOSE_phase_prefix, MOSE_io%ini_format, MOSE_io%sol_format, &
                          '[MOSE-IO]' )
      keys = '[MOSE-Parameters]'
    endif
    if ( iceRe ) then
      s    = s//Field_As_IC( trim(ICE_condinit), ICE_phase_prefix, ICE_io%ini_format, ICE_io%sol_format, '[ICE-IO]' )
      keys = keys//trim(merge(' and ', '     ', moseRe))//' [ICE-Parameters]'
    endif
    s = s//'set newrun = true in '//trim(adjustl(keys))

  end function Way_Out


  !> The items of Way_Out for one restored field: keep a copy of the initial condition it replaces, copy it there (a
  !  .vtm with the block files it lists, OUTPUT/vtk/<name>Block<b>.vts, which it names relative to its directory), and
  !  the ic-format of the field's format when the two differ.
  function Field_As_IC ( field, prefix, iniFmt, solFmt, ioSection ) result( s )
    implicit none
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


  ! ---------------------------------------------------------------------------------------------
  ! During the run
  ! ---------------------------------------------------------------------------------------------

  !> After a solver's postprocess: a new output puts that solver's field on disk at its current iteration.
  subroutine Note_Outputs ()
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use MOSE_Config_Types_m,     only: MOSE_sim => obj_sim_param, MOSE_io => obj_io
    use ICE_Config_Types_m,      only: ICE_sim => obj_sim_param, ICE_io => obj_io
    implicit none

    if ( MOSE_io%output_number /= MOSE_outNum ) then
      MOSE_outNum    = MOSE_io%output_number
      MOSE_fieldIter = MOSE_sim%iter_general
    endif
    if ( hydra_MI2_obj%run_ICE .and. ICE_io%output_number /= ICE_outNum ) then
      ICE_outNum    = ICE_io%output_number
      ICE_fieldIter = ICE_sim%iter_general
    endif

  end subroutine Note_Outputs


  !> At the end of a coupling step with ICE: when both fields on disk are of this step's end, ICE's hand-off goes with
  !  them, so that a run stopped later restarts from here. With sol-overwrite = false a restart reads the unnumbered
  !  fields of the last end of a run, which the hand-off of that end matches: no hand-off of the numbered outputs.
  subroutine Save_Handoff_If_Fields_Current ( ICE_simulation, ICE_has_sources )
    use MOSE_Config_Types_m,  only: MOSE_sim => obj_sim_param, MOSE_io => obj_io
    use ICE_Config_Types_m,   only: ICE_sim => obj_sim_param, ICE_io => obj_io
    use ICE_Advanced_Types_m, only: ICE_simulation_type
    implicit none
    type(ICE_simulation_type), intent(inout) :: ICE_simulation
    logical,                   intent(in) :: ICE_has_sources

    if ( MOSE_fieldIter == MOSE_sim%iter_general .and. ICE_fieldIter == ICE_sim%iter_general &
         .and. MOSE_io%sol_overwrite .and. ICE_io%sol_overwrite ) &
      call Write_ICE_Exchange( ICE_simulation, ICE_has_sources, MOSE_sim%iter_general, quiet=.true. )

  end subroutine Save_Handoff_If_Fields_Current


  !> A restarted solver ('MOSE' | 'ICE') that made no iteration in this run: its end-of-run output is the one it
  !  restarted from, not written again.
  logical function Made_No_Iteration ( solver )
    use MOSE_Config_Types_m, only: MOSE_sim => obj_sim_param
    use ICE_Config_Types_m,  only: ICE_sim => obj_sim_param
    implicit none
    character(len=*), intent(in) :: solver

    if ( solver == 'MOSE' ) then
      Made_No_Iteration = .not. MOSE_sim%newrun .and. MOSE_sim%iter_general == MOSE_iter0
    else
      Made_No_Iteration = .not. ICE_sim%newrun .and. ICE_sim%iter_general == ICE_iter0
    endif

  end function Made_No_Iteration


  ! ---------------------------------------------------------------------------------------------
  ! The hand-off file
  ! ---------------------------------------------------------------------------------------------

  !> What ICE's last coupling step hands the gas and leaves for its next step, kept for a restart: the source of its
  !  last step's stages (or block), and the relaxation times and smoke rates its next step starts from, are not
  !  functions of the state ICE writes, so its field file cannot give them back.
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
    integer :: g

    allocate(xv(0))
    do g = 1, ngroups
      call Add( 'S', XS, ncond(g) )
      call Add( 'H', XH, 0 )
      call Add( 'B', XB, 5 )
      call Add( 'G', XG, 5 )
      call Add( 'TAU', XU, 0 )
      if (srate(g)) call Add( 'R', XR, ncond(g) )
      if (rate)     call Add( 'GA', XA, 5 )
      if (rate)     call Add( 'GT', XT, 0 )
    enddo

  contains

    !> Variable tag_p<g>, or its components tag<n>_p<g>, n = 1..ncomp.
    subroutine Add ( tag, what, ncomp )
      character(len=*), intent(in) :: tag
      integer,          intent(in) :: what, ncomp
      character(len=32) :: nm
      integer :: n

      if (ncomp == 0) then
        write(nm,'(A,A,I0)') tag, '_p', g
        xv = [xv, xvar_t(nm, what, g, 0)]
      endif
      do n = 1, ncomp
        write(nm,'(A,I0,A,I0)') tag, n, '_p', g
        xv = [xv, xvar_t(nm, what, g, n)]
      enddo

    end subroutine Add

  end subroutine Exchange_Vars


  !> Variable v of the hand-off file in cell (i,j,k) of block b: read into x, or written from x (put). An infinite
  !  relaxation time (no drag) is held as -1; a smoke family without srate there holds 0.
  subroutine Exchange_Slot ( ICE_domain, b, v, i, j, k, x, put )
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use ICE_Advanced_Types_m,  only: ICE_domain_type
    use hydraMI2_Mod_Coupling, only: goneRate
    implicit none
    type(ICE_domain_type), intent(inout) :: ICE_domain
    integer,               intent(in)    :: b, i, j, k
    type(xvar_t),          intent(in)    :: v
    real(R8),              intent(inout) :: x
    logical,               intent(in)    :: put

    associate (cp => ICE_domain%blk(b)%cond_phase(v%g), n => v%n, g => v%g)
      if (put) then
        select case (v%what)
        case (XS); cp%source(n,i,j,k)      = x
        case (XH); cp%heat(i,j,k)          = x
        case (XB); cp%birth(n,i,j,k)       = x
        case (XG); cp%gone(n,i,j,k)        = x
        case (XU); cp%tau(i,j,k)           = merge(ieee_value(1._R8, ieee_positive_inf), x, x == -1._R8)
        case (XR); cp%srate(n,i,j,k)       = x
        case (XA); goneRate(b)%a(n,g,i,j,k) = x
        case (XT); goneRate(b)%t(g,i,j,k)   = x
        end select
      else
        x = 0._R8
        select case (v%what)
        case (XS); x = cp%source(n,i,j,k)
        case (XH); x = cp%heat(i,j,k)
        case (XB); x = cp%birth(n,i,j,k)
        case (XG); x = cp%gone(n,i,j,k)
        case (XU); x = merge(-1._R8, cp%tau(i,j,k), cp%tau(i,j,k) > huge(x))
        case (XR); if (allocated(cp%srate)) x = cp%srate(n,i,j,k)
        case (XA); x = goneRate(b)%a(n,g,i,j,k)
        case (XT); x = goneRate(b)%t(g,i,j,k)
        end select
      endif
    end associate

  end subroutine Exchange_Slot


  !> The arrays ICExMOSE and Open_Removals read and those ICE's next step starts from, to ICE_Exchange_File, with ICE's
  !  iteration (CYCLE), whether ICE has produced sources and MOSE's iteration. Tecplot ASCII through ORION: 17 digits,
  !  read back exactly. quiet: no line in the log (the writes of the coupling steps). Written to <file>.part, then
  !  renamed: a run killed during the write keeps the hand-off it had.
  subroutine Write_ICE_Exchange ( ICE_sim, hasSources, iterMOSE, quiet )
    use, intrinsic :: iso_c_binding, only: c_null_char
    use ICE_Advanced_Types_m, only: ICE_simulation_type
    use ICE_Global_m,         only: ngroups
    use ICE_Config_Types_m,   only: ice_sim_param => obj_sim_param
    use Lib_ORION_data,       only: orion_data, copyORION
    use Lib_Tecplot,          only: tec_write_structured_multiblock
    use hydraMI2_Mod_Coupling, only: goneRate
    implicit none
    type(ICE_simulation_type), intent(inout) :: ICE_sim
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
          call Exchange_Slot( ICE_sim%domain(1), b, xv(v), i, j, k, io%block(b)%vars(v,i,j,k), put=.false. )
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
    use ICE_Advanced_Types_m,  only: ICE_simulation_type
    use ICE_Global_m,          only: ngroups, ncond
    use ICE_Config_Types_m,    only: ice_sim_param => obj_sim_param
    use Lib_ORION_data,        only: orion_data
    use Lib_Tecplot,           only: tec_read_structured_multiblock
    use hydraMI2_Mod_Coupling, only: Gone_Rate_Setup
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
      if (.not. Same_Mesh( io%block(b)%mesh, ICE_sim%ODP(1)%block(b)%mesh )) then
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
          call Exchange_Slot( ICE_sim%domain(1), b, xv(v), i, j, k, io%block(b)%vars(col(v),i,j,k), put=.true. )
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

    !> Two meshes of a block, both there and equal.
    logical function Same_Mesh ( a, b )
      real(R8), allocatable, intent(in) :: a(:,:,:,:), b(:,:,:,:)

      Same_Mesh = allocated(a) .and. allocated(b)
      if (Same_Mesh) Same_Mesh = all(shape(a) == shape(b))
      if (Same_Mesh) Same_Mesh = all(a == b)

    end function Same_Mesh

    function Refusal ( what ) result( w )
      character(len=*), intent(in)  :: what
      character(len=:), allocatable :: w

      w = f//' does not match the ICE set-up of this run: '//what

    end function Refusal

  end subroutine Read_ICE_Exchange

end module hydraMI2_Restart
