module hydraMI2_Mod_Import
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Parameters_m, only: llen_

  implicit none
  private
  public :: Assign_Phases

  !> One entry per *phase.txt of INPUT/ found in filelist.txt.
  type :: phase_entry_t
    character(len=llen_) :: file     = ''   ! as listed, e.g. INPUT/partL-phase.txt
    character(len=llen_) :: name     = ''   ! partL   ('' for the unnamed INPUT/phase.txt)
    character(len=llen_) :: prefix   = ''   ! partL-  ('' when unnamed)
    character(len=16)    :: kind     = ''   ! 'ideal-gas' | 'condensed'
    character(len=24)    :: word     = ''   ! type word on line 1 as written by ATLAS (e.g. liquid-dispersed)
    character(len=16)    :: modeling = ''   ! token on line 1: 'lagrangian' | 'eulerian' | '' (none)
    character(len=8)     :: solver   = ''   ! 'MOSE' | 'ICE' | 'IGLOO' | '' (unused)
    character(len=8)     :: how      = ''   ! 'ini' (solver key) | 'legacy' (fallback)
  end type phase_entry_t

  !> INI section and key where each condensed solver names the ATLAS phase it solves.
  character(len=*), parameter :: ICE_SECTION   = 'ICE-Parameters'
  character(len=*), parameter :: IGLOO_SECTION = 'IGLOO-General'
  character(len=*), parameter :: PHASE_KEY     = 'phase'

contains

  !> Route every *phase.txt of INPUT/ to the solver that owns it and set the solver phase prefixes:
  !  ideal-gas -> MOSE; condensed -> ICE / IGLOO by [ICE-Parameters] / [IGLOO-General] phase
  !  (legacy fallback: the only candidate). Checks naming rules, the ATLAS modeling token and
  !  the files each solver will open.
  subroutine Assign_Phases()
    use hydra_Get_File_List_m,   only: Get_File_List
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use MOSE_Global_m,           only: MOSE_phase_prefix
    use ICE_Global_m,            only: ICE_phase_prefix, ICE_gas_by_coupler
    use IGLOO_variables,         only: IGLOO_phase_prefix
    implicit none
    character(len=llen_), allocatable :: file_list(:)
    type(phase_entry_t),  allocatable :: ph(:)
    character(len=llen_) :: keyICE, keyIGLOO
    integer :: i, j, num_files, n
    integer :: iGas, iICE, iIGLOO
    logical :: run_ICE, run_IGLOO

    run_ICE   = hydra_MI2_obj%run_ICE
    run_IGLOO = hydra_MI2_obj%run_IGLOO

    call Get_File_List(file_list, num_files)

    ! --- discover -------------------------------------------------------------------------
    allocate(ph(num_files))
    n = 0
    do i = 1, num_files
      if (index(file_list(i), 'phase.txt') == 0) cycle
      n = n + 1
      ph(n)%file = adjustl(file_list(i))
      call Split_Phase_Filename(ph(n)%file, ph(n)%name, ph(n)%prefix)
      call Read_Phase_Header(ph(n)%file, ph(n)%kind, ph(n)%word, ph(n)%modeling)
    enddo
    if (n == 0) then
      write(*,'(A)') '[ERROR] [MI2] no INPUT/'//'*phase.txt listed in filelist.txt'
      error stop 1
    endif
    ph = ph(1:n)

    ! --- naming rules ----------------------------------------------------------------------
    if (n > 1) then
      do i = 1, n
        if (len_trim(ph(i)%name) == 0) then
          write(*,'(A)') '[ERROR] [MI2] '//trim(ph(i)%file)//' is an UNNAMED phase but the case has '// &
                         'other phase files. ATLAS matches phases by substring and the empty name '// &
                         'matches every phase: give every [GPB-Phase*] a distinct name.'
          error stop 1
        endif
      enddo
      do i = 1, n
        do j = 1, n
          if (i == j) cycle
          if (trim(ph(i)%name) == trim(ph(j)%name)) then
            write(*,'(A)') '[ERROR] [MI2] two phase files share the name "'//trim(ph(i)%name)//'": '// &
                           trim(ph(i)%file)//' and '//trim(ph(j)%file)
            error stop 1
          endif
          if (index(trim(ph(i)%name), trim(ph(j)%name)) > 0) then
            write(*,'(A)') '[ERROR] [MI2] phase name "'//trim(ph(j)%name)//'" is a substring of "'// &
                           trim(ph(i)%name)//'". ATLAS (BCB/ICB) associates phases by substring '// &
                           'match, so these two would cross-associate silently. Rename one.'
            error stop 1
          endif
        enddo
      enddo
    endif

    ! --- gas --------------------------------------------------------------------------------
    iGas = 0
    do i = 1, n
      if (ph(i)%kind /= 'ideal-gas') cycle
      if (iGas /= 0) then
        write(*,'(A)') '[ERROR] [MI2] two ideal-gas phases: '//trim(ph(iGas)%file)//' and '// &
                       trim(ph(i)%file)//'. MOSE takes exactly one.'
        error stop 1
      endif
      iGas = i
      ph(i)%solver = 'MOSE'
    enddo
    if (iGas == 0) then
      write(*,'(A)') '[ERROR] [MI2] no ideal-gas phase file found (MOSE needs one)'
      error stop 1
    endif

    ! --- condensed: the solver keys -----------------------------------------------------------
    call Read_Solver_Phase_Keys(keyICE, keyIGLOO)
    iICE = 0; iIGLOO = 0
    if (run_ICE)   call Claim_By_Key(keyICE,   'ICE',   ICE_SECTION,   'eulerian',   iICE)
    if (run_IGLOO) call Claim_By_Key(keyIGLOO, 'IGLOO', IGLOO_SECTION, 'lagrangian', iIGLOO)
    if (run_ICE .and. run_IGLOO) then
      if (iICE == iIGLOO) then   ! same file claimed twice: both keys name it, or both fell back to it
        write(*,'(A)') '[ERROR] [MI2] ICE and IGLOO both resolve to '//trim(ph(iICE)%file)//'. One condensed '// &
                       'phase cannot feed both solvers (it would be injected into the gas twice): give each '// &
                       'solver its own phase via ['//ICE_SECTION//'] / ['//IGLOO_SECTION//'] '//PHASE_KEY//'.'
        error stop 1
      endif
    endif

    ! --- hand the prefixes over -----------------------------------------------------------------
    call Fits(ph(iGas), len(MOSE_phase_prefix), 'MOSE')
    MOSE_phase_prefix = trim(ph(iGas)%prefix)
    if (iICE /= 0) then
      call Fits(ph(iICE), len(ICE_phase_prefix), 'ICE')
      ICE_phase_prefix = trim(ph(iICE)%prefix)
      ! the gas phase and the species a burnt mass becomes are the coupler's
      ICE_gas_by_coupler = .true.
    endif
    if (iIGLOO /= 0) then
      call Fits(ph(iIGLOO), len(IGLOO_phase_prefix), 'IGLOO')
      IGLOO_phase_prefix = trim(ph(iIGLOO)%prefix)
    endif

    ! --- the files each solver will open ------------------------------------------------------
    if (iICE   /= 0) call Require_Files(ph(iICE),   'ICE')
    if (iIGLOO /= 0) call Require_Files(ph(iIGLOO), 'IGLOO')

    ! --- report ----------------------------------------------------------------------------------
    write(*,*)
    write(*,'(A)') ' HYDRA-MI2 phase routing'
    do i = 1, n
      if (len_trim(ph(i)%solver) == 0) then
        write(*,'(A)') ' -  '//trim(ph(i)%file)//'  ('//trim(ph(i)%word)//trim(Token_Tag(ph(i)))// &
                       ')  --> UNUSED (no enabled solver names it)'
      elseif (ph(i)%kind == 'ideal-gas') then
        write(*,'(A)') ' -  '//trim(ph(i)%file)//'  (ideal-gas)  --> MOSE   prefix "'//trim(ph(i)%prefix)//'"'
      else
        write(*,'(A)') ' -  '//trim(ph(i)%file)//'  ('//trim(ph(i)%word)//trim(Token_Tag(ph(i)))//')  --> '// &
                       trim(ph(i)%solver)//'  prefix "'//trim(ph(i)%prefix)//'"  ('// &
                       trim(merge('from its phase key   ', 'LEGACY fallback      ', ph(i)%how == 'ini'))//')'
      endif
    enddo
    write(*,*)

    deallocate(file_list, ph)

  contains

    !> Resolve the phase of `solver` from its INI key; without a key, the legacy fallback.
    subroutine Claim_By_Key(key, solver, section, expect, owner)
      character(len=*), intent(in)  :: key, solver, section, expect
      integer,          intent(out) :: owner
      integer :: k, ncand, icand

      owner = 0
      if (len_trim(key) > 0) then
        do k = 1, n
          if (ph(k)%kind /= 'condensed') cycle
          if (trim(ph(k)%name) == trim(key)) owner = k
        enddo
        if (owner == 0) then
          write(*,'(A)') '[ERROR] [MI2] ['//section//'] '//PHASE_KEY//' = '//trim(key)// &
                         ' names no condensed phase file. Condensed phases found:'
          do k = 1, n
            if (ph(k)%kind == 'condensed') write(*,'(A)') '          - '//trim(ph(k)%file)// &
                                                          '  (name "'//trim(ph(k)%name)//'")'
          enddo
          error stop 1
        endif
        ph(owner)%how = 'ini'
      else
        ! legacy: no key, take the only condensed phase whose token does not contradict the solver
        ncand = 0; icand = 0
        do k = 1, n
          if (ph(k)%kind /= 'condensed') cycle
          if (len_trim(ph(k)%modeling) > 0 .and. trim(ph(k)%modeling) /= expect) cycle
          ncand = ncand + 1; icand = k
        enddo
        if (ncand /= 1) then
          write(*,'(A,I0,A)') '[ERROR] [MI2] '//solver//' is enabled but ['//section//'] '//PHASE_KEY// &
                              ' is not set, and ', ncand, ' condensed phase file(s) could be meant. '// &
                              'Set '//PHASE_KEY//' = <ATLAS phase name> in ['//section//'].'
          error stop 1
        endif
        owner = icand
        ph(owner)%how = 'legacy'
        write(*,'(A)') ' [MI2] DEPRECATED: ['//section//'] '//PHASE_KEY//' not set; '//solver//' takes '// &
                       trim(ph(owner)%file)//' because it is the only candidate. Set '//PHASE_KEY// &
                       ' = '//trim(ph(owner)%name)//' explicitly.'
      endif

      ! the ATLAS modeling token, if any, must agree with the solver
      if (len_trim(ph(owner)%modeling) > 0 .and. trim(ph(owner)%modeling) /= expect) then
        write(*,'(A)') '[ERROR] [MI2] '//trim(ph(owner)%file)//' is tagged "modeling='// &
                       trim(ph(owner)%modeling)//'" by ATLAS but ['//section//'] '//PHASE_KEY// &
                       ' hands it to '//solver//' (expects "'//expect//'"). Fix the ATLAS [GPB-Phase*] '// &
                       'modeling entry or the solver key so they correspond.'
        error stop 1
      endif
      ph(owner)%solver = solver
    end subroutine Claim_By_Key

    !> Refuse a prefix longer than the solver's prefix buffer.
    subroutine Fits(p, buflen, solver)
      type(phase_entry_t), intent(in) :: p
      integer,             intent(in) :: buflen
      character(len=*),    intent(in) :: solver
      if (len_trim(p%prefix) > buflen) then
        write(*,'(A,I0,A)') '[ERROR] [MI2] phase name "'//trim(p%name)//'" does not fit '//solver// &
                            '''s prefix buffer (max ', buflen-1, ' characters)'
        error stop 1
      endif
    end subroutine Fits

    !> Modeling-token tag for the routing report.
    function Token_Tag(p) result(s)
      type(phase_entry_t), intent(in) :: p
      character(len=32) :: s
      if (len_trim(p%modeling) > 0) then
        s = ', modeling='//trim(p%modeling)
      else
        s = ', no modeling token'
      endif
    end function Token_Tag

  end subroutine Assign_Phases


  !> [ICE-Parameters] phase and [IGLOO-General] phase from input.ini ('' when absent).
  subroutine Read_Solver_Phase_Keys(keyICE, keyIGLOO)
    use Finer, only: file_ini
    implicit none
    character(len=*), intent(out) :: keyICE, keyIGLOO
    type(file_ini) :: fini
    integer :: error

    keyICE = ''; keyIGLOO = ''
    call fini%load(filename='input.ini')
    call fini%get(section_name=ICE_SECTION,   option_name=PHASE_KEY, val=keyICE,   error=error)
    if (error /= 0) keyICE = ''
    call fini%get(section_name=IGLOO_SECTION, option_name=PHASE_KEY, val=keyIGLOO, error=error)
    if (error /= 0) keyIGLOO = ''
    keyICE   = adjustl(keyICE)
    keyIGLOO = adjustl(keyIGLOO)
    call fini%free

  end subroutine Read_Solver_Phase_Keys


  !> INPUT/partL-phase.txt -> name 'partL', prefix 'partL-'; INPUT/phase.txt -> '', ''.
  subroutine Split_Phase_Filename(file, name, prefix)
    implicit none
    character(len=*), intent(in)  :: file
    character(len=*), intent(out) :: name, prefix
    character(len=llen_) :: base
    integer :: islash, ilen

    islash = scan(trim(file), '/', back=.true.)
    base   = file(islash+1:)
    ilen   = len_trim(base)

    if (trim(base) == 'phase.txt') then
      name = ''; prefix = ''
      return
    endif
    if (ilen <= 10 .or. base(ilen-9:ilen) /= '-phase.txt') then
      write(*,'(A)') '[ERROR] [MI2] cannot parse phase file name "'//trim(file)// &
                     '": expected INPUT/<name>-phase.txt or INPUT/phase.txt'
      error stop 1
    endif
    name = base(1:ilen-10)
    if (index(trim(name), '-') > 0) then
      write(*,'(A)') '[ERROR] [MI2] phase name "'//trim(name)//'" contains "-": ATLAS and hydra split '// &
                     'the file name on the first dash, so the prefix would be wrong. Rename the phase.'
      error stop 1
    endif
    prefix = trim(name)//'-'

  end subroutine Split_Phase_Filename


  !> Line 1 of a phase file: kind, the ATLAS type word and the optional modeling token.
  !  condensed-dispersed | liquid-dispersed | solid-dispersed -> ICE/IGLOO; solid-bulk is refused.
  subroutine Read_Phase_Header(file, kind, word, modeling)
    use strings, only: lowercase
    implicit none
    character(len=*), intent(in)  :: file
    character(len=*), intent(out) :: kind, word, modeling
    character(len=llen_) :: line, low
    integer :: u, ios, p, q

    kind = ''; word = ''; modeling = ''
    open(newunit=u, file=trim(file), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') '[ERROR] [MI2] cannot open phase file '//trim(file)
      error stop 1
    endif
    read(u, '(A)', iostat=ios) line
    close(u)
    if (ios /= 0) then
      write(*,'(A)') '[ERROR] [MI2] phase file '//trim(file)//' is empty'
      error stop 1
    endif
    low = lowercase(line)

    ! the type word is the first blank-delimited token; 'dispersed' must be tested before 'solid'
    q = 1
    do while (q <= len_trim(low) .and. low(q:q) /= ' ')
      q = q + 1
    enddo
    word = low(1:q-1)
    if (index(low, 'dispersed') > 0) then
      kind = 'condensed'
    elseif (index(low, 'ideal-gas') > 0) then
      kind = 'ideal-gas'
    elseif (index(low, 'solid') > 0) then
      write(*,'(A)') '[ERROR] [MI2] '//trim(file)//' is a SOLID BULK phase ("'//trim(word)//'"): hydra-MI2 '// &
                     'has no solid solver (use hydra-MF)'
      error stop 1
    else
      write(*,'(A)') '[ERROR] [MI2] unknown phase type on line 1 of '//trim(file)//': "'//trim(line)//'"'
      error stop 1
    endif

    ! modeling=<word>, spaces around '=' tolerated
    p = index(low, 'modeling')
    if (p > 0) then
      q = p + len('modeling')
      do while (q <= len_trim(low) .and. low(q:q) == ' ')
        q = q + 1
      enddo
      if (q > len_trim(low) .or. low(q:q) /= '=') then
        write(*,'(A)') '[ERROR] [MI2] malformed modeling token on line 1 of '//trim(file)//': "'//trim(line)//'"'
        error stop 1
      endif
      q = q + 1
      do while (q <= len_trim(low) .and. low(q:q) == ' ')
        q = q + 1
      enddo
      p = q
      do while (q <= len_trim(low) .and. low(q:q) /= ' ')
        q = q + 1
      enddo
      modeling = low(p:q-1)
      if (trim(modeling) /= 'lagrangian' .and. trim(modeling) /= 'eulerian') then
        write(*,'(A)') '[ERROR] [MI2] modeling = "'//trim(modeling)//'" on line 1 of '//trim(file)// &
                       ': expected "lagrangian" (IGLOO) or "eulerian" (ICE)'
        error stop 1
      endif
    endif

  end subroutine Read_Phase_Header


  !> Check that the files the solver will open under this prefix exist; a missing ICE initial
  !  condition only warns.
  subroutine Require_Files(p, solver)
    implicit none
    type(phase_entry_t), intent(in) :: p
    character(len=*),    intent(in) :: solver
    character(len=llen_) :: f
    logical :: ex

    f = 'INPUT/'//trim(p%prefix)//'bc.txt'
    inquire(file=trim(f), exist=ex)
    if (.not. ex) call Missing(f)

    f = 'INPUT/'//trim(p%prefix)//'properties.dat'
    inquire(file=trim(f), exist=ex)
    if (.not. ex) then
      if (solver == 'IGLOO') then
        call Missing(f)
      else
        write(*,'(A)') ' [MI2] note: '//trim(f)//' absent -- ICE uses the constant [ICE-Physics] '// &
                       'density and specific-heat with one material, and stops with several'
      endif
    endif

    if (solver == 'ICE') then
      inquire(file='INPUT/'//trim(p%prefix)//'ic.tec',   exist=ex)
      if (.not. ex) inquire(file='INPUT/'//trim(p%prefix)//'ic.szplt', exist=ex)
      if (.not. ex) inquire(file='INPUT/'//trim(p%prefix)//'ic.plt',   exist=ex)
      if (.not. ex) write(*,'(A)') ' [MI2] WARNING: no INPUT/'//trim(p%prefix)//'ic.{tec,szplt,plt} -- ICE '// &
                                   'needs it unless [ICE-Parameters] newrun = false'
    endif

  contains

    !> Abort naming the missing file.
    subroutine Missing(f)
      character(len=*), intent(in) :: f
      write(*,'(A)') '[ERROR] [MI2] '//solver//' owns phase "'//trim(p%name)//'" but '//trim(f)//' does not exist'
      error stop 1
    end subroutine Missing

  end subroutine Require_Files

end module hydraMI2_Mod_Import
