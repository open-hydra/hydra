!> hydra-MI2's mass channels: the mass a law moves between two solvers' entities, read from [MI2-Channels].
!
!  One line per channel:
!      <name> = <from> <to> <law> [key=value ...] [birth-d=<m>]
!  with <from> and <to> each gas:<species>, ICE:<material> or IGLOO:<material>, and <law> a law of hydraMI2_Laws with
!  its constants. The receiver takes the mass into the drops already in a cell (absorption), or, with birth-d, as new
!  drops of that diameter. Routes: gas -> ICE, gas -> IGLOO (condensation, collection of a species), ICE -> IGLOO and
!  IGLOO -> ICE (a transfer between the two condensed solvers).
!
!  The gas keeps the books, as for every route of the coupler (hydraMI2_Mod_Coupling): the mass leaves the donor at
!  its state and enters the gas's ledger at the donor's energy plus its datum offset; it leaves the ledger into the
!  receiver at the receiver's energy plus its own offset. The ledger's difference (the latent heat of a condensation,
!  the datum of a transfer) goes into the gas. A gas donor loses the species itself.
!
!  IGLOO takes in, and gives up, its share one sweep after it is asked (it splits a cell's mass among the parcels it
!  held there). So where IGLOO is one end, the other end follows what IGLOO reports it took or gave: the channel is
!  exact cell by cell. Per coupling step:
!    Law_Step_Start   before MOSE's block: the laws' rates on the state the step starts from; ICE's and IGLOO's
!                     imports set; the exact ends booked
!    Law_After_IGLOO  after an IGLOO sweep: what IGLOO took or gave, handed to the other end and booked (these books
!                     stay in force, as IGLOO's sources do, until its next sweep)
!    LawxMOSE         MOSE's callback: the books into MOSE's residual
!  Rates are per second. A channel with an ICE end needs a steady MOSE or ICE (the coupler's rates: a converged field
!  holds per unit time); gas -> IGLOO runs time-accurate too, as IGLOO's own sources do.
module hydraMI2_Law_Channels
  use iso_fortran_env, only: R8 => real64
  use hydraMI2_Laws,   only: law_t, law_cell_t

  implicit none
  private

  public :: Read_Law_Channels, Law_Step_Start, Law_After_IGLOO, LawxMOSE, Laws_Handed, Report_Channels, lawsOn

  integer, parameter :: AT_GAS = 0, AT_ICE = 1, AT_IGLOO = 2
  character(len=5), parameter :: SOLVER(0:2) = ['gas  ', 'ICE  ', 'IGLOO']

  !> One end of a channel: where (AT_*), which entity (gas: species; ICE: family; IGLOO: material), the ICE family's
  !  material, its datum offset there.
  type :: end_t
    integer  :: at = AT_GAS, id = 0, mat = 0
    real(R8) :: dh = 0._R8
    character(len=64) :: name = ''
  end type end_t

  !> A channel: its name, ends, law, birth diameter (0: absorption), and its books per block, cell by cell, per unit
  !  volume and time: g(0) the datum part, g(1:np) what the gas gains in the residual's layout (a gas donor's species,
  !  momentum, energy), g(np+1) the mass the channel moves.
  type :: book_t
    real(R8), allocatable :: g(:,:,:,:)
  end type book_t
  type :: law_channel_t
    character(len=64) :: name = ''
    type(end_t)  :: from, to
    type(law_t)  :: law
    real(R8)     :: d = 0._R8
    type(book_t), allocatable :: blk(:)
  end type law_channel_t

  type(law_channel_t), allocatable :: ch(:)
  logical, protected :: lawsOn = .false.

  !> The volume of every cell of MOSE's finest grid, per block (the books' cells), for a coarse level's sums.
  type(book_t), allocatable :: volFine(:)

  character(len=*), parameter :: E_ = ' [ERROR] MI2: [MI2-Channels] '

contains

  ! ---------------------------------------------------------------------------------------------
  ! Setup
  ! ---------------------------------------------------------------------------------------------

  !> Read [MI2-Channels] after every solver's setup: each line's ends, law and birth diameter; refuse what the coupler
  !  cannot move; the books opened over MOSE's blocks.
  subroutine Read_Law_Channels ( MOSE_domain, IGLOO )
    use Finer,                   only: file_ini
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: np
    use IGLOO_module,            only: obj_IGLOO
    use hydraMI2_Laws,           only: Law_Setup
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    type(obj_IGLOO),        intent(in) :: IGLOO
    type(file_ini) :: fini
    character(len=:), allocatable :: pair(:), why
    character(len=256) :: tok(32)
    character(len=64)  :: keys(32), values(32)
    type(law_channel_t) :: c
    integer :: n, k, nk, e, b

    allocate(ch(0))
    ! the books' cells, the finest grid's, and their volume
    allocate(volFine(size(MOSE_domain%blk)))
    do b = 1, size(MOSE_domain%blk)
      associate (d => MOSE_domain%blk(b)%dim)
      allocate(volFine(b)%g(1, d(1), d(2), d(3)))
      volFine(b)%g(1,:,:,:) = MOSE_domain%blk(b)%vol(1:d(1), 1:d(2), 1:d(3))
      end associate
    enddo
    call fini%load(filename='input.ini')
    if (.not. fini%has_section(section_name='MI2-Channels')) return
    do while (fini%loop(section_name='MI2-Channels', option_pairs=pair))
      c = law_channel_t()
      c%name = trim(pair(1))
      tok = ''
      read(pair(2),*,iostat=k) tok
      n = count(len_trim(tok) > 0)
      if (n < 3) call Refuse('"'//trim(pair(2))//'" does not read "<from> <to> <law> [key=value ...]"')
      c%from = End_Of(tok(1), IGLOO)
      c%to   = End_Of(tok(2), IGLOO)
      nk = 0
      do k = 4, n
        e = index(tok(k), '=')
        if (e < 2) call Refuse('"'//trim(tok(k))//'" is not key=value')
        nk = nk + 1
        keys(nk) = tok(k)(1:e-1); values(nk) = tok(k)(e+1:)
      enddo
      k = findloc(keys(1:nk), 'birth-d', dim=1)
      if (k > 0) then
        read(values(k),*,iostat=e) c%d
        if (e /= 0 .or. .not. c%d > 0._R8) call Refuse('birth-d='//trim(values(k))//' is not a diameter > 0')
      endif
      c%law = Law_Setup( trim(tok(3)), keys(1:nk), values(1:nk), why )
      if (len(why) > 0) call Refuse(why)
      call Check_Route( c, IGLOO )
      allocate(c%blk(size(MOSE_domain%blk)))
      do b = 1, size(MOSE_domain%blk)
        allocate(c%blk(b)%g(0:np+1, MOSE_domain%blk(b)%dim(1), MOSE_domain%blk(b)%dim(2), MOSE_domain%blk(b)%dim(3)), &
                 source=0._R8)
      enddo
      ch = [ch, c]
      write(*,'(A)') ' [MI2] channel '//trim(c%name)//': '//trim(SOLVER(c%from%at))//' '//trim(c%from%name)//' -> '// &
        trim(SOLVER(c%to%at))//' '//trim(c%to%name)//', law '//trim(c%law%name)// &
        trim(merge(' (born as new drops)', ' (absorbed)         ', c%d > 0._R8))
      write(*,'(A,ES15.7,A,ES15.7,A)') '       datum offsets: donor ', c%from%dh, ', receiver ', c%to%dh, ' J/kg'
    enddo
    lawsOn = size(ch) > 0

  contains

    subroutine Refuse ( what )
      character(len=*), intent(in) :: what
      write(*,'(A)') E_//trim(c%name)//': '//what
      error stop 1
    end subroutine Refuse

    !> The end named t (gas:<species>, ICE:<material>, IGLOO:<material>).
    function End_Of ( t, IGLOO ) result(x)
      use FLINT_Lib_Thermodynamic, only: ns, species_names
      use ICE_Global_m,            only: ngroups, mat_of
      use ICE_Config_Types_m,      only: ICE_condensed => obj_condensed
      use IGLOO_Mod_Exchange,      only: material_id
      use hydraMI2_Config_Types_m, only: hydra_MI2_obj
      character(len=*), intent(in) :: t
      type(obj_IGLOO),  intent(in) :: IGLOO
      type(end_t) :: x
      integer :: col, s, g

      col = index(t, ':')
      if (col < 2) call Refuse('"'//trim(t)//'" is not gas:<species>, ICE:<material> or IGLOO:<material>')
      x%name = t(col+1:)
      select case (t(1:col-1))
      case ('gas')
        x%at = AT_GAS
        do s = 1, ns
          if (trim(adjustl(species_names(s))) == trim(x%name)) x%id = s
        enddo
        if (x%id == 0) call Refuse(trim(x%name)//' is not a species of the gas')
      case ('ICE')
        x%at = AT_ICE
        if (.not. hydra_MI2_obj%run_ICE) call Refuse(trim(t)//': ICE is not run')
        do g = 1, ngroups
          if (trim(ICE_condensed(mat_of(g))%name) /= trim(x%name)) cycle
          if (x%id > 0) call Refuse(trim(x%name)//' has several ICE families: a channel goes to one')
          x%id = g; x%mat = mat_of(g)
        enddo
        if (x%id == 0) call Refuse(trim(x%name)//' has no ICE family')
      case ('IGLOO')
        x%at = AT_IGLOO
        if (.not. hydra_MI2_obj%run_IGLOO) call Refuse(trim(t)//': IGLOO is not run')
        x%id = material_id(IGLOO, x%name)
        if (x%id == 0) call Refuse(trim(x%name)//' is not an IGLOO material')
      case default
        call Refuse('"'//trim(t)//'" names no solver: gas, ICE or IGLOO')
      end select
    end function End_Of

    !> The routes the coupler moves, the law's needs, and each condensed end's datum offset.
    subroutine Check_Route ( c, IGLOO )
      use hydraMI2_Config_Types_m, only: hydra_MI2_obj
      use hydraMI2_Mod_Coupling,   only: handRates
      use ICE_Global_m,            only: solid_of
      type(law_channel_t), intent(inout) :: c
      type(obj_IGLOO),     intent(in)    :: IGLOO
      if (c%to%at == AT_GAS) call Refuse('a channel ends in ICE or IGLOO: the gas receives a condensed mass through '// &
                                         'the record of its process (products.txt)')
      if (c%from%at == c%to%at) call Refuse('a channel joins two solvers')
      if (c%law%needsDonor .and. c%from%at == AT_IGLOO) &
        call Refuse('the law '//trim(c%law%name)//' needs the donor''s mass in a cell, which IGLOO does not give')
      if (hydra_MI2_obj%relaxation .and. (c%from%at == AT_IGLOO .or. c%to%at == AT_IGLOO)) &
        call Refuse('relaxation = true blends two sweeps of IGLOO''s sources, and a channel''s books are not blended')
      if (c%to%at == AT_ICE .and. c%d == 0._R8) then
        if (solid_of(c%to%id)) call Refuse('the ICE family of '//trim(c%to%name)//' solidifies: an absorbed mass '// &
                                           'would change its solid fraction (give birth-d: new drops)')
      endif
      if (c%from%at /= AT_GAS) c%from%dh = Offset(c%from, IGLOO)
      c%to%dh = Offset(c%to, IGLOO)
      ! MOSE and ICE never step in the same cycle when both are time-accurate: a rate in force would not be the
      ! same amount at the two ends
      if ((c%from%at == AT_ICE .or. c%to%at == AT_ICE) .and. .not. handRates) &
        call Refuse('a channel with an ICE end needs a steady MOSE or ICE (time-accurate = false): its rate is per '// &
                    'second, and a time-accurate MOSE and ICE do not step together')
    end subroutine Check_Route

    !> A condensed end's datum offset: the material's formation enthalpy (a record of products.txt) less the energy its
    !  solver credits it at T_datum_ref.
    real(R8) function Offset ( x, IGLOO )
      use hydraMI2_Chemistry,      only: Formation_Enthalpy, T_datum_ref
      use hydraMI2_Config_Types_m, only: hydra_MI2_obj
      use IGLOO_Mod_Exchange,      only: absorption_packet
      use ICE_Mod_Exchange,        only: ICE_absorption_packet => absorption_packet
      type(end_t),     intent(in) :: x
      type(obj_IGLOO), intent(in) :: IGLOO
      real(R8) :: h, pk(5)
      Offset = 0._R8
      if (.not. Formation_Enthalpy(x%name, h)) &
        call Refuse(trim(x%name)//' has no formation enthalpy in products.txt (no record names it): the coupler '// &
                    'cannot book its energy against the gas''s')
      if (.not. hydra_MI2_obj%datum_correction) return
      if (x%at == AT_ICE) then
        pk = ICE_absorption_packet(x%id, 1._R8, [0._R8, 0._R8, 0._R8], T_datum_ref)
      else
        pk = absorption_packet(IGLOO, x%id, 1._R8, [0._R8, 0._R8, 0._R8], T_datum_ref)
      endif
      Offset = h - pk(5)
    end function Offset

  end subroutine Read_Law_Channels


  ! ---------------------------------------------------------------------------------------------
  ! Per coupling step, outside any parallel region
  ! ---------------------------------------------------------------------------------------------

  !> Before MOSE's block: every ICE family's and IGLOO's imports cleared (the generated mass adds its own after the
  !  block); then each channel's rate on the state the step starts from, cell by cell, set at its ends:
  !    gas -> ICE     ICE takes it in now: booked now
  !    gas -> IGLOO   IGLOO is offered it (or born now: booked now); what it takes is booked after its sweep
  !    ICE -> IGLOO   IGLOO is offered it at the drops' state (or born now, ICE giving it up now: booked now)
  !    IGLOO -> ICE   IGLOO is asked to give it up; what it gives goes to ICE after its sweep
  !  A cell whose receiver must hold drops to take a mass in (ICE absorption) and holds none moves nothing.
  subroutine Law_Step_Start ( MOSE_domain, ICE_sim, IGLOO )
    use MOSE_Advanced_Types_m,   only: MOSE_domain_type
    use MOSE_Global_m,           only: nsc, np
    use ICE_Advanced_Types_m,    only: ICE_simulation_type
    use ICE_Global_m,            only: ngroups, nbase, rho_empty
    use ICE_Mod_Exchange,        only: ICE_clear => clear_import_rates, ICE_remove => add_removal_rate
    use IGLOO_module,            only: obj_IGLOO
    use IGLOO_Mod_Exchange,      only: IGLOO_clear => clear_import_rates, add_absorption_rate, add_removal_rate, &
                                       add_birth_rate, absorption_packet
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    use hydraMI2_Mod_Coupling,   only: Gas_State_From_Prim
    use hydraMI2_Laws,           only: Law_Rate
    implicit none
    type(MOSE_domain_type),    intent(in)    :: MOSE_domain
    type(ICE_simulation_type), intent(inout) :: ICE_sim
    type(obj_IGLOO),           intent(inout) :: IGLOO
    character(len=:), allocatable :: why
    type(law_cell_t) :: s
    real(R8) :: rhoi(nsc), R, gam, mu, kappa, rate, u(3), T, vol, taken(6), pd(5), pr(5)
    integer  :: n, b, i, j, k, g

    if (hydra_MI2_obj%run_ICE) then
      do g = 1, ngroups
        call ICE_clear( ICE_sim%domain(1), g )
      enddo
    endif
    if (hydra_MI2_obj%run_IGLOO) call IGLOO_clear( IGLOO )

    do n = 1, size(ch)
      associate (c => ch(n))
      do b = 1, size(MOSE_domain%blk)
        if (.not. Lagged(c)) c%blk(b)%g = 0._R8
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          call Gas_State_From_Prim( MOSE_domain%blk(b)%P(:,i,j,k), rhoi, s%u, s%rho, R, gam, s%T, mu, kappa )
          s%p   = MOSE_domain%blk(b)%P(np,i,j,k)
          vol   = MOSE_domain%blk(b)%vol(i,j,k)
          s%vol = vol
          s%donor = -1._R8; s%receiver = -1._R8
          if (c%from%at == AT_GAS) s%donor = rhoi(c%from%id)
          if (c%from%at == AT_ICE) s%donor = ICE_sim%domain(1)%blk(b)%cond_phase(c%from%id)%prim(1,i,j,k)
          if (c%to%at == AT_ICE)   s%receiver = ICE_sim%domain(1)%blk(b)%cond_phase(c%to%id)%prim(1,i,j,k)
          rate = Law_Rate( c%law, s )
          if (.not. rate > 0._R8) cycle
          ! ICE absorbs only into drops it holds
          if (c%to%at == AT_ICE .and. c%d == 0._R8 .and. .not. s%receiver > rho_empty) cycle
          ! the state the mass leaves at: the gas's, or the ICE drops'
          u = s%u; T = s%T
          if (c%from%at == AT_ICE) then
            if (.not. s%donor > rho_empty) cycle
            u = ICE_sim%domain(1)%blk(b)%cond_phase(c%from%id)%prim(2:4,i,j,k)
            T = ICE_sim%domain(1)%blk(b)%cond_phase(c%from%id)%prim(nbase(c%from%id)-1,i,j,k)
          endif

          select case (c%to%at)
          case (AT_ICE)
            if (c%from%at == AT_GAS) then
              pr = ICE_Take( c, ICE_sim%domain(1)%blk(b), i, j, k, rate, u, T )
              call Book( c, b, i, j, k, -rate, -pr, -c%to%dh*rate, rate )
            else
              ! IGLOO -> ICE: IGLOO gives up its share in its sweep (Law_After_IGLOO)
              call add_removal_rate( IGLOO, c%from%id, b, i, j, k, rate*vol, why )
              call Stop_If( why )
            endif
          case (AT_IGLOO)
            if (c%d > 0._R8) then
              call add_birth_rate( IGLOO, c%to%id, b, i, j, k, rate*vol, u, T, c%d, why )
              call Stop_If( why )
              pr = absorption_packet( IGLOO, c%to%id, rate, u, T )
              if (c%from%at == AT_GAS) then
                call Book( c, b, i, j, k, -rate, -pr, -c%to%dh*rate, rate )
              else
                call ICE_remove( ICE_sim%domain(1)%blk(b), c%from%id, i, j, k, rate, taken, why )
                if (len(why) > 0) cycle
                pd = [taken(1:4), taken(5)]
                call Book( c, b, i, j, k, 0._R8, pd - pr, (c%from%dh - c%to%dh)*rate, rate )
              endif
            else
              call add_absorption_rate( IGLOO, c%to%id, b, i, j, k, absorption_packet(IGLOO, c%to%id, rate*vol, u, T), &
                                        why )
              call Stop_If( why )
            endif
          end select
        enddo; enddo; enddo
      enddo
      end associate
    enddo

  end subroutine Law_Step_Start


  !> After IGLOO's sweep: what it took in and gave up for each channel, cell by cell (exchange_rate's absorbed and
  !  removed, the cell's [kg/s, N, W]). gas -> IGLOO: the gas loses it; ICE -> IGLOO: ICE gives up the same mass at
  !  its drops' state; IGLOO -> ICE: ICE takes in what IGLOO gave, at its state. Booked.
  subroutine Law_After_IGLOO ( MOSE_domain, ICE_sim, IGLOO )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use ICE_Advanced_Types_m,  only: ICE_simulation_type
    use ICE_Mod_Exchange,      only: ICE_remove => add_removal_rate
    use IGLOO_module,          only: obj_IGLOO
    use IGLOO_Mod_Exchange,    only: exchange_t, exchange_rate
    implicit none
    type(MOSE_domain_type),    intent(in)    :: MOSE_domain
    type(ICE_simulation_type), intent(inout) :: ICE_sim
    type(obj_IGLOO),           intent(in)    :: IGLOO
    character(len=:), allocatable :: why
    type(exchange_t) :: x
    real(R8) :: vol, m, a(5), taken(6), pr(5), u(3), T
    integer  :: n, b, i, j, k

    do n = 1, size(ch)
      associate (c => ch(n))
      if (.not. Lagged(c)) cycle
      do b = 1, size(MOSE_domain%blk)
        c%blk(b)%g = 0._R8
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          if (c%to%at == AT_IGLOO) then
            x = exchange_rate( IGLOO, c%to%id, b, i, j, k )
            a = x%absorbed/vol
            m = a(1)
            if (.not. m > 0._R8) cycle
            if (c%from%at == AT_GAS) then
              call Book( c, b, i, j, k, -m, -a, -c%to%dh*m, m )
            else
              call ICE_remove( ICE_sim%domain(1)%blk(b), c%from%id, i, j, k, m, taken, why )
              if (len(why) > 0) cycle
              call Book( c, b, i, j, k, 0._R8, taken(1:5) - a, (c%from%dh - c%to%dh)*m, m )
            endif
          elseif (c%from%at == AT_IGLOO) then
            x = exchange_rate( IGLOO, c%from%id, b, i, j, k )
            a = x%removed/vol
            m = a(1)
            if (.not. m > 0._R8) cycle
            u = a(2:4)/m
            T = IGLOO_Temperature( IGLOO, c%from%id, a(5)/m - 0.5_R8*dot_product(u, u) )
            pr = ICE_Take( c, ICE_sim%domain(1)%blk(b), i, j, k, m, u, T )
            call Book( c, b, i, j, k, 0._R8, a - pr, (c%from%dh - c%to%dh)*m, m )
          endif
        enddo; enddo; enddo
      enddo
      end associate
    enddo

  end subroutine Law_After_IGLOO


  !> Channel c follows what IGLOO reports after its sweep (IGLOO absorbs, or gives up): its books are made after the
  !  sweep and kept until the next one, the books of IGLOO's sources in force.
  logical function Lagged ( c )
    type(law_channel_t), intent(in) :: c
    Lagged = (c%to%at == AT_IGLOO .and. c%d == 0._R8) .or. c%from%at == AT_IGLOO
  end function Lagged


  !> ICE family c%to%id of cell (i,j,k) takes in m [kg/(m^3 s)] at (u, T): absorbed, or born at c%d. Its packet.
  function ICE_Take ( c, blk, i, j, k, m, u, T ) result(pk)
    use ICE_Advanced_Types_m, only: ICE_block_type
    use ICE_Mod_Exchange,     only: absorption_packet, add_absorption_rate, birth_packet, add_birth_rate
    implicit none
    type(law_channel_t),  intent(in)    :: c
    type(ICE_block_type), intent(inout) :: blk
    integer,              intent(in)    :: i, j, k
    real(R8),             intent(in)    :: m, u(3), T
    real(R8) :: pk(5), bp(6)
    character(len=:), allocatable :: why

    if (c%d > 0._R8) then
      bp = birth_packet( c%to%id, m, u, T, c%d )
      call add_birth_rate( blk, c%to%id, i, j, k, bp )
      pk = bp(1:5)
    else
      pk = absorption_packet( c%to%id, m, u, T )
      call add_absorption_rate( blk, c%to%id, i, j, k, pk, why )
      call Stop_If( why )
    endif

  end function ICE_Take


  !> The temperature of IGLOO material mat whose enthalpy, on its datum, is h [J/kg]: cp T + hOff, or its table.
  real(R8) function IGLOO_Temperature ( IGLOO, mat, h ) result(T)
    use IGLOO_module, only: obj_IGLOO
    implicit none
    type(obj_IGLOO), intent(in) :: IGLOO
    integer,         intent(in) :: mat
    real(R8),        intent(in) :: h
    integer :: it

    associate (p => IGLOO%material(mat))
    if (.not. p%cpVariable) then
      T = (h - p%hOff)/p%cp
      return
    endif
    ! the table holds h at every integer kelvin: the last row below h, then linear within it
    it = lbound(p%hTab, 1)
    do while (it < ubound(p%hTab, 1) - 1 .and. p%hTab(it+1) <= h)
      it = it + 1
    enddo
    T = real(it, R8) + (h - p%hTab(it))/(p%hTab(it+1) - p%hTab(it))
    end associate

  end function IGLOO_Temperature


  !> Book what the gas gains in cell (i,j,k) of block b from channel c, per unit volume and time: the donor species'
  !  mass dmass (a gas donor), momentum and energy pk(2:5), the datum part datum (added to the energy too); and the
  !  mass the channel moved there.
  subroutine Book ( c, b, i, j, k, dmass, pk, datum, moved )
    use MOSE_Global_m, only: nu, nw, np
    implicit none
    type(law_channel_t), intent(inout) :: c
    integer,             intent(in)    :: b, i, j, k
    real(R8),            intent(in)    :: dmass, pk(5), datum, moved

    if (c%from%at == AT_GAS) c%blk(b)%g(c%from%id,i,j,k) = c%blk(b)%g(c%from%id,i,j,k) + dmass
    c%blk(b)%g(nu:nw,i,j,k) = c%blk(b)%g(nu:nw,i,j,k) + pk(2:4)
    c%blk(b)%g(np,i,j,k)    = c%blk(b)%g(np,i,j,k) + pk(5) + datum
    c%blk(b)%g(0,i,j,k)     = c%blk(b)%g(0,i,j,k) + datum
    c%blk(b)%g(np+1,i,j,k)  = c%blk(b)%g(np+1,i,j,k) + moved

  end subroutine Book


  subroutine Stop_If ( why )
    character(len=*), intent(in) :: why
    if (len(why) == 0) return
    write(*,'(A)') E_//why
    error stop 1
  end subroutine Stop_If


  ! ---------------------------------------------------------------------------------------------
  ! The books into MOSE, and summed for the report
  ! ---------------------------------------------------------------------------------------------

  !> MOSE's callback: each channel's books into the residual (R = R - S), times the cell volume; on a coarse level
  !  (a start on MOSE's coarse grids) the books of each coarse cell's sub-cells times their volume, summed. Inside
  !  MOSE's parallel region: orphaned !$OMP DO.
  subroutine LawxMOSE ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: np
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer  :: n, b, i, j, k, r(3), lo(3), hi(3), fi, fj, fk
    real(R8) :: s(np)

    do n = 1, size(ch)
      do b = 1, size(MOSE_domain%blk)
        if (MOSE_domain%mg_level == 1) then
          !$OMP DO COLLAPSE(3) PRIVATE(i,j,k)
          do k = 1, MOSE_domain%blk(b)%dim(3)
          do j = 1, MOSE_domain%blk(b)%dim(2)
          do i = 1, MOSE_domain%blk(b)%dim(1)
            MOSE_domain%blk(b)%R(1:np,i,j,k) = MOSE_domain%blk(b)%R(1:np,i,j,k) &
                                               - ch(n)%blk(b)%g(1:np,i,j,k)*MOSE_domain%blk(b)%vol(i,j,k)
          enddo; enddo; enddo
          !$OMP END DO
        else
          ! the sub-cells per coarse cell: the finest grid's cells (the books') over the level's
          r = [size(ch(n)%blk(b)%g, 2), size(ch(n)%blk(b)%g, 3), size(ch(n)%blk(b)%g, 4)] / MOSE_domain%blk(b)%dim
          !$OMP DO COLLAPSE(3) PRIVATE(i,j,k,lo,hi,fi,fj,fk,s)
          do k = 1, MOSE_domain%blk(b)%dim(3)
          do j = 1, MOSE_domain%blk(b)%dim(2)
          do i = 1, MOSE_domain%blk(b)%dim(1)
            lo = ([i, j, k] - 1)*r + 1
            hi = [i, j, k]*r
            s  = 0._R8
            do fk = lo(3), hi(3)
            do fj = lo(2), hi(2)
            do fi = lo(1), hi(1)
              s = s + ch(n)%blk(b)%g(1:np,fi,fj,fk)*volFine(b)%g(1,fi,fj,fk)
            enddo; enddo; enddo
            MOSE_domain%blk(b)%R(1:np,i,j,k) = MOSE_domain%blk(b)%R(1:np,i,j,k) - s
          enddo; enddo; enddo
          !$OMP END DO
        endif
      enddo
    enddo

  end subroutine LawxMOSE


  !> What the channels hand the gas, summed over the cells [per s]: the species' mass, momentum, energy, datum part.
  subroutine Laws_Handed ( MOSE_domain, mass, mom, energy, datum )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: nsc, nu, np
    implicit none
    type(MOSE_domain_type), intent(in)  :: MOSE_domain
    real(R8),               intent(out) :: mass(nsc), mom(3), energy, datum
    integer :: n, b, s, i, j, k
    real(R8) :: vol

    mass = 0._R8; mom = 0._R8; energy = 0._R8; datum = 0._R8
    do n = 1, size(ch)
      do b = 1, size(MOSE_domain%blk)
        do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
        do i = 1, MOSE_domain%blk(b)%dim(1)
          vol = MOSE_domain%blk(b)%vol(i,j,k)
          do s = 1, nsc
            mass(s) = mass(s) + ch(n)%blk(b)%g(s,i,j,k)*vol
          enddo
          mom    = mom    + ch(n)%blk(b)%g(nu:nu+2,i,j,k)*vol
          energy = energy + ch(n)%blk(b)%g(np,i,j,k)*vol
          datum  = datum  + ch(n)%blk(b)%g(0,i,j,k)*vol
        enddo; enddo; enddo
      enddo
    enddo

  end subroutine Laws_Handed


  !> One report line per channel: the mass it moves [kg/s], summed over the cells.
  subroutine Report_Channels ( MOSE_domain )
    use MOSE_Advanced_Types_m, only: MOSE_domain_type
    use MOSE_Global_m,         only: np
    implicit none
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer  :: n, b
    real(R8) :: moved

    do n = 1, size(ch)
      moved = 0._R8
      do b = 1, size(MOSE_domain%blk)
        moved = moved + sum(ch(n)%blk(b)%g(np+1,:,:,:)*MOSE_domain%blk(b)%vol(1:size(ch(n)%blk(b)%g,2), &
                                                                         1:size(ch(n)%blk(b)%g,3), 1:size(ch(n)%blk(b)%g,4)))
      enddo
      write(*,'(A,ES24.16)') ' [MI2-report] channel '//trim(ch(n)%name)//' moved [kg/s] ', moved
    enddo

  end subroutine Report_Channels

end module hydraMI2_Law_Channels
