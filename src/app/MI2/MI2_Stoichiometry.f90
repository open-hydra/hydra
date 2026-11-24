!> What the mass a condensed solver hands the gas becomes there. ATLAS GPB writes, per condensed phase,
!  INPUT/<prefix>products.txt: for each material the record of the process it runs (evaporation, combustion), the cap
!  record of an IGLOO material that burns with an oxide cap, the release record of a cap oxide handed to the gas at
!  burnout, and the smoke a burning material makes. This module turns that file into one table per condensed solver
!  (productsIGLOO, productsICE): the gas-species rows and the enthalpy-datum offsets. Built at setup, read-only after.
module hydraMI2_Stoichiometry
  use iso_fortran_env, only: R8 => real64

  implicit none
  private

  public :: gas_products_t, productsIGLOO, productsICE
  public :: MAP_IGLOO, MAP_ICE, T_datum_ref
  public :: Load_Species_Map, Report_Stoichiometry, Build_Datum_Offsets
  public :: Smoke_Record, Set_Smoke_Datum
  public :: Cap_Rates

  !> The two condensed solvers a table is loaded for.
  integer,          parameter :: MAP_IGLOO = 1, MAP_ICE = 2
  character(len=5), parameter :: SOLVER(2) = ['IGLOO', 'ICE  ']

  !> Reference temperature of the datum offsets (FLINT standard state).
  real(R8), parameter :: T_datum_ref = 298.15_R8

  !> The process a material runs, and the energy rule of its record.
  integer, parameter :: PROC_NONE = 0, PROC_EVAP = 1, PROC_COMB = 2
  integer, parameter :: EN_NONE = 0, EN_LATENT = 1, EN_FORMATION = 2

  !> What one condensed solver's materials hand the gas, materials in the order of the solver's phase file.
  type :: gas_products_t
    !> Some material evaporates or burns: the solver hands the gas mass.
    logical :: active = .false.
    !> mapped(m): material m has a record (it evaporates or burns).
    logical,  allocatable :: mapped(:)
    !> nu(m,s) = kg of gas species s exchanged per kg of material m converted (> 0 produced, < 0 consumed); with
    !  the material's smoke the row sums to 1.
    real(R8), allocatable :: nu(:,:)
    !> dh(m) = constant enthalpy-datum offset of material m [J per kg converted]: its enthalpy on the gas's
    !  formation-inclusive datum less the one its solver credits, at T_datum_ref.
    real(R8), allocatable :: dh(:)
    !> An IGLOO material that burns with an oxide cap: capped(m). Per kg of oxide its cap gains the gas exchanges
    !  nuCap(m,s) (the row sums to 1/nuOxide(m) - 1); per kg of cap oxide handed to the gas at burnout it receives
    !  nuRel(m,s) (the row sums to 1). nuOxide(m): kg of oxide per kg of the material its cap takes. dhOxide(m): the
    !  datum offset of the oxide, as dh.
    logical,  allocatable :: capped(:)
    real(R8), allocatable :: nuCap(:,:), nuRel(:,:), nuOxide(:), dhOxide(:)
    !> The smoke a material burns into (the tail of its combustion record): kg of smoke per kg burnt, and the smoke
    !  material's datum offset (its formation enthalpy at T_datum_ref less ICE's energy of it there); 0 without one.
    !  dhSmoke is set for ICE's materials only (Set_Smoke_Datum).
    real(R8), allocatable :: nuSmoke(:), dhSmoke(:)

    !> What the products file says, kept for the setup that follows the load.
    character(len=256),             private :: file = ''
    integer,           allocatable, private :: process(:), energy(:)
    real(R8),          allocatable, private :: h298(:)     ! the material's formation-inclusive enthalpy at T_datum_ref
    integer,           allocatable, private :: prod(:)     ! the material a capped material's oxide is, else 0
    real(R8),          allocatable, private :: hprod(:)    ! that oxide's formation-inclusive enthalpy at T_datum_ref
    character(len=64), allocatable, private :: smoke(:), sphase(:)    ! the smoke a burning material makes, its phase
    real(R8),          allocatable, private :: hsm(:)                 ! the smoke's h298
    character(len=64), allocatable, private :: rsmoke(:), rsphase(:)  ! the smoke a cap oxide handed at burnout becomes
    real(R8),          allocatable, private :: rnusm(:), rhsm(:)      ! kg of it per kg of oxide handed, its h298
  end type gas_products_t

  type(gas_products_t), protected, target :: productsIGLOO, productsICE

contains

  !> The table of solver which (MAP_IGLOO | MAP_ICE).
  function Table ( which ) result(p)
    integer, intent(in)           :: which
    type(gas_products_t), pointer :: p

    if (which == MAP_IGLOO) then
      p => productsIGLOO
    else
      p => productsICE
    endif

  end function Table


  !> Read INPUT/<prefix>products.txt and fill the table of solver which: for each material the record of the process
  !  it runs (combustion if it burns, else evaporation if it evaporates); for a material that burns with an oxide cap
  !  (product > 0) also its cap record and, when its cap is handed to the gas (releases), the release record of its
  !  oxide. Stops on a material that converts without a record, a file of another gas phase, an unknown material or
  !  species, a row that does not close the mass.
  subroutine Load_Species_Map ( which, prefix, material_names, evaporates, burns, product, releases, oxide_ratio )
    use Finer,                   only: file_ini
    use FLINT_Lib_Thermodynamic, only: ns
    use MOSE_Global_m,           only: MOSE_phase_prefix
    implicit none
    integer,          intent(in) :: which                 !> MAP_IGLOO | MAP_ICE
    character(len=*), intent(in) :: prefix                !> the solver's phase prefix
    character(len=*), intent(in) :: material_names(:)
    logical,          intent(in) :: evaporates(:), burns(:)
    integer,  intent(in), optional :: product(:)          !> the material that is the oxide of a capped material, else 0
    logical,  intent(in), optional :: releases(:)         !> the cap of the material is handed to the gas at burnout
    real(R8), intent(in), optional :: oxide_ratio         !> kg of oxide per kg of material the solver's cap law holds
    character(len=*), parameter  :: E = ' [ERROR] MI2: '
    type(gas_products_t), pointer :: p
    character(len=512)  :: fname, line
    character(len=1024) :: raw
    character(len=64)   :: key, mname, proc, enword, spname, gpref, pname, sname, sphase
    integer  :: nmat, n1, ns1, u, ios, m, s, k, nterm, nsfile, im, en, error, row
    real(R8) :: h298, cmass, cmole, rowsum, cprod, hprod, snu, sh
    logical  :: exists, gasline
    ! the rows as read: (1) the process each material runs, (2) its cap, (3) the release of the oxide it is
    real(R8), allocatable :: nu_s(:,:), cap_s(:,:), rel_s(:,:)
    logical,  allocatable :: listed(:), found(:), capfound(:), relfound(:), relwanted(:)
    real(R8), allocatable :: capnu(:), smnu(:)
    character(len=64), allocatable :: capname(:)
    type(file_ini) :: fini

    nmat = size(material_names)
    n1   = max(nmat,1)
    ns1  = max(ns,1)
    allocate(nu_s(n1,ns1), cap_s(n1,ns1), rel_s(n1,ns1), source=0._R8)
    allocate(listed(n1), found(n1), capfound(n1), relfound(n1), relwanted(n1), source=.false.)
    allocate(capnu(n1), smnu(n1), source=0._R8)
    allocate(capname(n1)); capname = ''

    p => Table(which)
    call Reset_Book()
    if (present(product)) p%prod(1:nmat) = product(1:nmat)
    if (present(releases)) then
      do m = 1, nmat
        if (p%prod(m) > 0 .and. releases(m)) relwanted(p%prod(m)) = .true.
      enddo
    endif
    do m = 1, nmat
      if (burns(m)) then
        p%process(m) = PROC_COMB
      elseif (evaporates(m)) then
        p%process(m) = PROC_EVAP
      endif
    enddo
    fname  = 'INPUT/'//trim(prefix)//'products.txt'
    p%file = fname

    ! the former key of input.ini names nothing any more
    call fini%load(filename='input.ini')
    do m = 1, nmat
      raw = ''
      call fini%get(section_name='MI2-Coupling', option_name='products-'//trim(adjustl(material_names(m))), &
                    val=raw, error=error)
      if (error /= 0 .or. len_trim(raw) == 0) cycle
      write(*,'(A)') E//'[MI2-Coupling] products-'//trim(adjustl(material_names(m)))//' is no longer read: '// &
        'the gas species of a material come from '//trim(fname)//', which ATLAS GPB writes; remove the key'
      error stop 1
    enddo

    inquire(file=fname, exist=exists)
    if (.not. exists) then
      do m = 1, nmat
        if (p%process(m) == PROC_NONE) cycle
        write(*,'(A)') E//trim(SOLVER(which))//' '//Process_Verb(m)//' material '//trim(adjustl(material_names(m)))// &
          ' and '//trim(fname)//' is missing: run ATLAS GPB, which writes it from the [GPB-Phase*] sections'
        error stop 1
      enddo
      call Store()
      return
    endif

    ! --- the records ---------------------------------------------------------------------------
    gasline = .false.
    open(newunit=u, file=fname, status='old', action='read')
    do
      call Next_Record(u, line, ios); if (ios /= 0) exit
      read(line,*) key
      select case (trim(key))
      case ('gas')
        read(line,*,iostat=ios) key, gpref, nsfile
        if (ios == 0 .and. trim(gpref) == '-') gpref = ''
        if (ios /= 0 .or. trim(gpref) /= trim(MOSE_phase_prefix) .or. nsfile /= ns) then
          write(*,'(A)') E//trim(fname)//' was written for another gas phase: run ATLAS GPB again'
          error stop 1
        endif
        gasline = .true.
      case ('material')
        read(line,*,iostat=ios) key, mname, proc, nterm, enword, h298
        if (ios /= 0) call Bad('a record that does not read "material <name> <process> <terms> <energy> <h298>"')
        pname = ''; cprod = 0._R8; hprod = 0._R8
        if (trim(proc) == 'cap') then
          read(line,*,iostat=ios) key, mname, proc, nterm, enword, h298, pname, cprod, hprod
          if (ios /= 0) call Bad('a cap record that does not end with "<product> <kg per kg> <its h298>"')
        endif
        ! a combustion or release record whose reaction makes a condensed product: the smoke, the ICE phase it is born in
        sname = ''; snu = 0._R8; sh = 0._R8; sphase = ''
        if (trim(proc) == 'combustion' .or. trim(proc) == 'release') then
          read(line,*,iostat=ios) key, mname, proc, nterm, enword, h298, sname
          if (ios == 0) then
            read(line,*,iostat=ios) key, mname, proc, nterm, enword, h298, sname, snu, sh, sphase
            if (ios /= 0) call Bad('a '//trim(proc)//' record whose smoke does not read "<material> <kg per kg> '// &
                                   '<its h298> <ICE phase>"')
            if (.not. snu > 0._R8) call Bad('a '//trim(proc)//' record whose smoke forms no mass')
          else
            sname = ''
          endif
        endif
        im = 0
        do m = 1, nmat
          if (trim(adjustl(material_names(m))) == trim(mname)) im = m
        enddo
        if (im == 0) then
          write(*,'(A)') E//trim(fname)//' names the material '//trim(mname)//', which '//trim(SOLVER(which))// &
            ' does not have: run ATLAS GPB again'
          error stop 1
        endif
        listed(im) = .true.
        select case (trim(enword))
        case ('none');      en = EN_NONE
        case ('latent');    en = EN_LATENT
        case ('formation'); en = EN_FORMATION
        case default;       en = -1
        end select
        ! row: 1 the process the material runs, 2 its cap, 3 the release of the oxide it is, 0 a record not used
        row = 0
        select case (trim(proc))
        case ('none');        if (nterm /= 0 .or. en /= EN_NONE) call Bad('a record "none" with terms or an energy word')
        case ('evaporation'); if (nterm < 1 .or. en /= EN_LATENT) call Bad('an evaporation record that is not "latent"')
                              if (p%process(im) == PROC_EVAP) row = 1
        case ('combustion');  if ((nterm < 1 .and. len_trim(sname) == 0) .or. en /= EN_FORMATION) &
                                call Bad('a combustion record that is not "formation"')
                              if (p%process(im) == PROC_COMB) row = 1
        case ('cap');         if (nterm < 0 .or. en /= EN_FORMATION .or. .not. cprod > 0._R8) &
                                call Bad('a cap record that is not "formation" or forms no product')
                              if (p%prod(im) > 0) row = 2
        case ('release');     if ((nterm < 1 .and. len_trim(sname) == 0) .or. en /= EN_FORMATION) &
                                call Bad('a release record that is not "formation"')
                              if (relwanted(im)) row = 3
        case default;         call Bad('the process "'//trim(proc)//'"')
        end select
        if ((row == 1 .and. found(im)) .or. (row == 2 .and. capfound(im)) .or. (row == 3 .and. relfound(im))) &
          call Bad('two records of the same process for '//trim(mname))
        do k = 1, nterm
          call Next_Record(u, line, ios)
          if (ios /= 0) call Bad('fewer terms than its last record announces')
          read(line,*,iostat=ios) spname, cmass, cmole
          if (ios /= 0) call Bad('a term that does not read "<species> <nu_mass> <nu_mole>"')
          ! every gas species of a record is a species of the gas phase, whether this run uses the record or not
          s = Species_Index(spname)
          if (s == 0) then
            write(*,'(A)') E//trim(fname)//' names the gas species '//trim(spname)//', which the gas phase does not '// &
              'have: run ATLAS GPB again'
            error stop 1
          endif
          select case (row)
          case (1); nu_s(im,s)  = nu_s(im,s)  + cmass
          case (2); cap_s(im,s) = cap_s(im,s) + cmass
          case (3); rel_s(im,s) = rel_s(im,s) + cmass
          end select
        enddo
        select case (row)
        case (1)
          found(im) = .true.; p%energy(im) = en; p%h298(im) = h298
          p%smoke(im) = sname; smnu(im) = snu; p%hsm(im) = sh; p%sphase(im) = sphase
        case (2)
          capfound(im) = .true.; capnu(im) = cprod; p%hprod(im) = hprod; capname(im) = pname
        case (3)
          relfound(im) = .true.
          p%rsmoke(im) = sname; p%rnusm(im) = snu; p%rhsm(im) = sh; p%rsphase(im) = sphase
        end select
      case default
        call Bad('the record "'//trim(key)//'"')
      end select
    enddo
    close(u)
    if (.not. gasline) call Bad('no "gas" line')

    call Check_Process_Rows()
    call Check_Cap_Rows()
    call Store()

  contains

    !> The book of the table, emptied at the size of this solver's materials.
    subroutine Reset_Book()
      if (allocated(p%process)) deallocate(p%process, p%energy, p%h298, p%prod, p%hprod)
      if (allocated(p%smoke))   deallocate(p%smoke, p%sphase, p%hsm)
      if (allocated(p%rsmoke))  deallocate(p%rsmoke, p%rsphase, p%rnusm, p%rhsm)
      allocate(p%process(n1), p%energy(n1), p%prod(n1))
      p%process = PROC_NONE; p%energy = EN_NONE; p%prod = 0
      allocate(p%h298(n1), p%hprod(n1), p%hsm(n1), p%rnusm(n1), p%rhsm(n1), source=0._R8)
      allocate(p%smoke(n1), p%sphase(n1), p%rsmoke(n1), p%rsphase(n1))
      p%smoke = ''; p%sphase = ''; p%rsmoke = ''; p%rsphase = ''
    end subroutine Reset_Book

    !> Every material is listed; one that converts has the record of its process, whose gas species and smoke
    !  together take the mass the material loses.
    subroutine Check_Process_Rows()
      do m = 1, nmat
        if (.not. listed(m)) then
          write(*,'(A)') E//'material '//trim(adjustl(material_names(m)))//' of '//trim(SOLVER(which))// &
            ' is missing from '//trim(fname)//': run ATLAS GPB again'
          error stop 1
        endif
        if (p%process(m) == PROC_NONE) cycle
        if (.not. found(m)) then
          write(*,'(A)') E//trim(SOLVER(which))//' '//Process_Verb(m)//' material '// &
            trim(adjustl(material_names(m)))//' and '//trim(fname)//' has no '// &
            trim(merge('combustion ', 'evaporation', p%process(m) == PROC_COMB))//' record for it: name its vapour '// &
            '(psat-vapour) or its reaction (particle-mechanism) in its [GPB-Phase*] section and run ATLAS GPB again'
          error stop 1
        endif
        rowsum = sum(nu_s(m,1:ns)) + smnu(m)
        if (abs(rowsum - 1._R8) > 1.e-8_R8) then
          write(*,'(A,ES23.15)') E//trim(fname)//': the row of '//trim(adjustl(material_names(m)))// &
            ' with its smoke does not sum to 1 (the gas and the smoke would not gain the mass the material loses): ', &
            rowsum
          error stop 1
        endif
      enddo
    end subroutine Check_Process_Rows

    !> A material that burns with an oxide cap: its cap record, the oxide it names, the ratio the solver's law holds,
    !  and the release record of the oxide when the cap is handed to the gas.
    subroutine Check_Cap_Rows()
      do m = 1, nmat
        if (p%prod(m) == 0) cycle
        if (.not. capfound(m)) then
          write(*,'(A)') E//trim(SOLVER(which))//' material '//trim(adjustl(material_names(m)))//' burns with an '// &
            'oxide cap and '//trim(fname)//' has no cap record for it: add the reaction that makes '// &
            trim(adjustl(material_names(p%prod(m))))//' from it to particle-mechanism and run ATLAS GPB again'
          error stop 1
        endif
        if (trim(capname(m)) /= trim(adjustl(material_names(p%prod(m))))) then
          write(*,'(A)') E//trim(fname)//': the cap record of '//trim(adjustl(material_names(m)))//' makes '// &
            trim(capname(m))//', and '//trim(SOLVER(which))//' holds its product as '// &
            trim(adjustl(material_names(p%prod(m))))//': run ATLAS GPB again'
          error stop 1
        endif
        rowsum = sum(cap_s(m,1:ns)) + capnu(m)
        if (abs(rowsum - 1._R8) > 1.e-8_R8) then
          write(*,'(A,ES23.15)') E//trim(fname)//': the cap record of '//trim(adjustl(material_names(m)))// &
            ' does not close the mass (its gas species and its oxide per kg of the material): ', rowsum
          error stop 1
        endif
        if (present(oxide_ratio)) then
          ! 1e-7: one formula with two tables of atomic weights that differ in their last digit stays within it
          if (abs(capnu(m) - oxide_ratio) > 1.e-7_R8*oxide_ratio) then
            write(*,'(A,ES15.7,A,ES15.7,A)') E//trim(fname)//': the cap record of '// &
              trim(adjustl(material_names(m)))//' forms ', capnu(m), ' kg of oxide per kg, and the cap law of '// &
              trim(SOLVER(which))//' holds ', oxide_ratio, ' (alumina from aluminium): the gas would not receive '// &
              'what the cap does not take'
            error stop 1
          endif
        endif
        if (relwanted(p%prod(m))) then
          if (.not. relfound(p%prod(m))) then
            write(*,'(A)') E//trim(SOLVER(which))//' material '//trim(adjustl(material_names(m)))//' hands its cap '// &
              'to the gas at burnout (residue = off) and '//trim(fname)//' has no release record for '// &
              trim(adjustl(material_names(p%prod(m))))//': add the reaction that turns it into gas species to '// &
              'particle-mechanism and run ATLAS GPB again, or keep the residue (residue = on)'
            error stop 1
          endif
          rowsum = sum(rel_s(p%prod(m),1:ns)) + p%rnusm(p%prod(m))
          if (abs(rowsum - 1._R8) > 1.e-8_R8) then
            write(*,'(A,ES23.15)') E//trim(fname)//': the release record of '// &
              trim(adjustl(material_names(p%prod(m))))//' does not sum to 1: ', rowsum
            error stop 1
          endif
        endif
      enddo
    end subroutine Check_Cap_Rows

    !> The rows the coupler reads, from those read: the process rows, the cap of a capped material per kg of its
    !  oxide, the release row of that oxide, the smoke per kg burnt.
    subroutine Store()
      integer :: mm
      if (allocated(p%nu)) deallocate(p%nu, p%mapped, p%capped, p%nuCap, p%nuRel, p%nuOxide, p%nuSmoke, p%dhSmoke)
      allocate(p%nu, source=nu_s)
      allocate(p%mapped(n1)); p%mapped = p%process /= PROC_NONE
      p%active = any(p%process(1:nmat) /= PROC_NONE)
      allocate(p%capped(n1)); p%capped = .false.
      allocate(p%nuCap(n1,ns1), p%nuRel(n1,ns1), source=0._R8)
      allocate(p%nuOxide(n1)); p%nuOxide = 1._R8
      do mm = 1, nmat
        if (p%prod(mm) == 0 .or. .not. capfound(mm)) cycle
        p%capped(mm)  = .true.
        p%nuOxide(mm) = capnu(mm)
        p%nuCap(mm,:) = cap_s(mm,:) / capnu(mm)
        if (relfound(p%prod(mm))) p%nuRel(mm,:) = rel_s(p%prod(mm),:)
      enddo
      allocate(p%nuSmoke, source=smnu)
      allocate(p%dhSmoke(n1)); p%dhSmoke = 0._R8
      if (.not. p%active) &
        write(*,'(A)') ' [MI2] '//trim(SOLVER(which))//': no material evaporates or burns: no mass is handed to the gas.'
    end subroutine Store

    !> 'burns' or 'evaporates', padded alike, for material m.
    function Process_Verb(m) result(w)
      integer, intent(in) :: m
      character(len=:), allocatable :: w
      w = trim(merge('burns     ', 'evaporates', p%process(m) == PROC_COMB))
    end function Process_Verb

    subroutine Bad(what)
      character(len=*), intent(in) :: what
      write(*,'(A)') E//trim(fname)//' holds '//what//': run ATLAS GPB again'
      error stop 1
    end subroutine Bad

    !> The next line that is neither blank nor a comment.
    subroutine Next_Record(unit, buf, stat)
      integer,          intent(in)  :: unit
      character(len=*), intent(out) :: buf
      integer,          intent(out) :: stat
      do
        read(unit,'(A)',iostat=stat) buf; if (stat /= 0) return
        buf = adjustl(buf)
        if (len_trim(buf) > 0 .and. buf(1:1) /= '#') return
      enddo
    end subroutine Next_Record

    integer function Species_Index(name)
      use FLINT_Lib_Thermodynamic, only: species_names
      character(len=*), intent(in) :: name
      integer :: i
      Species_Index = 0
      do i = 1, ns
        if (trim(adjustl(species_names(i))) == trim(name)) then
          Species_Index = i; return
        endif
      enddo
    end function Species_Index

  end subroutine Load_Species_Map


  !> The smoke a record of material m of a solver makes: its combustion record or, release, the release record of m as
  !  a cap oxide handed to the gas at burnout. False without one; else the smoke's material, kg per kg, formation
  !  enthalpy at T_datum_ref and the ICE phase (prefix) it is born in.
  logical function Smoke_Record ( which, m, name, nus, h298, phase, release )
    integer,          intent(in)  :: which, m
    character(len=*), intent(out) :: name, phase
    real(R8),         intent(out) :: nus, h298
    logical, optional, intent(in) :: release
    type(gas_products_t), pointer :: p
    logical :: rel

    name = ''; phase = ''; nus = 0._R8; h298 = 0._R8
    Smoke_Record = .false.
    rel = .false.
    if (present(release)) rel = release
    p => Table(which)
    if (.not. allocated(p%smoke)) return
    if (m < 1 .or. m > size(p%smoke)) return
    if (rel) then
      if (len_trim(p%rsmoke(m)) == 0) return
      name = p%rsmoke(m); phase = p%rsphase(m); nus = p%rnusm(m); h298 = p%rhsm(m)
    else
      if (len_trim(p%smoke(m)) == 0) return
      name = p%smoke(m); phase = p%sphase(m); nus = p%nuSmoke(m); h298 = p%hsm(m)
    endif
    Smoke_Record = .true.

  end function Smoke_Record


  !> The datum offset of the smoke of ICE material m: its formation enthalpy less ICE's energy of it at T_datum_ref.
  subroutine Set_Smoke_Datum ( m, dh )
    integer,  intent(in) :: m
    real(R8), intent(in) :: dh

    productsICE%dhSmoke(m) = dh

  end subroutine Set_Smoke_Datum


  !> The rates of a capped IGLOO material [kg/s] from its net mass column and its cap's two: (1) the material
  !  consumed, (2) the part of it burnt into gas species, (3) the oxide its cap gained, (4) the cap oxide handed
  !  to the gas.
  pure function Cap_Rates ( m, net, gained, handed ) result(cr)
    integer,  intent(in) :: m
    real(R8), intent(in) :: net, gained, handed
    real(R8) :: cr(4)

    cr(3) = gained
    cr(4) = handed
    cr(1) = net + gained - handed
    cr(2) = cr(1) - gained/productsIGLOO%nuOxide(m)

  end function Cap_Rates


  !> Echo the rows a solver's materials are handed to the gas with, and the process of each.
  subroutine Report_Stoichiometry ( which, material_names, nmat )
    use FLINT_Lib_Thermodynamic, only: ns, species_names
    implicit none
    integer,          intent(in) :: which
    character(len=*), intent(in) :: material_names(:)
    integer,          intent(in) :: nmat
    character(len=11), parameter :: PROCW(0:2) = ['inert      ', 'evaporation', 'combustion ']
    type(gas_products_t), pointer :: p
    integer  :: m, s
    real(R8) :: rowsum

    p => Table(which)
    if (.not. allocated(p%process)) return
    if (all(p%process(1:nmat) == PROC_NONE)) return

    write(*,*)
    write(*,'(A)') ' [MI2] '//trim(SOLVER(which))//' material -> gas species [kg of gas per kg of material converted], '// &
                   trim(p%file)
    do m = 1, nmat
      if (p%process(m) == PROC_NONE) then
        write(*,'(A,A,A)') '   - ', trim(material_names(m)), ' : (inert: no mass handed to the gas)'
        cycle
      endif
      write(*,'(A,A,A,A,A)') '   - ', trim(material_names(m)), ' [', trim(PROCW(p%process(m))), '] :'
      rowsum = 0._R8
      do s = 1, ns
        if (p%nu(m,s) == 0._R8) cycle
        call Put_Species( s, p%nu(m,s) )
        rowsum = rowsum + p%nu(m,s)
      enddo
      write(*,'(A,F10.5)') '        SUM : ', rowsum
      if (.not. p%capped(m)) cycle
      write(*,'(A,F10.5,A)') '        its oxide cap takes '//trim(material_names(p%prod(m)))//' : ', &
        p%nuOxide(m), '  kg per kg of the material; per kg of oxide the gas exchanges'
      do s = 1, ns
        if (p%nuCap(m,s) /= 0._R8) call Put_Species( s, p%nuCap(m,s) )
      enddo
      if (all(p%nuRel(m,1:ns) == 0._R8)) then
        write(*,'(A)') '        the cap stays condensed at burnout (residue = on)'
        cycle
      endif
      write(*,'(A)') '        the cap is handed to the gas at burnout (residue = off), per kg of it'
      do s = 1, ns
        if (p%nuRel(m,s) /= 0._R8) call Put_Species( s, p%nuRel(m,s) )
      enddo
    enddo
    write(*,*)

  contains

    !> One species of a row: its name, kg per kg, produced or consumed.
    subroutine Put_Species ( s, x )
      integer,  intent(in) :: s
      real(R8), intent(in) :: x
      write(*,'(A,A12,A,F10.5,A)') '        ', trim(adjustl(species_names(s))), ' : ', x, &
        merge('  (consumed)', '  (produced)', x < 0._R8)
    end subroutine Put_Species

  end subroutine Report_Stoichiometry


  !> Build the datum offset of each converting material of solver which from the energy rule of its record, href the
  !  enthalpy its solver credits at T_datum_ref on its own datum:
  !    latent    : SUM_s nu*h_FLINT(s) - (href + Lv)   the vapour's enthalpy less the latent heat is the liquid's
  !    formation : h298 - href                         the record gives the material's formation-inclusive enthalpy
  !  and the same of the oxide of a capped material. Thot is the temperature the material converts at (boiling,
  !  burning) and hhot the enthalpy its solver credits there: they size two messages, never an offset.
  !  [MI2-Coupling] datum-correction = false zeroes the offsets.
  subroutine Build_Datum_Offsets ( which, material_names, nmat, beck, beta, href, Lv, qcomb, Thot, hhot, Lvconst )
    use FLINT_Lib_Thermodynamic, only: ns, species_names, H0i
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    implicit none
    integer,          intent(in) :: which
    character(len=*), intent(in) :: material_names(:)
    integer,          intent(in) :: nmat
    logical,          intent(in) :: beck(:)    !> the material burns on Beckstead's law: q-comb is a share of its heat
    real(R8),         intent(in) :: beta(:)    !> beta-part, the share of q-comb the drop keeps (Beckstead) [-]
    real(R8),         intent(in) :: href(:)    !> enthalpy the solver credits at T_ref, its own datum [J/kg]
    real(R8),         intent(in) :: Lv(:)      !> latent heat at T_ref [J/kg]
    real(R8),         intent(in) :: qcomb(:)   !> the solver's q-comb [J/kg]
    real(R8),         intent(in) :: Thot(:)    !> boiling temperature (evaporation), burning temperature (combustion) [K]
    real(R8),         intent(in) :: hhot(:)    !> enthalpy the solver credits at Thot, its own datum [J/kg]
    logical,          intent(in) :: Lvconst(:) !> the solver holds the latent heat as a constant, not as a column
    type(gas_products_t), pointer :: p
    integer  :: m, s, ip
    real(R8) :: hprod, hcredited, qr, hph, off, qhot
    real(R8), allocatable :: dh(:), dhp(:)
    logical  :: apply

    p => Table(which)
    allocate(dh(1:max(nmat,1)), dhp(1:max(nmat,1)), source=0._R8)

    if (p%active) then
      ! datum-correction = false computes the offsets but does not apply them
      apply = hydra_MI2_obj%datum_correction

      write(*,*)
      write(*,'(A)') ' [MI2] enthalpy-datum reconciliation  [J per kg of material converted]'
      write(*,'(A,F8.2,A)') '       anchored at T_ref = ', T_datum_ref, ' K (FLINT standard state)'
      do m = 1, nmat
        if ( p%process(m) == PROC_NONE ) cycle

        hprod = 0._R8; hph = 0._R8
        do s = 1, ns
          if (p%nu(m,s) == 0._R8) cycle
          hprod = hprod + p%nu(m,s) * H0i(s, T_datum_ref)
          hph   = hph   + p%nu(m,s) * H0i(s, Thot(m))
        enddo

        write(*,'(A,A)')          '   - ', trim(material_names(m))
        do s = 1, ns
          if (p%nu(m,s) == 0._R8) cycle
          write(*,'(A,A12,A,ES15.7)') '        h_FLINT  ', trim(adjustl(species_names(s))), ' : ', &
                                      H0i(s, T_datum_ref)
        enddo
        write(*,'(A,ES15.7)')     '        products (SUM nu*h_FLINT)          : ', hprod

        select case (p%energy(m))
        case (EN_LATENT)
          hcredited = href(m) + Lv(m)
          dh(m)     = hprod - hcredited
          write(*,'(A,ES15.7)')   '        h_'//SOLVER(which)//'(liquid table)              : ', href(m)
          write(*,'(A,ES15.7)')   '        Lv (already inside the source)     : ', Lv(m)
          write(*,'(A,ES15.7)')   '        credited by '//SOLVER(which)//' (h + Lv)         : ', hcredited
          ! a constant latent heat: what it misses of the one the two heat capacities give, at the boiling temperature
          if (Lvconst(m)) then
            off = (hhot(m) - href(m)) - (hph - hprod)
            write(*,'(A,F8.2,A,ES12.4,A,F6.2,A)') '        [WARNING] the latent heat is a constant: at ', Thot(m), &
              ' K it differs by ', off, ' J/kg (', 100._R8*abs(off)/max(abs(Lv(m)), tiny(1._R8)), ' % of Lv) from '// &
              'the one the heat capacities of the liquid and of the vapour give. The energy the gas and the drops '// &
              'hold together is exact; its share between them is not: give the latent heat as the Lv column of '// &
              'the property table, which ATLAS writes'
          endif
        case (EN_FORMATION)
          dh(m) = p%h298(m) - href(m)
          qr    = p%h298(m) - hprod
          write(*,'(A,ES15.7)')   '        h298 of the material (ATLAS)       : ', p%h298(m)
          write(*,'(A,ES15.7)')   '        h_'//SOLVER(which)//'(table)                     : ', href(m)
          write(*,'(A,ES15.7)')   '        heat of reaction, in the species   : ', qr
          ! like with like: from the condensed material where it burns to the record's species there, beside q-comb
          qhot  = p%h298(m) + (hhot(m) - href(m)) - hph
          write(*,'(A,F8.2,A,ES15.7)') '        the same from the material at ', Thot(m), ' K, where it burns : ', qhot
          write(*,'(A,ES15.7)')   '        q-comb of the solver               : ', qcomb(m)
          if (beck(m)) then
            if ((qcomb(m) > 0._R8 .or. beta(m) > 0._R8) .and. abs(qcomb(m) - qhot) > 0.05_R8*abs(qhot)) then
              write(*,'(A,ES15.7,A)') '        [WARNING] q-comb = ', qcomb(m), ' differs from it by more than 5 %: '// &
                'the share beta-part*q-comb the drop keeps is taken from another heat than the gas receives'
            endif
          else
            write(*,'(A)')        '        (combustion = direct: q-comb sets the rate; the gas receives the heat above)'
            if (qcomb(m) == 0._R8) write(*,'(A)') '        (q-comb = 0: the rate is pure vaporisation, '// &
              'B = cp (T_g - T_p)/L_v; it ignores the heat of reaction the gas receives)'
          endif
          ! the oxide of a capped material, on the same two datums
          if (p%capped(m)) then
            ip     = p%prod(m)
            dhp(m) = p%hprod(m) - href(ip)
            qr     = p%h298(m) - p%nuOxide(m)*p%hprod(m)
            do s = 1, ns
              qr = qr - p%nuOxide(m)*p%nuCap(m,s) * H0i(s, T_datum_ref)
            enddo
            write(*,'(A,ES15.7)') '        h298 of its oxide (ATLAS)          : ', p%hprod(m)
            write(*,'(A,ES15.7)') '        h_'//SOLVER(which)//'(table) of its oxide        : ', href(ip)
            write(*,'(A,ES15.7)') '        heat of the cap reaction, per kg   : ', qr
            if (.not. apply) dhp(m) = 0._R8
            write(*,'(A,ES15.7)') '        => datum offset of the oxide       : ', dhp(m)
          endif
        end select
        if (.not. apply) dh(m) = 0._R8
        write(*,'(A,ES15.7)')     '        => datum offset applied by MI2     : ', dh(m)
        if ( apply .and. abs(dh(m)) < 1.e-6_R8*max(abs(hprod),1._R8) ) then
          write(*,'(A)')          '        (offset ~ 0: the two datums already agree)'
        endif
      enddo
      if (.not. apply) then
        write(*,'(A)') '   [WARNING] datum-correction = false: the offsets above are computed but NOT'
        write(*,'(A)') '             applied. The gas receives energy on the condensed solver''s datum,'
        write(*,'(A)') '             which is KNOWN WRONG. Diagnostic use only.'
      endif
      write(*,*)
    endif

    if (allocated(p%dh)) deallocate(p%dh, p%dhOxide)
    allocate(p%dh, source=dh)
    allocate(p%dhOxide, source=dhp)

  end subroutine Build_Datum_Offsets

end module hydraMI2_Stoichiometry
