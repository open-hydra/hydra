module hydraMI2_Stoichiometry
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Load_Species_Map, Report_Stoichiometry, Build_Datum_Offsets
  public :: nu, has_stoichiometry, dhDatum, T_datum_ref
  public :: nuICE, has_stoichiometry_ICE, dhDatumICE, mappedICE
  public :: capped, nuCap, nuRel, nuOxide, dhProduct
  public :: Smoke_Record, Set_Smoke_Datum, nuSmokeICE, dhSmokeICE
  public :: MAP_IGLOO, MAP_ICE

  !> The two condensed solvers a map is loaded for.
  integer, parameter :: MAP_IGLOO = 1, MAP_ICE = 2

  !> nu(m,s) = kg of gas species s exchanged per kg of IGLOO material m converted (nu > 0 produced,
  !  nu < 0 consumed); a row sums to 1. nuICE: the same for ICE's materials.
  real(R8), allocatable, protected :: nu(:,:), nuICE(:,:)
  logical,               protected :: has_stoichiometry = .false., has_stoichiometry_ICE = .false.

  !> Reference temperature of the datum offset (FLINT standard state).
  real(R8), parameter :: T_datum_ref = 298.15_R8

  !> dhDatum(m) = constant enthalpy-datum offset of IGLOO material m [J per kg converted]: its enthalpy on the
  !  gas's formation-inclusive datum less the one IGLOO credits, at T_datum_ref. dhDatumICE: ICE's materials.
  real(R8), allocatable, protected :: dhDatum(:), dhDatumICE(:)

  !> mappedICE(m): ICE material m has a record (it evaporates or burns).
  logical,  allocatable, protected :: mappedICE(:)

  !> An IGLOO material that burns with an oxide cap: capped(m). Per kg of oxide its cap gains the gas exchanges
  !  nuCap(m,s) (the row sums to 1/nuOxide(m) - 1); per kg of cap oxide handed to the gas at burnout it receives
  !  nuRel(m,s) (the row sums to 1). nuOxide(m): kg of oxide per kg of the material its cap takes.
  !  dhProduct(m): the datum offset of the oxide, as dhDatum.
  logical,  allocatable, protected :: capped(:)
  real(R8), allocatable, protected :: nuCap(:,:), nuRel(:,:), nuOxide(:), dhProduct(:)
  !> The smoke of an ICE material that burns into a condensed product (the tail of its combustion record): kg of
  !  smoke per kg burnt, and the smoke material's datum offset (its formation enthalpy at T_datum_ref less ICE's
  !  energy of it there), 0 without a smoke.
  real(R8), allocatable, protected :: nuSmokeICE(:), dhSmokeICE(:)

  !> Process a material runs, and the energy rule of its record.
  integer, parameter :: PROC_NONE = 0, PROC_EVAP = 1, PROC_COMB = 2
  integer, parameter :: EN_NONE = 0, EN_LATENT = 1, EN_FORMATION = 2

  !> What the products file says of one solver's materials.
  type :: map_book_t
    character(len=256)    :: file = ''
    integer,  allocatable :: process(:), energy(:)
    real(R8), allocatable :: h298(:)      !> the material's formation-inclusive enthalpy at T_datum_ref [J/kg]
    integer,  allocatable :: prod(:)      !> the material a capped material's oxide is, else 0
    real(R8), allocatable :: hprod(:)     !> that oxide's formation-inclusive enthalpy at T_datum_ref [J/kg]
    character(len=64), allocatable :: smoke(:), sphase(:)   !> the smoke a burning material makes, the phase it is of
    real(R8), allocatable :: nusm(:), hsm(:)                !> kg of smoke per kg burnt, its h298 [J/kg]
  end type map_book_t
  type(map_book_t) :: book(2)

  character(len=5), parameter :: SOLVER(2) = ['IGLOO', 'ICE  ']

contains

  !> Read INPUT/<prefix>products.txt (ATLAS GPB) and keep, for each material of the solver, the record of the
  !  process it runs: combustion if it burns, else evaporation if it evaporates; for a material that burns with an
  !  oxide cap (product > 0) also its cap record and, when its cap is handed to the gas (releases), the release
  !  record of its oxide. Stops on a material that converts without a record, a file of another gas phase, an
  !  unknown material or species, a row that does not close the mass.
  subroutine Load_Species_Map ( which, prefix, material_names, evaporates, burns, product, releases, oxide_ratio )
    use Finer,                   only: file_ini
    use FLINT_Lib_Thermodynamic, only: ns, species_names
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
    character(len=512)  :: fname, line
    character(len=1024) :: raw
    character(len=64)   :: key, mname, proc, enword, spname, gpref, pname, sname, sphase
    integer  :: nmat, u, ios, m, s, k, nterm, nsfile, im, en, error, row
    real(R8) :: h298, cmass, cmole, rowsum, cprod, hprod, snu, sh
    logical  :: exists, gasline
    logical,  allocatable :: listed(:), found(:), capfound(:), relfound(:), relwanted(:)
    integer,  allocatable :: prod(:)
    real(R8), allocatable :: nu_s(:,:), cap_s(:,:), rel_s(:,:), capnu(:), caph(:), smnu(:), smh(:)
    character(len=64), allocatable :: capname(:), smn(:), smph(:)
    type(file_ini) :: fini

    nmat = size(material_names)
    allocate(nu_s(max(nmat,1), max(ns,1)), listed(max(nmat,1)), found(max(nmat,1)))
    nu_s = 0._R8; listed = .false.; found = .false.
    allocate(cap_s(max(nmat,1), max(ns,1)), rel_s(max(nmat,1), max(ns,1)), capnu(max(nmat,1)), caph(max(nmat,1)), &
             capname(max(nmat,1)), capfound(max(nmat,1)), relfound(max(nmat,1)), relwanted(max(nmat,1)), &
             prod(max(nmat,1)))
    cap_s = 0._R8; rel_s = 0._R8; capnu = 0._R8; caph = 0._R8; capname = ''
    allocate(smnu(max(nmat,1)), smh(max(nmat,1)), smn(max(nmat,1)), smph(max(nmat,1)))
    smnu = 0._R8; smh = 0._R8; smn = ''; smph = ''
    capfound = .false.; relfound = .false.; relwanted = .false.; prod = 0
    if (present(product)) prod(1:nmat) = product(1:nmat)
    if (present(releases)) then
      do m = 1, nmat
        if (prod(m) > 0 .and. releases(m)) relwanted(prod(m)) = .true.
      enddo
    endif
    associate ( bk => book(which) )
    if (allocated(bk%process)) deallocate(bk%process, bk%energy, bk%h298, bk%prod, bk%hprod)
    if (allocated(bk%smoke)) deallocate(bk%smoke, bk%sphase, bk%nusm, bk%hsm)
    allocate(bk%process(max(nmat,1)), bk%energy(max(nmat,1)), bk%h298(max(nmat,1)), bk%prod(max(nmat,1)), &
             bk%hprod(max(nmat,1)))
    allocate(bk%smoke(max(nmat,1)), bk%sphase(max(nmat,1)), bk%nusm(max(nmat,1)), bk%hsm(max(nmat,1)))
    bk%smoke = ''; bk%sphase = ''; bk%nusm = 0._R8; bk%hsm = 0._R8
    bk%process = PROC_NONE; bk%energy = EN_NONE; bk%h298 = 0._R8; bk%prod = prod; bk%hprod = 0._R8
    do m = 1, nmat
      if (burns(m)) then
        bk%process(m) = PROC_COMB
      elseif (evaporates(m)) then
        bk%process(m) = PROC_EVAP
      endif
    enddo
    fname = 'INPUT/'//trim(prefix)//'products.txt'
    bk%file = fname

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
        if (bk%process(m) == PROC_NONE) cycle
        write(*,'(A)') E//trim(SOLVER(which))//' '//trim(merge('burns     ', 'evaporates', bk%process(m) == PROC_COMB))// &
          ' material '//trim(adjustl(material_names(m)))//' and '//trim(fname)//' is missing: run ATLAS GPB, '// &
          'which writes it from the [GPB-Phase*] sections'
        error stop 1
      enddo
      call Store()
      return
    endif

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
        ! a combustion record whose reaction makes a condensed product: the smoke, the ICE phase it is born in
        sname = ''; snu = 0._R8; sh = 0._R8; sphase = ''
        if (trim(proc) == 'combustion') then
          read(line,*,iostat=ios) key, mname, proc, nterm, enword, h298, sname
          if (ios == 0) then
            read(line,*,iostat=ios) key, mname, proc, nterm, enword, h298, sname, snu, sh, sphase
            if (ios /= 0) call Bad('a combustion record whose smoke does not read "<material> <kg per kg> <its h298> '// &
                                   '<ICE phase>"')
            if (.not. snu > 0._R8) call Bad('a combustion record whose smoke forms no mass')
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
                              if (bk%process(im) == PROC_EVAP) row = 1
        case ('combustion');  if (nterm < 1 .or. en /= EN_FORMATION) &
                                call Bad('a combustion record that is not "formation"')
                              if (bk%process(im) == PROC_COMB) row = 1
        case ('cap');         if (nterm < 0 .or. en /= EN_FORMATION .or. .not. cprod > 0._R8) &
                                call Bad('a cap record that is not "formation" or forms no product')
                              if (prod(im) > 0) row = 2
        case ('release');     if (nterm < 1 .or. en /= EN_FORMATION) call Bad('a release record that is not "formation"')
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
          found(im) = .true.; bk%energy(im) = en; bk%h298(im) = h298
          smn(im) = sname; smnu(im) = snu; smh(im) = sh; smph(im) = sphase
        case (2)
          capfound(im) = .true.; capnu(im) = cprod; caph(im) = hprod; capname(im) = pname
        case (3)
          relfound(im) = .true.
        end select
      case default
        call Bad('the record "'//trim(key)//'"')
      end select
    enddo
    close(u)
    if (.not. gasline) call Bad('no "gas" line')

    do m = 1, nmat
      if (.not. listed(m)) then
        write(*,'(A)') E//'material '//trim(adjustl(material_names(m)))//' of '//trim(SOLVER(which))// &
          ' is missing from '//trim(fname)//': run ATLAS GPB again'
        error stop 1
      endif
      if (bk%process(m) == PROC_NONE) cycle
      if (.not. found(m)) then
        write(*,'(A)') E//trim(SOLVER(which))//' '//trim(merge('burns     ', 'evaporates', bk%process(m) == PROC_COMB))// &
          ' material '//trim(adjustl(material_names(m)))//' and '//trim(fname)//' has no '// &
          trim(merge('combustion ', 'evaporation', bk%process(m) == PROC_COMB))//' record for it: name its vapour '// &
          '(psat-vapour) or its reaction (particle-mechanism) in its [GPB-Phase*] section and run ATLAS GPB again'
        error stop 1
      endif
      ! the gas's species and the smoke together take the mass the material loses
      rowsum = sum(nu_s(m,1:ns)) + smnu(m)
      if (abs(rowsum - 1._R8) > 1.e-8_R8) then
        write(*,'(A,ES23.15)') E//trim(fname)//': the row of '//trim(adjustl(material_names(m)))// &
          ' with its smoke does not sum to 1 (the gas and the smoke would not gain the mass the material loses): ', rowsum
        error stop 1
      endif
    enddo

    ! a material that burns with an oxide cap: its cap record, the oxide it names, the ratio the solver's law holds
    do m = 1, nmat
      if (prod(m) == 0) cycle
      if (.not. capfound(m)) then
        write(*,'(A)') E//trim(SOLVER(which))//' material '//trim(adjustl(material_names(m)))//' burns with an oxide '// &
          'cap and '//trim(fname)//' has no cap record for it: add the reaction that makes '// &
          trim(adjustl(material_names(prod(m))))//' from it to particle-mechanism and run ATLAS GPB again'
        error stop 1
      endif
      if (trim(capname(m)) /= trim(adjustl(material_names(prod(m))))) then
        write(*,'(A)') E//trim(fname)//': the cap record of '//trim(adjustl(material_names(m)))//' makes '// &
          trim(capname(m))//', and '//trim(SOLVER(which))//' holds its product as '// &
          trim(adjustl(material_names(prod(m))))//': run ATLAS GPB again'
        error stop 1
      endif
      rowsum = sum(cap_s(m,1:ns)) + capnu(m)
      if (abs(rowsum - 1._R8) > 1.e-8_R8) then
        write(*,'(A,ES23.15)') E//trim(fname)//': the cap record of '//trim(adjustl(material_names(m)))// &
          ' does not close the mass (its gas species and its oxide per kg of the material): ', rowsum
        error stop 1
      endif
      if (present(oxide_ratio)) then
        ! one formula, two tables of atomic weights that differ in their last digit (rev7/derive/oxide_ratio.py)
        if (abs(capnu(m) - oxide_ratio) > 1.e-7_R8*oxide_ratio) then
          write(*,'(A,ES15.7,A,ES15.7,A)') E//trim(fname)//': the cap record of '//trim(adjustl(material_names(m)))// &
            ' forms ', capnu(m), ' kg of oxide per kg, and the cap law of '//trim(SOLVER(which))//' holds ', &
            oxide_ratio, ' (alumina from aluminium): the gas would not receive what the cap does not take'
          error stop 1
        endif
      endif
      if (relwanted(prod(m))) then
        if (.not. relfound(prod(m))) then
          write(*,'(A)') E//trim(SOLVER(which))//' material '//trim(adjustl(material_names(m)))//' hands its cap to '// &
            'the gas at burnout (residue = off) and '//trim(fname)//' has no release record for '// &
            trim(adjustl(material_names(prod(m))))//': add the reaction that turns it into gas species to '// &
            'particle-mechanism and run ATLAS GPB again, or keep the residue (residue = on)'
          error stop 1
        endif
        rowsum = sum(rel_s(prod(m),1:ns))
        if (abs(rowsum - 1._R8) > 1.e-8_R8) then
          write(*,'(A,ES23.15)') E//trim(fname)//': the release record of '//trim(adjustl(material_names(prod(m))))// &
            ' does not sum to 1: ', rowsum
          error stop 1
        endif
      endif
      bk%hprod(m) = caph(m)
    enddo
    bk%smoke = smn; bk%sphase = smph; bk%nusm = smnu; bk%hsm = smh
    end associate
    call Store()

  contains

    !> Hand the rows to the solver's module arrays.
    subroutine Store()
      integer :: mm
      if (which == MAP_IGLOO) then
        if (allocated(nu)) deallocate(nu)
        allocate(nu, source=nu_s)
        has_stoichiometry = any(book(which)%process(1:nmat) /= PROC_NONE)
        if (allocated(capped)) deallocate(capped, nuCap, nuRel, nuOxide)
        allocate(capped(max(nmat,1)), nuCap(max(nmat,1), max(ns,1)), nuRel(max(nmat,1), max(ns,1)), nuOxide(max(nmat,1)))
        capped = .false.; nuCap = 0._R8; nuRel = 0._R8; nuOxide = 1._R8
        do mm = 1, nmat
          if (prod(mm) == 0 .or. .not. capfound(mm)) cycle
          capped(mm)  = .true.
          nuOxide(mm) = capnu(mm)
          nuCap(mm,:) = cap_s(mm,:) / capnu(mm)
          if (relfound(prod(mm))) nuRel(mm,:) = rel_s(prod(mm),:)
        enddo
      else
        if (allocated(nuICE)) deallocate(nuICE)
        allocate(nuICE, source=nu_s)
        has_stoichiometry_ICE = any(book(which)%process(1:nmat) /= PROC_NONE)
        if (allocated(mappedICE)) deallocate(mappedICE)
        allocate(mappedICE(max(nmat,1)))
        mappedICE = book(which)%process /= PROC_NONE
        if (allocated(nuSmokeICE)) deallocate(nuSmokeICE, dhSmokeICE)
        allocate(nuSmokeICE(max(nmat,1)), dhSmokeICE(max(nmat,1)))
        nuSmokeICE = smnu; dhSmokeICE = 0._R8
      endif
      if (.not. any(book(which)%process(1:nmat) /= PROC_NONE)) &
        write(*,'(A)') ' [MI2] '//trim(SOLVER(which))//': no material evaporates or burns: no mass is handed to the gas.'
    end subroutine Store

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


  !> The smoke the combustion record of material m of a solver makes: false without one; else its material, kg per kg
  !  burnt, formation enthalpy at T_datum_ref and the ICE phase (prefix) it is born in.
  logical function Smoke_Record ( which, m, name, nus, h298, phase )
    integer,          intent(in)  :: which, m
    character(len=*), intent(out) :: name, phase
    real(R8),         intent(out) :: nus, h298
    name = ''; phase = ''; nus = 0._R8; h298 = 0._R8
    Smoke_Record = .false.
    if (.not. allocated(book(which)%smoke)) return
    if (m > size(book(which)%smoke)) return
    if (len_trim(book(which)%smoke(m)) == 0) return
    name = book(which)%smoke(m); phase = book(which)%sphase(m); nus = book(which)%nusm(m); h298 = book(which)%hsm(m)
    Smoke_Record = .true.
  end function Smoke_Record


  !> The datum offset of the smoke of ICE material m: its formation enthalpy less ICE's energy of it at T_datum_ref.
  subroutine Set_Smoke_Datum ( m, dh )
    integer,  intent(in) :: m
    real(R8), intent(in) :: dh
    dhSmokeICE(m) = dh
  end subroutine Set_Smoke_Datum


  !> Echo the rows a solver's materials are handed to the gas with, and the process of each.
  subroutine Report_Stoichiometry ( which, material_names, nmat )
    use FLINT_Lib_Thermodynamic, only: ns, species_names
    implicit none
    integer,          intent(in) :: which
    character(len=*), intent(in) :: material_names(:)
    integer,          intent(in) :: nmat
    integer  :: m, s
    real(R8) :: rowsum
    character(len=11), parameter :: PROCW(0:2) = ['inert      ', 'evaporation', 'combustion ']

    if (.not. allocated(book(which)%process)) return
    if (all(book(which)%process(1:nmat) == PROC_NONE)) return

    write(*,*)
    write(*,'(A)') ' [MI2] '//trim(SOLVER(which))//' material -> gas species [kg of gas per kg of material converted], '// &
                   trim(book(which)%file)
    do m = 1, nmat
      if (book(which)%process(m) == PROC_NONE) then
        write(*,'(A,A,A)') '   - ', trim(material_names(m)), ' : (inert: no mass handed to the gas)'
        cycle
      endif
      write(*,'(A,A,A,A,A)') '   - ', trim(material_names(m)), ' [', trim(PROCW(book(which)%process(m))), '] :'
      rowsum = 0._R8
      do s = 1, ns
        if (Row(m,s) == 0._R8) cycle
        write(*,'(A,A12,A,F10.5,A)') '        ', trim(adjustl(species_names(s))), ' : ', Row(m,s), &
          merge('  (consumed)', '  (produced)', Row(m,s) < 0._R8)
        rowsum = rowsum + Row(m,s)
      enddo
      write(*,'(A,F10.5)') '        SUM : ', rowsum
      if (which /= MAP_IGLOO) cycle
      if (.not. capped(m)) cycle
      write(*,'(A,F10.5,A)') '        its oxide cap takes '//trim(material_names(book(which)%prod(m)))//' : ', &
        nuOxide(m), '  kg per kg of the material; per kg of oxide the gas exchanges'
      do s = 1, ns
        if (nuCap(m,s) == 0._R8) cycle
        write(*,'(A,A12,A,F10.5,A)') '        ', trim(adjustl(species_names(s))), ' : ', nuCap(m,s), &
          merge('  (consumed)', '  (produced)', nuCap(m,s) < 0._R8)
      enddo
      if (all(nuRel(m,1:ns) == 0._R8)) then
        write(*,'(A)') '        the cap stays condensed at burnout (residue = on)'
        cycle
      endif
      write(*,'(A)') '        the cap is handed to the gas at burnout (residue = off), per kg of it'
      do s = 1, ns
        if (nuRel(m,s) == 0._R8) cycle
        write(*,'(A,A12,A,F10.5,A)') '        ', trim(adjustl(species_names(s))), ' : ', nuRel(m,s), &
          merge('  (consumed)', '  (produced)', nuRel(m,s) < 0._R8)
      enddo
    enddo
    write(*,*)

  contains

    real(R8) function Row(m, s)
      integer, intent(in) :: m, s
      if (which == MAP_IGLOO) then
        Row = nu(m,s)
      else
        Row = nuICE(m,s)
      endif
    end function Row

  end subroutine Report_Stoichiometry


  !> Build the datum offset of each converting material from the energy rule of its record, href the enthalpy its
  !  solver credits at T_datum_ref on its own datum:
  !    latent    : SUM_s nu*h_FLINT(s) - (href + Lv)   the vapour's enthalpy less the latent heat is the liquid's
  !    formation : h298 - href                         the record gives the material's formation-inclusive enthalpy
  !  and the same of the oxide of a capped material. Thot is the temperature the material converts at (boiling,
  !  burning) and hhot the enthalpy its solver credits there: they size two messages, never an offset.
  !  [MI2-Coupling] datum-correction = false zeroes the offsets.
  subroutine Build_Datum_Offsets ( which, material_names, nmat, beck, beta, href, Lv, qcomb, Thot, hhot, Lvconst )
    use Finer,                   only: file_ini
    use FLINT_Lib_Thermodynamic, only: ns, species_names, H0i
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
    integer  :: m, s, error, ip
    real(R8) :: hprod, hcredited, qr, hph, off
    real(R8) :: qhot
    real(R8), allocatable :: dh(:), dhp(:), nu_s(:,:)
    type(file_ini) :: fini
    logical  :: apply

    allocate(dh(1:max(nmat,1)), dhp(1:max(nmat,1)))
    dh = 0._R8; dhp = 0._R8

    if (Active()) then
      if (which == MAP_IGLOO) then
        allocate(nu_s, source=nu)
      else
        allocate(nu_s, source=nuICE)
      endif

      ! datum-correction switch: false computes the offsets but does not apply them
      apply = .true.
      call fini%load(filename='input.ini')
      call fini%get(section_name='MI2-Coupling', option_name='datum-correction', val=apply, error=error)
      if (error /= 0) apply = .true.

      write(*,*)
      write(*,'(A)') ' [MI2] enthalpy-datum reconciliation  [J per kg of material converted]'
      write(*,'(A,F8.2,A)') '       anchored at T_ref = ', T_datum_ref, ' K (FLINT standard state)'
      do m = 1, nmat
        if ( book(which)%process(m) == PROC_NONE ) cycle

        hprod = 0._R8; hph = 0._R8
        do s = 1, ns
          if (nu_s(m,s) == 0._R8) cycle
          hprod = hprod + nu_s(m,s) * H0i(s, T_datum_ref)
          hph   = hph   + nu_s(m,s) * H0i(s, Thot(m))
        enddo

        write(*,'(A,A)')          '   - ', trim(material_names(m))
        do s = 1, ns
          if (nu_s(m,s) == 0._R8) cycle
          write(*,'(A,A12,A,ES15.7)') '        h_FLINT  ', trim(adjustl(species_names(s))), ' : ', &
                                      H0i(s, T_datum_ref)
        enddo
        write(*,'(A,ES15.7)')     '        products (SUM nu*h_FLINT)          : ', hprod

        select case (book(which)%energy(m))
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
          dh(m) = book(which)%h298(m) - href(m)
          qr    = book(which)%h298(m) - hprod
          write(*,'(A,ES15.7)')   '        h298 of the material (ATLAS)       : ', book(which)%h298(m)
          write(*,'(A,ES15.7)')   '        h_'//SOLVER(which)//'(table)                     : ', href(m)
          write(*,'(A,ES15.7)')   '        heat of reaction, in the species   : ', qr
          ! like with like: from the condensed material where it burns to the record's species there, beside q-comb
          qhot  = book(which)%h298(m) + (hhot(m) - href(m)) - hph
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
          if (which == MAP_IGLOO) then
            if (capped(m)) then
              ip     = book(which)%prod(m)
              dhp(m) = book(which)%hprod(m) - href(ip)
              qr     = book(which)%h298(m) - nuOxide(m)*book(which)%hprod(m)
              do s = 1, ns
                qr = qr - nuOxide(m)*nuCap(m,s) * H0i(s, T_datum_ref)
              enddo
              write(*,'(A,ES15.7)') '        h298 of its oxide (ATLAS)          : ', book(which)%hprod(m)
              write(*,'(A,ES15.7)') '        h_'//SOLVER(which)//'(table) of its oxide        : ', href(ip)
              write(*,'(A,ES15.7)') '        heat of the cap reaction, per kg   : ', qr
              if (.not. apply) dhp(m) = 0._R8
              write(*,'(A,ES15.7)') '        => datum offset of the oxide       : ', dhp(m)
            endif
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

    if (which == MAP_IGLOO) then
      if (allocated(dhDatum)) deallocate(dhDatum)
      allocate(dhDatum, source=dh)
      if (allocated(dhProduct)) deallocate(dhProduct)
      allocate(dhProduct, source=dhp)
    else
      if (allocated(dhDatumICE)) deallocate(dhDatumICE)
      allocate(dhDatumICE, source=dh)
    endif

  contains

    logical function Active()
      if (which == MAP_IGLOO) then
        Active = has_stoichiometry
      else
        Active = has_stoichiometry_ICE
      endif
    end function Active

  end subroutine Build_Datum_Offsets

end module hydraMI2_Stoichiometry
