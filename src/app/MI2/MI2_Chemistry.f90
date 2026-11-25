!> What the mass a condensed material loses becomes. ATLAS GPB writes, for each condensed phase,
!  INPUT/<prefix>products.txt: per material the records of the processes it can run, each a reaction:
!    evaporation  the material -> its vapour                                    (energy word "latent")
!    combustion   the material [+ gas species] -> gas species [+ a generated mass, born in an ICE family]
!    carry        the material [+ gas species] -> gas species + the carried material, which stays on the drop
!    release      a carried material handed to the gas at burnout -> gas species [+ a generated mass]
!  This module reads the records a solver's materials use into one table per solver (chemIGLOO, chemICE) and holds
!  each material's enthalpy-datum offset. Built at setup, read-only after.
module hydraMI2_Chemistry
  use iso_fortran_env, only: R8 => real64

  implicit none
  private

  public :: record_t, material_chem_t, solver_chem_t
  public :: chemIGLOO, chemICE, IS_IGLOO, IS_ICE, T_datum_ref
  public :: Load_Chemistry, Set_Datum_Offsets, Set_Product_Datum, Report_Chemistry, Gas_Of, Formation_Enthalpy

  integer,          parameter :: IS_IGLOO = 1, IS_ICE = 2
  character(len=5), parameter :: SOLVER(2) = ['IGLOO', 'ICE  ']

  !> Reference temperature of the datum offsets (FLINT standard state).
  real(R8), parameter :: T_datum_ref = 298.15_R8

  !> One record: per kg of the material it consumes, nu(s) kg of gas species s (< 0 consumed) and nuProduct kg of the
  !  condensed product (the carried material of a carry record, else a generated mass born in the ICE phase `phase`).
  !  h298 and hProduct: the formation-inclusive enthalpies at T_datum_ref of the material and of the product. dhProduct:
  !  the product's datum offset in the solver that receives it (set by the coupler).
  type :: record_t
    logical  :: on = .false., latent = .false.
    real(R8), allocatable :: nu(:)
    character(len=64) :: product = '', phase = ''
    real(R8) :: nuProduct = 0._R8, hProduct = 0._R8, h298 = 0._R8, dhProduct = 0._R8
  end type record_t

  !> A material: the record of the process it runs, of what it carries, of its release as a carried material; dh its
  !  datum offset (its formation-inclusive enthalpy less the one its solver credits, at T_datum_ref); h298 that
  !  enthalpy, when a record of its own gives it (known).
  type :: material_chem_t
    character(len=64) :: name = ''
    type(record_t)    :: process, carry, release
    real(R8)          :: dh = 0._R8, h298 = 0._R8
    logical           :: known = .false.
  end type material_chem_t

  !> A solver's table, materials in the order of its phase file. active: some material hands the gas mass.
  type :: solver_chem_t
    logical :: active = .false.
    character(len=256) :: file = ''
    type(material_chem_t), allocatable :: mat(:)
  end type solver_chem_t

  type(solver_chem_t), protected, target :: chemIGLOO, chemICE

  character(len=*), parameter :: E_ = ' [ERROR] MI2: '

contains

  function Table ( which ) result(t)
    integer, intent(in)          :: which
    type(solver_chem_t), pointer :: t
    if (which == IS_IGLOO) then
      t => chemIGLOO
    else
      t => chemICE
    endif
  end function Table


  !> Read INPUT/<prefix>products.txt for solver which. Each material keeps the records it uses: its process (combustion
  !  if it burns, else evaporation if it evaporates), its carry record when it carries a material, its release record
  !  when a material carrying it hands it to the gas at burnout (released). Stops on a record that is missing or does
  !  not close the mass, a file of another gas phase, an unknown material or species.
  subroutine Load_Chemistry ( which, prefix, names, evaporates, burns, carries, released )
    use FLINT_Lib_Thermodynamic, only: ns
    use MOSE_Global_m,           only: MOSE_phase_prefix
    integer,          intent(in) :: which
    character(len=*), intent(in) :: prefix, names(:)
    logical,          intent(in) :: evaporates(:), burns(:), carries(:), released(:)
    type(solver_chem_t), pointer :: t
    type(record_t)      :: r
    character(len=512)  :: line
    character(len=64)   :: key, mname, proc, enword, gpref, spname
    character(len=256)  :: fname
    integer  :: u, ios, m, s, k, nterm, nsfile
    real(R8) :: cmass, cmole
    logical  :: exists, gasline

    t => Table(which)
    fname  = 'INPUT/'//trim(prefix)//'products.txt'
    t%file = fname
    if (allocated(t%mat)) deallocate(t%mat)
    allocate(t%mat(size(names)))
    do m = 1, size(names)
      t%mat(m)%name = adjustl(names(m))
    enddo

    inquire(file=fname, exist=exists)
    if (.not. exists) then
      do m = 1, size(names)
        if (evaporates(m) .or. burns(m) .or. carries(m)) then
          write(*,'(A)') E_//trim(SOLVER(which))//' material '//trim(t%mat(m)%name)//' converts and '//trim(fname)// &
            ' is missing: run ATLAS GPB, which writes it'
          error stop 1
        endif
      enddo
      return
    endif

    gasline = .false.
    open(newunit=u, file=fname, status='old', action='read')
    do
      call Next_Line(u, line, ios); if (ios /= 0) exit
      read(line,*) key
      select case (trim(key))
      case ('gas')
        read(line,*,iostat=ios) key, gpref, nsfile
        if (ios == 0 .and. trim(gpref) == '-') gpref = ''
        if (ios /= 0 .or. trim(gpref) /= trim(MOSE_phase_prefix) .or. nsfile /= ns) &
          call Bad('was written for another gas phase')
        gasline = .true.
      case ('material')
        r = record_t()
        allocate(r%nu(ns), source=0._R8)
        read(line,*,iostat=ios) key, mname, proc, nterm, enword, r%h298
        if (ios /= 0) call Bad('holds a record that does not read "material <name> <process> <terms> <energy> <h298>"')
        ! the tail: the condensed product, then the phase a generated mass is born in
        read(line,*,iostat=ios) key, mname, proc, nterm, enword, r%h298, r%product, r%nuProduct, r%hProduct, r%phase
        if (ios /= 0) then
          r%phase = ''
          read(line,*,iostat=ios) key, mname, proc, nterm, enword, r%h298, r%product, r%nuProduct, r%hProduct
          if (ios /= 0) then
            r%product = ''; r%nuProduct = 0._R8; r%hProduct = 0._R8
          endif
        endif
        r%latent = trim(enword) == 'latent'
        m = Material_Index(mname)
        do k = 1, nterm
          call Next_Line(u, line, ios)
          if (ios /= 0) call Bad('ends before the terms of the record of '//trim(mname))
          read(line,*,iostat=ios) spname, cmass, cmole
          if (ios /= 0) call Bad('holds a term that does not read "<species> <nu_mass> <nu_mole>"')
          s = Species_Index(spname)
          r%nu(s) = r%nu(s) + cmass
        enddo
        r%on = .true.
        select case (trim(proc))
        case ('none')
        case ('evaporation')
          if (evaporates(m) .and. .not. burns(m)) call Keep(t%mat(m)%process)
        case ('combustion')
          if (burns(m)) call Keep(t%mat(m)%process)
        case ('carry')
          if (carries(m)) call Keep(t%mat(m)%carry)
        case ('release')
          if (released(m)) call Keep(t%mat(m)%release)
        case default
          call Bad('names the process "'//trim(proc)//'"')
        end select
      case default
        call Bad('holds the line "'//trim(key)//'"')
      end select
    enddo
    close(u)
    if (.not. gasline) call Bad('has no "gas" line')

    do m = 1, size(names)
      if ((evaporates(m) .or. burns(m)) .and. .not. t%mat(m)%process%on) &
        call Missing(m, trim(merge('combustion ', 'evaporation', burns(m))))
      if (carries(m) .and. .not. t%mat(m)%carry%on) call Missing(m, 'carry')
      if (released(m) .and. .not. t%mat(m)%release%on) call Missing(m, 'release')
    enddo
    t%active = any(evaporates .or. burns)

  contains

    !> The record read becomes slot rec, once, and closes the mass.
    subroutine Keep ( rec )
      type(record_t), intent(inout) :: rec
      if (rec%on) call Bad('holds two '//trim(proc)//' records for '//trim(mname))
      if (abs(sum(r%nu) + r%nuProduct - 1._R8) > 1.e-8_R8) &
        call Bad('holds a '//trim(proc)//' record of '//trim(mname)//' that does not close the mass')
      rec = r
    end subroutine Keep

    integer function Material_Index ( name )
      character(len=*), intent(in) :: name
      integer :: i
      Material_Index = 0
      do i = 1, size(names)
        if (trim(t%mat(i)%name) == trim(name)) Material_Index = i
      enddo
      if (Material_Index == 0) call Bad('names the material '//trim(name)//', which '//trim(SOLVER(which))// &
                                        ' does not have')
    end function Material_Index

    integer function Species_Index ( name )
      use FLINT_Lib_Thermodynamic, only: species_names
      character(len=*), intent(in) :: name
      integer :: i
      Species_Index = 0
      do i = 1, ns
        if (trim(adjustl(species_names(i))) == trim(name)) Species_Index = i
      enddo
      if (Species_Index == 0) call Bad('names the gas species '//trim(name)//', which the gas phase does not have')
    end function Species_Index

    subroutine Missing ( mm, what )
      integer,          intent(in) :: mm
      character(len=*), intent(in) :: what
      write(*,'(A)') E_//trim(fname)//' has no '//what//' record for the '//trim(SOLVER(which))//' material '// &
        trim(t%mat(mm)%name)//', which needs one: give its reaction (particle-mechanism) or its vapour (psat-vapour) '// &
        'in its [GPB-Phase*] section and run ATLAS GPB again'
      error stop 1
    end subroutine Missing

    subroutine Bad ( what )
      character(len=*), intent(in) :: what
      write(*,'(A)') E_//trim(fname)//' '//what//': run ATLAS GPB again'
      error stop 1
    end subroutine Bad

  end subroutine Load_Chemistry


  !> The next line that is neither blank nor a comment.
  subroutine Next_Line ( unit, buf, stat )
    integer,          intent(in)  :: unit
    character(len=*), intent(out) :: buf
    integer,          intent(out) :: stat
    do
      read(unit,'(A)',iostat=stat) buf; if (stat /= 0) return
      buf = adjustl(buf)
      if (len_trim(buf) > 0 .and. buf(1:1) /= '#') return
    enddo
  end subroutine Next_Line


  !> The gas species that F kg of material through record rec makes (< 0 consumed), added to r(1:ns).
  pure subroutine Gas_Of ( rec, F, r )
    type(record_t), intent(in)    :: rec
    real(R8),       intent(in)    :: F
    real(R8),       intent(inout) :: r(:)
    if (F == 0._R8 .or. .not. rec%on) return
    r(1:size(rec%nu)) = r(1:size(rec%nu)) + rec%nu*F
  end subroutine Gas_Of


  !> Each material's datum offset in solver which, from href(m), the enthalpy its solver credits it at T_datum_ref on
  !  its own datum, and Lv(m) its latent heat there: an evaporation record, the vapour's enthalpy less (href + Lv); a
  !  reaction, h298 - href. hcarried(m): the solver's enthalpy at T_datum_ref of the material m carries, for its
  !  offset. [MI2-Coupling] datum-correction = false sets every offset to 0 (diagnostic use only).
  subroutine Set_Datum_Offsets ( which, href, Lv, hcarried )
    use FLINT_Lib_Thermodynamic, only: H0i
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    integer,  intent(in) :: which
    real(R8), intent(in) :: href(:), Lv(:), hcarried(:)
    type(solver_chem_t), pointer :: t
    integer  :: m, s

    t => Table(which)
    do m = 1, size(t%mat)
      associate (c => t%mat(m))
      c%dh = 0._R8
      if (c%process%on) then
        if (c%process%latent) then
          c%dh = - href(m) - Lv(m)
          do s = 1, size(c%process%nu)
            c%dh = c%dh + c%process%nu(s)*H0i(s, T_datum_ref)
          enddo
        else
          c%dh = c%process%h298 - href(m)
        endif
      elseif (c%carry%on) then
        c%dh = c%carry%h298 - href(m)
      elseif (c%release%on) then
        c%dh = c%release%h298 - href(m)
      endif
      c%known = c%process%on .or. c%carry%on .or. c%release%on
      if (c%known) c%h298 = c%dh + href(m)
      if (c%carry%on) c%carry%dhProduct = c%carry%hProduct - hcarried(m)
      if (.not. hydra_MI2_obj%datum_correction) then
        c%dh = 0._R8
        c%carry%dhProduct = 0._R8
      endif
      end associate
    enddo

  end subroutine Set_Datum_Offsets


  !> The datum offset dh of the generated mass of material m's process record (release: of its release record), in
  !  the ICE family that receives it. datum-correction = false keeps 0.
  subroutine Set_Product_Datum ( which, m, release, dh )
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    integer,  intent(in) :: which, m
    logical,  intent(in) :: release
    real(R8), intent(in) :: dh
    type(solver_chem_t), pointer :: t

    t => Table(which)
    if (.not. hydra_MI2_obj%datum_correction) return
    if (release) then
      t%mat(m)%release%dhProduct = dh
    else
      t%mat(m)%process%dhProduct = dh
    endif

  end subroutine Set_Product_Datum


  !> The formation-inclusive enthalpy at T_datum_ref of the material called name [J/kg], from either solver's table: a
  !  record of its own, else a record whose product it is. False when no record gives it.
  logical function Formation_Enthalpy ( name, h ) result(known)
    character(len=*), intent(in)  :: name
    real(R8),         intent(out) :: h
    type(solver_chem_t), pointer :: t
    integer :: w, m

    known = .false.
    h = 0._R8
    do w = IS_IGLOO, IS_ICE
      t => Table(w)
      if (.not. allocated(t%mat)) cycle
      do m = 1, size(t%mat)
        if (trim(t%mat(m)%name) == trim(name) .and. t%mat(m)%known) then
          h = t%mat(m)%h298; known = .true.; return
        endif
      enddo
    enddo
    do w = IS_IGLOO, IS_ICE
      t => Table(w)
      if (.not. allocated(t%mat)) cycle
      do m = 1, size(t%mat)
        if (Product_Is(t%mat(m)%process) .or. Product_Is(t%mat(m)%carry) .or. Product_Is(t%mat(m)%release)) return
      enddo
    enddo

  contains

    logical function Product_Is ( rec )
      type(record_t), intent(in) :: rec
      Product_Is = rec%on .and. trim(rec%product) == trim(name)
      if (Product_Is) then
        h = rec%hProduct; known = .true.
      endif
    end function Product_Is

  end function Formation_Enthalpy


  !> Echo, per material that converts, the gas species of its records and its datum offset.
  subroutine Report_Chemistry ( which )
    use FLINT_Lib_Thermodynamic, only: species_names
    use hydraMI2_Config_Types_m, only: hydra_MI2_obj
    integer, intent(in) :: which
    type(solver_chem_t), pointer :: t
    integer :: m

    t => Table(which)
    if (.not. allocated(t%mat)) return
    if (.not. any(t%mat(:)%process%on .or. t%mat(:)%carry%on)) return
    write(*,'(A)') ' [MI2] '//trim(SOLVER(which))//' material -> gas species, per kg of material ('//trim(t%file)//')'
    do m = 1, size(t%mat)
      call Echo(t%mat(m)%process, 'process')
      call Echo(t%mat(m)%carry,   'carry  ')
      call Echo(t%mat(m)%release, 'release')
      if (t%mat(m)%process%on .or. t%mat(m)%carry%on) &
        write(*,'(A,ES15.7,A)') '       datum offset ', t%mat(m)%dh, ' J/kg'
    enddo
    if (.not. hydra_MI2_obj%datum_correction) write(*,'(A)') ' [WARNING] MI2: datum-correction = false: every offset '// &
      'is 0 and the gas receives the energy on the condensed solvers'' datums, which is wrong. Diagnostic use only'

  contains

    subroutine Echo ( rec, what )
      type(record_t),   intent(in) :: rec
      character(len=*), intent(in) :: what
      integer :: s
      character(len=1024) :: buf
      if (.not. rec%on) return
      buf = ''
      do s = 1, size(rec%nu)
        if (rec%nu(s) == 0._R8) cycle
        write(buf,'(A,1X,A,ES12.4)') trim(buf), trim(adjustl(species_names(s))), rec%nu(s)
      enddo
      if (len_trim(rec%product) > 0) &
        write(buf,'(A,A,ES12.4,A)') trim(buf)//' | ', trim(rec%product), rec%nuProduct, &
          trim(merge(' (born in '//trim(rec%phase)//')', ' (carried)               ', len_trim(rec%phase) > 0))
      write(*,'(A)') '   - '//trim(t%mat(m)%name)//' '//what//':'//trim(buf)
    end subroutine Echo

  end subroutine Report_Chemistry

end module hydraMI2_Chemistry
