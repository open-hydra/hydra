module hydraMI2_Stoichiometry
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Build_Species_Stoichiometry
  public :: nu, has_stoichiometry, Report_Stoichiometry
  public :: Build_Datum_Offsets, dhDatum, T_datum_ref

  !> nu(m,s) = kg of gas species s exchanged per kg of IGLOO material m converted
  !  (nu > 0 produced, nu < 0 consumed). Row sum is 1 for pure gasification.
  real(R8), allocatable, protected :: nu(:,:)
  logical,               protected :: has_stoichiometry = .false.

  !> Reference temperature of the datum offset (FLINT standard state).
  real(R8), parameter :: T_datum_ref = 298.15_R8

  !> dhDatum(m) = constant enthalpy-datum offset for material m [J per kg converted]: FLINT's
  !  formation-inclusive enthalpy of the products minus what IGLOO already credits (h_liq + Lv).
  real(R8), allocatable, protected :: dhDatum(:)

contains

  !> Read [MI2-Coupling] products-<material> = <species>:<coefficient>, ... from input.ini and
  !  resolve it by species name into nu(nmat,nsc).
  subroutine Build_Species_Stoichiometry ( material_names, nmat )
    use Finer,                    only: file_ini
    use FLINT_Lib_Thermodynamic, only: ns, species_names
    implicit none
    character(len=*), intent(in) :: material_names(:)
    integer,          intent(in) :: nmat
    ! Local
    type(file_ini)       :: fini
    character(len=1024)  :: raw
    character(len=256)   :: token, spname
    integer              :: m, s, error, p0, p1, icolon, nterm
    real(R8)             :: coeff
    logical              :: found

    if (allocated(nu)) deallocate(nu)
    allocate(nu(1:max(nmat,1), 1:max(ns,1)))
    nu = 0._R8
    has_stoichiometry = .false.

    call fini%load(filename='input.ini')

    do m = 1, nmat
      raw = ''
      call fini%get(section_name='MI2-Coupling',                                  &
                    option_name='products-'//trim(adjustl(material_names(m))),    &
                    val=raw, error=error)
      if (error /= 0 .or. len_trim(raw) == 0) cycle

      has_stoichiometry = .true.
      nterm = 0
      p0 = 1
      do
        if (p0 > len_trim(raw)) exit
        p1 = index(raw(p0:), ',')
        if (p1 == 0) then
          token = adjustl(raw(p0:len_trim(raw)))
          p0 = len_trim(raw) + 1
        else
          token = adjustl(raw(p0:p0+p1-2))
          p0 = p0 + p1
        endif
        if (len_trim(token) == 0) cycle

        icolon = index(token, ':')
        if (icolon <= 1) then
          write(*,'(A)') '[ERROR] [MI2-Coupling] products-'//trim(material_names(m))// &
                         ': expected <species>:<coefficient>, got "'//trim(token)//'"'
          error stop 1
        endif
        spname = adjustl(token(1:icolon-1))
        read(token(icolon+1:),*,iostat=error) coeff
        if (error /= 0) then
          write(*,'(A)') '[ERROR] [MI2-Coupling] products-'//trim(material_names(m))// &
                         ': cannot parse coefficient in "'//trim(token)//'"'
          error stop 1
        endif

        ! resolve the species by name against the gas table
        found = .false.
        do s = 1, ns
          if (trim(adjustl(species_names(s))) == trim(spname)) then
            nu(m,s) = nu(m,s) + coeff
            found = .true.
            exit
          endif
        enddo
        if (.not. found) then
          write(*,'(A)') '[ERROR] [MI2-Coupling] products-'//trim(material_names(m))// &
                         ': unknown gas species "'//trim(spname)//'". Available:'
          do s = 1, ns
            write(*,'(A,A)') '          - ', trim(adjustl(species_names(s)))
          enddo
          error stop 1
        endif
        nterm = nterm + 1
      enddo

      if (nterm == 0) then
        write(*,'(A)') '[ERROR] [MI2-Coupling] products-'//trim(material_names(m))//' is empty'
        error stop 1
      endif
    enddo

  end subroutine Build_Species_Stoichiometry


  !> Echo the resolved stoichiometry matrix with the mass-balance reading of each row.
  subroutine Report_Stoichiometry ( material_names, nmat )
    use FLINT_Lib_Thermodynamic, only: ns, species_names
    implicit none
    character(len=*), intent(in) :: material_names(:)
    integer,          intent(in) :: nmat
    integer  :: m, s
    real(R8) :: rowsum

    if (.not. has_stoichiometry) then
      write(*,'(A)') ' [MI2] no [MI2-Coupling] products-<material> given: mass coupling is OFF.'
      return
    endif

    write(*,*)
    write(*,'(A)') ' [MI2] material -> gas-species stoichiometry [kg gas per kg material converted]'
    do m = 1, nmat
      rowsum = sum(nu(m,1:ns))
      if (all(nu(m,1:ns) == 0._R8)) then
        write(*,'(A,A,A)') '   - ', trim(material_names(m)), ' : (none -- inert, no mass coupling)'
        cycle
      endif
      write(*,'(A,A,A)') '   - ', trim(material_names(m)), ' :'
      do s = 1, ns
        if (nu(m,s) == 0._R8) cycle
        write(*,'(A,A12,A,F10.5,A)') '        ', trim(adjustl(species_names(s))), ' : ', nu(m,s), &
          merge('  (consumed)', '  (produced)', nu(m,s) < 0._R8)
      enddo
      write(*,'(A,F10.5)') '        SUM (=1 for pure gasification) : ', rowsum
      if (rowsum > 1._R8 + 1.e-10_R8) then
        write(*,'(A)') '        [WARNING] sum > 1: the gas gains more mass than the material loses.'
        write(*,'(A)') '                  Only valid if a consumed species is missing from the list.'
      endif
    enddo
    write(*,*)

  end subroutine Report_Stoichiometry


  !> Build dhDatum(m) = SUM_s nu(m,s)*h_FLINT(s,T_ref) - (hIGLOOref(m) + Lv(m)) and print the
  !  per-material datum gap. [MI2-Coupling] datum-correction = false zeroes the offsets.
  subroutine Build_Datum_Offsets ( material_names, nmat, hIGLOOref, Lv )
    use Finer,                   only: file_ini
    use FLINT_Lib_Thermodynamic, only: ns, species_names, H0i
    implicit none
    character(len=*), intent(in) :: material_names(:)
    integer,          intent(in) :: nmat
    real(R8),         intent(in) :: hIGLOOref(:)   !> h IGLOO uses at T_ref, its own datum [J/kg]
    real(R8),         intent(in) :: Lv(:)          !> per-material latent heat [J/kg]
    integer  :: m, s, error
    real(R8) :: hprod, hcredited
    type(file_ini) :: fini
    logical  :: apply

    if (allocated(dhDatum)) deallocate(dhDatum)
    allocate(dhDatum(1:max(nmat,1)))
    dhDatum = 0._R8

    if (.not. has_stoichiometry) return

    ! datum-correction switch: false computes the offsets but does not apply them
    apply = .true.
    call fini%load(filename='input.ini')
    call fini%get(section_name='MI2-Coupling', option_name='datum-correction', val=apply, error=error)
    if (error /= 0) apply = .true.

    write(*,*)
    write(*,'(A)') ' [MI2] enthalpy-datum reconciliation  [J per kg of material converted]'
    write(*,'(A,F8.2,A)') '       anchored at T_ref = ', T_datum_ref, ' K (FLINT standard state)'
    do m = 1, nmat
      if ( all(nu(m,1:ns) == 0._R8) ) cycle

      hprod = 0._R8
      do s = 1, ns
        if (nu(m,s) == 0._R8) cycle
        hprod = hprod + nu(m,s) * H0i(s, T_datum_ref)
      enddo

      hcredited  = hIGLOOref(m) + Lv(m)
      dhDatum(m) = hprod - hcredited
      if (.not. apply) dhDatum(m) = 0._R8

      write(*,'(A,A)')          '   - ', trim(material_names(m))
      do s = 1, ns
        if (nu(m,s) == 0._R8) cycle
        write(*,'(A,A12,A,ES15.7)') '        h_FLINT  ', trim(adjustl(species_names(s))), ' : ', &
                                    H0i(s, T_datum_ref)
      enddo
      write(*,'(A,ES15.7)')     '        products (SUM nu*h_FLINT)          : ', hprod
      write(*,'(A,ES15.7)')     '        h_IGLOO(liquid table)              : ', hIGLOOref(m)
      write(*,'(A,ES15.7)')     '        Lv (already inside sourceEn)       : ', Lv(m)
      write(*,'(A,ES15.7)')     '        credited by IGLOO (h_IGLOO + Lv)   : ', hcredited
      write(*,'(A,ES15.7)')     '        => datum offset applied by MI2     : ', dhDatum(m)
      if ( apply .and. abs(dhDatum(m)) < 1.e-6_R8*max(abs(hprod),1._R8) ) then
        write(*,'(A)')          '        (offset ~ 0: the two datums already agree)'
      endif
    enddo
    if (.not. apply) then
      write(*,'(A)') '   [WARNING] datum-correction = false: the offsets above are computed but NOT'
      write(*,'(A)') '             applied. The gas receives energy on IGLOO''s liquid-table datum,'
      write(*,'(A)') '             which is KNOWN WRONG. Diagnostic use only.'
    endif
    write(*,*)

  end subroutine Build_Datum_Offsets

end module hydraMI2_Stoichiometry
