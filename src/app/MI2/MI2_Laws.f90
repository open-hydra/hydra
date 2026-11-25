!> The laws of hydra-MI2's mass channels ([MI2-Channels]): how much mass a channel moves in a cell.
!
!  A law is one pure function of what the channel sees in a cell (law_cell_t) and of its constants (the key=value
!  pairs of its line), returning a mass rate per unit volume [kg/(m^3 s)], >= 0. The channel moves that mass from its
!  donor to its receiver; the coupler does the rest (the donor gives it up, the receiver takes it in, the gas keeps
!  the books). To add a law:
!    1. write its function at the end of this module;
!    2. add its name to LAW_NAMES and its constants to Law_Setup;
!    3. add a case to Law_Rate.
!  Nothing else of hydra changes.
module hydraMI2_Laws
  use iso_fortran_env, only: R8 => real64

  implicit none
  private

  public :: law_t, law_cell_t, Law_Setup, Law_Rate

  !> What a law sees in a cell: the gas (density, velocity, temperature, pressure, and the density of the donor's
  !  species when the donor is the gas), the donor's and the receiver's mass per unit volume [kg/m^3] (< 0: not known,
  !  an IGLOO end), the cell volume.
  type :: law_cell_t
    real(R8) :: rho = 0._R8, u(3) = 0._R8, T = 0._R8, p = 0._R8
    real(R8) :: donor = -1._R8, receiver = -1._R8
    real(R8) :: vol = 0._R8
  end type law_cell_t

  !> A law: its name, its index in LAW_NAMES, its constants.
  type :: law_t
    character(len=32) :: name = ''
    integer  :: id = 0
    real(R8) :: c(4) = 0._R8
    logical  :: needsDonor = .false.     ! its rate depends on the donor's mass (refused for an IGLOO donor)
  end type law_t

  integer,           parameter :: LAW_CONSTANT = 1, LAW_RELAXATION = 2
  character(len=32), parameter :: LAW_NAMES(2) = [character(len=32) :: 'constant', 'relaxation']

contains

  !> The law called `name` with the constants of keys(:) = values(:). why is empty when it is taken, else says why
  !  not.
  function Law_Setup ( name, keys, values, why ) result(law)
    character(len=*),              intent(in)  :: name, keys(:), values(:)
    character(len=:), allocatable, intent(out) :: why
    type(law_t) :: law

    why = ''
    law%name = name
    law%id   = findloc(LAW_NAMES, name, dim=1)
    select case (law%id)
    case (LAW_CONSTANT)
      ! constant rate=<kg/(m^3 s)>: the same rate in every cell
      law%c(1) = Constant_Of('rate')
    case (LAW_RELAXATION)
      ! relaxation k=<1/s>: the donor's mass relaxes to zero at the rate k
      law%c(1) = Constant_Of('k')
      law%needsDonor = .true.
    case default
      why = 'no law is called '//trim(name)//' (the laws: '//trim(Names())//')'
    end select

  contains

    !> The value of constant `key`, which must be given, finite and >= 0.
    real(R8) function Constant_Of ( key ) result(x)
      character(len=*), intent(in) :: key
      integer :: k, ios
      x = 0._R8
      k = findloc(keys, key, dim=1)
      if (k == 0) then
        why = 'the law '//trim(name)//' needs '//key//'='
        return
      endif
      read(values(k),*,iostat=ios) x
      if (ios /= 0 .or. .not. (x >= 0._R8 .and. x <= huge(x))) why = key//'='//trim(values(k))//' is not a real >= 0'
    end function Constant_Of

    function Names () result(s)
      character(len=:), allocatable :: s
      integer :: i
      s = ''
      do i = 1, size(LAW_NAMES)
        s = s//' '//trim(LAW_NAMES(i))
      enddo
    end function Names

  end function Law_Setup


  !> The mass rate law moves in a cell of state s [kg/(m^3 s)].
  pure real(R8) function Law_Rate ( law, s ) result(rate)
    type(law_t),      intent(in) :: law
    type(law_cell_t), intent(in) :: s

    select case (law%id)
    case (LAW_CONSTANT);   rate = law%c(1)
    case (LAW_RELAXATION); rate = law%c(1)*max(s%donor, 0._R8)
    case default;          rate = 0._R8
    end select

  end function Law_Rate

end module hydraMI2_Laws
