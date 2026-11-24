!> The channels of hydra-MI2: what one solver's entity hands another's across their steps.
!
!  A channel carries amounts that a sender makes over its own clock into a receiver that steps on another one. It
!  keeps them in its store, cell by cell, with the time they were made over, and hands them on in one of two ways
!  (As_Rate):
!    as a rate    when either end steps in pseudo-time: the store over the time it was made in, kept in force until
!                 the sender makes more (a converged steady field holds per unit time)
!    as amounts   when both ends are time-accurate: what the sender made, whole, over the receiver's next steps
!  This module knows no solver: the coupler (hydraMI2_Mod_Coupling) fills a channel from its sender and empties it
!  into its receiver.
module hydraMI2_Channels
  use iso_fortran_env, only: R8 => real64

  implicit none
  private

  public :: channel_t, store_t
  public :: As_Rate, Open_Store, Clock, Store_To_Rate, Open_Share, Close_Share, Keep, Empty

  !> A channel's store in one block, cell by cell: a(component, i, j, k) the amount made and not yet handed, t the
  !  time it was made over, r the rate in force (a channel handed as a rate).
  type :: store_t
    real(R8), allocatable :: a(:,:,:,:), t(:,:,:), r(:,:,:,:)
  end type store_t

  !> One route from a sender's entity to a receiver's (the solver's family or material; 0: the gas). What it
  !  carries per kg of what triggers it (nu), its datum offset per kg carried (dh), and the diameter a packet is born
  !  with in a condensed receiver (d). asRate: how it is handed; share: under amounts, the part of the store each
  !  iteration of the receiver's current block takes (0 outside one). Its store, per block.
  type :: channel_t
    character(len=40) :: name    = ''
    integer           :: from    = 0, to = 0, trigger = 0
    real(R8)          :: nu      = 0._R8, dh = 0._R8, d = 0._R8
    logical           :: asRate  = .false.
    real(R8)          :: share   = 0._R8
    type(store_t), allocatable :: blk(:)
  end type channel_t

contains

  !> The clock rule: a channel goes as a rate unless both its ends are time-accurate.
  pure logical function As_Rate ( senderTA, receiverTA )
    logical, intent(in) :: senderTA, receiverTA

    As_Rate = .not. (senderTA .and. receiverTA)

  end function As_Rate


  !> The store of ch over blocks of interior dimensions dims(1:3, b), ncomp components, zero: the amounts, the time
  !  when timed, the rate when rated. What is already there is kept: a store opened again gains what it lacks.
  subroutine Open_Store ( ch, dims, ncomp, timed, rated )
    type(channel_t), intent(inout) :: ch
    integer,         intent(in)    :: dims(:,:), ncomp
    logical,         intent(in)    :: timed, rated
    integer :: b, n(3)

    if (.not. allocated(ch%blk)) allocate(ch%blk(size(dims, 2)))
    do b = 1, size(ch%blk)
      n = dims(1:3, b)
      if (.not. allocated(ch%blk(b)%a)) allocate(ch%blk(b)%a(ncomp, n(1), n(2), n(3)), source=0._R8)
      if (timed .and. .not. allocated(ch%blk(b)%t)) allocate(ch%blk(b)%t(n(1), n(2), n(3)), source=0._R8)
      if (rated .and. .not. allocated(ch%blk(b)%r)) allocate(ch%blk(b)%r(ncomp, n(1), n(2), n(3)), source=0._R8)
    enddo

  end subroutine Open_Store


  !> After a sender's step of dt (cell by cell) in block b: the time its store is being made over.
  subroutine Clock ( ch, b, dt )
    type(channel_t), intent(inout) :: ch
    integer,         intent(in)    :: b
    real(R8),        intent(in)    :: dt(:,:,:)
    integer :: i, j, k

    do k = 1, size(ch%blk(b)%t, 3)
    do j = 1, size(ch%blk(b)%t, 2)
    do i = 1, size(ch%blk(b)%t, 1)
      ch%blk(b)%t(i,j,k) = ch%blk(b)%t(i,j,k) + dt(i,j,k)
    enddo; enddo; enddo

  end subroutine Clock


  !> As a rate: where the sender made something since the last call (time > 0), the store over that time becomes the
  !  rate in force, and the store and its time are emptied. Elsewhere the rate in force stays.
  subroutine Store_To_Rate ( ch )
    type(channel_t), intent(inout) :: ch
    integer  :: b, i, j, k
    real(R8) :: t

    do b = 1, size(ch%blk)
      do k = 1, size(ch%blk(b)%t, 3)
      do j = 1, size(ch%blk(b)%t, 2)
      do i = 1, size(ch%blk(b)%t, 1)
        t = ch%blk(b)%t(i,j,k)
        if (.not. t > 0._R8) cycle
        ch%blk(b)%r(:,i,j,k) = ch%blk(b)%a(:,i,j,k)/t
        ch%blk(b)%a(:,i,j,k) = 0._R8
        ch%blk(b)%t(i,j,k)   = 0._R8
      enddo; enddo; enddo
    enddo

  end subroutine Store_To_Rate


  !> As amounts, at the start of a receiver's block of nIter iterations: each iteration takes 1/nIter of the store.
  subroutine Open_Share ( ch, nIter )
    type(channel_t), intent(inout) :: ch
    integer,         intent(in)    :: nIter

    if (nIter > 0) ch%share = 1._R8 / real(nIter, R8)

  end subroutine Open_Share


  !> After a receiver's block that made nRun of the iterations Open_Share shared the store over: the part of the store
  !  they took (0 when no share was open); the share closed. The caller books that part, then keeps the rest (Keep).
  real(R8) function Close_Share ( ch, nRun ) result( frac )
    type(channel_t), intent(inout) :: ch
    integer,         intent(in)    :: nRun

    frac = 0._R8
    if (.not. ch%share > 0._R8) return
    frac     = real(nRun, R8)*ch%share
    ch%share = 0._R8

  end function Close_Share


  !> The store times f: what it keeps once a part of it has been handed.
  subroutine Keep ( ch, f )
    type(channel_t), intent(inout) :: ch
    real(R8),        intent(in)    :: f
    integer :: b

    do b = 1, size(ch%blk)
      ch%blk(b)%a = ch%blk(b)%a*f
    enddo

  end subroutine Keep


  !> The store and its time emptied (what it held has been handed whole).
  subroutine Empty ( ch )
    type(channel_t), intent(inout) :: ch
    integer :: b

    do b = 1, size(ch%blk)
      ch%blk(b)%a = 0._R8
      if (allocated(ch%blk(b)%t)) ch%blk(b)%t = 0._R8
    enddo

  end subroutine Empty

end module hydraMI2_Channels
