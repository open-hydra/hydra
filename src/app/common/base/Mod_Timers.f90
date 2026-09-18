!>@brief Wall-clock instrumentation of the hydra coupling loop.
!>
!> The solver libraries time only their own `Explicit_Step`: the fluid solver's
!> own timers (MOSE_Mod_Timers, ARES_Mod_Timers) report the fluid and
!> FUSS_Mod_Timers the solid, and neither can see the
!> interface exchange that runs between them. This module closes that gap by
!> timing one coupling step end to end and splitting it into the three parts
!> that scale differently:
!>
!>   * FUSS     -- the solid sub-iterations
!>   * fluid    -- the fluid sub-iterations   (MOSE in hydra-MF, ARES in hydra-AF)
!>   * coupling -- `BC_Coupled_Fluxes`, of which the MPI exchange is reported
!>                 separately because that is the term that grows with rank
!>                 count as the fluid and solid partitions drift apart
!>
!> Times are reported as the **maximum** over the ranks, i.e. the critical path,
!> matching the convention of the two solver timers so the three logs can be
!> read together. Reductions go through FUSS_Mod_MPI, the one parallel module
!> linked by both hydra-MF and hydra-AF.
module hydra_Mod_Timers
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I8 => int64
  use iso_c_binding,   only: c_int, c_int64_t

  implicit none
  private

  !> Core-cycles, so that scaling can be reported with the clock divided out.
  !> A single active core turbos well above the all-core clock, so an efficiency
  !> built from seconds mixes the chip's power management into the solver's
  !> parallel behaviour. Perfect parallelism keeps TOTAL core-cycles per step
  !> constant however many cores are used, so cycles are the clock-free measure.
  interface
    subroutine hydra_cyc_open_thread(slot) bind(C, name='hydra_cyc_open_thread')
      import :: c_int
      integer(c_int), intent(in) :: slot
    end subroutine hydra_cyc_open_thread
    subroutine hydra_cyc_count(n) bind(C, name='hydra_cyc_count')
      import :: c_int
      integer(c_int), intent(out) :: n
    end subroutine hydra_cyc_count
    subroutine hydra_cyc_read_all(total) bind(C, name='hydra_cyc_read_all')
      import :: c_int64_t
      integer(c_int64_t), intent(out) :: total
    end subroutine hydra_cyc_read_all
  end interface

  logical     :: cyc_ok   = .false.   !< counter available on this build/kernel
  integer(I8) :: c_win_beg = 0_I8     !< cycles at the start of the window

  !> Reduced timings for one window of coupling steps.
  type :: stats_type
    real(R8) :: tmax  = 0.0_R8   !< coupling step on the slowest rank
    real(R8) :: fuss  = 0.0_R8   !< solid solve, max over ranks
    real(R8) :: fluid = 0.0_R8   !< fluid solve, max over ranks
    real(R8) :: coup  = 0.0_R8   !< interface exchange incl. local work, max
    real(R8) :: xchg  = 0.0_R8   !< the MPI part of `coup`, max
  end type stats_type

  real(R8) :: t_step_beg = 0.0_R8, t_fuss_beg = 0.0_R8, t_fluid_beg = 0.0_R8
  real(R8) :: t_coup_beg = 0.0_R8, t_xchg_beg = 0.0_R8

  real(R8) :: a_step = 0.0_R8, a_fuss = 0.0_R8, a_fluid = 0.0_R8
  real(R8) :: a_coup = 0.0_R8, a_xchg = 0.0_R8
  integer  :: n_acc = 0

  real(R8) :: t_run_beg = 0.0_R8
  real(R8) :: r_step = 0.0_R8, r_fuss = 0.0_R8, r_fluid = 0.0_R8
  real(R8) :: r_coup = 0.0_R8, r_xchg = 0.0_R8
  integer  :: n_run = 0

  !> Interface entries owned by this rank whose partner is owned by another
  !> rank, and the total. Filled once by the coupling schedule builder.
  integer, public :: n_iface_mpi   = 0
  integer, public :: n_iface_total = 0

  !> Name of the fluid solver, used only in the reports. Each executable sets
  !> it at start-up: 'MOSE' in hydra-MF, 'ARES' in hydra-AF.
  character(len=4), public :: timer_fluid_name = 'MOSE'

  public :: timer_wtime
  public :: timer_run_begin, timer_summary, timer_report
  public :: timer_step_begin, timer_step_end
  public :: timer_fuss_begin, timer_fuss_end
  public :: timer_fluid_begin, timer_fluid_end
  public :: timer_coup_begin, timer_coup_end
  public :: timer_xchg_begin, timer_xchg_end
  public :: timer_set_interface_counts

contains

  !> Wall-clock reading, in seconds. MPI_WTIME when available (monotonic and
  !> consistent across the ranks of a job), SYSTEM_CLOCK otherwise.
  real(R8) function timer_wtime()
#ifdef USE_MPI
    timer_wtime = MPI_WTIME()
#else
    integer(I8) :: count, rate
    call system_clock(count, rate)
    timer_wtime = real(count, R8) / real(rate, R8)
#endif
  end function timer_wtime


  !> Mark the start of the coupling loop, i.e. the end of set-up.
  subroutine timer_run_begin()
    ! The `!$` sentinel, not #ifdef USE_OPENMP: USE_OPENMP is a CMake option and
    ! is never passed to the compiler as a definition, so an #ifdef on it is
    ! always false. That silently gave every thread slot 0 and counted one.
    !$ use omp_lib, only: omp_get_thread_num
    integer(c_int)     :: slot, nc
    integer(c_int64_t) :: c0

    t_run_beg = timer_wtime()

    ! Attach a counter to every thread of the pool. Done here, after set-up, so
    ! the pool already has its full complement; the counters then persist for
    ! the run and each window is read as a difference.
    !$omp parallel private(slot)
    slot = 0_c_int
    !$ slot = int(omp_get_thread_num(), c_int)
    call hydra_cyc_open_thread(slot)
    !$omp end parallel

    call hydra_cyc_count(nc)
    cyc_ok = (nc > 0)
    if (cyc_ok) then
      call hydra_cyc_read_all(c0)
      c_win_beg = int(c0, I8)
    end if
    r_step = 0.0_R8; r_fuss = 0.0_R8; r_fluid = 0.0_R8
    r_coup = 0.0_R8; r_xchg = 0.0_R8
    n_run  = 0

    ! Window accumulators too: Setup_Coupled_Metrics runs the metric-node
    ! interface exchange before the coupling loop starts, and that exchange is
    ! timed. Without this it is charged to the first reported window.
    a_step = 0.0_R8; a_fuss = 0.0_R8; a_fluid = 0.0_R8
    a_coup = 0.0_R8; a_xchg = 0.0_R8
    n_acc  = 0
  end subroutine timer_run_begin

  subroutine timer_step_begin()
    t_step_beg = timer_wtime()
  end subroutine timer_step_begin

  subroutine timer_step_end()
    real(R8) :: dt
    dt = timer_wtime() - t_step_beg
    a_step = a_step + dt;  r_step = r_step + dt
    n_acc  = n_acc + 1;    n_run  = n_run + 1
  end subroutine timer_step_end

  subroutine timer_fuss_begin()
    t_fuss_beg = timer_wtime()
  end subroutine timer_fuss_begin

  subroutine timer_fuss_end()
    real(R8) :: dt
    dt = timer_wtime() - t_fuss_beg
    a_fuss = a_fuss + dt; r_fuss = r_fuss + dt
  end subroutine timer_fuss_end

  subroutine timer_fluid_begin()
    t_fluid_beg = timer_wtime()
  end subroutine timer_fluid_begin

  subroutine timer_fluid_end()
    real(R8) :: dt
    dt = timer_wtime() - t_fluid_beg
    a_fluid = a_fluid + dt; r_fluid = r_fluid + dt
  end subroutine timer_fluid_end

  !> Whole of `BC_Coupled_Fluxes`: pack, unpack, the local-local ghost fill and
  !> the flux assembly, as well as the MPI below.
  subroutine timer_coup_begin()
    t_coup_beg = timer_wtime()
  end subroutine timer_coup_begin

  subroutine timer_coup_end()
    real(R8) :: dt
    dt = timer_wtime() - t_coup_beg
    a_coup = a_coup + dt; r_coup = r_coup + dt
  end subroutine timer_coup_end

  !> The MPI part alone: the IRECV/ISEND/WAITALL of the four exchange phases.
  subroutine timer_xchg_begin()
    t_xchg_beg = timer_wtime()
  end subroutine timer_xchg_begin

  subroutine timer_xchg_end()
    real(R8) :: dt
    dt = timer_wtime() - t_xchg_beg
    a_xchg = a_xchg + dt; r_xchg = r_xchg + dt
  end subroutine timer_xchg_end


  !> Record how much of the interface this rank has to reach over MPI.
  !> The fluid solver and FUSS partition their own phase independently, so a fluid block
  !> and the solid block under it need not share a rank; this is the fraction
  !> that does not, which is what turns interface work into messages.
  subroutine timer_set_interface_counts(n_mpi, n_total)
    integer, intent(in) :: n_mpi, n_total
    n_iface_mpi   = n_mpi
    n_iface_total = n_total
  end subroutine timer_set_interface_counts


  !> Reduce the accumulators over the ranks. Collective.
  function reduce_times(t_step, t_solid, t_fluid, t_coup, t_xchg, n) result(s)
    use FUSS_Mod_MPI, only: mpi_allreduce_max_r8
    real(R8), intent(in) :: t_step, t_solid, t_fluid, t_coup, t_xchg
    integer,  intent(in) :: n
    type(stats_type)     :: s
    real(R8) :: rn

    rn = real(max(n, 1), R8)
    call mpi_allreduce_max_r8(t_step,  s%tmax)
    call mpi_allreduce_max_r8(t_solid, s%fuss)
    call mpi_allreduce_max_r8(t_fluid, s%fluid)
    call mpi_allreduce_max_r8(t_coup,  s%coup)
    call mpi_allreduce_max_r8(t_xchg,  s%xchg)

    s%tmax = s%tmax / rn; s%fuss = s%fuss / rn; s%fluid = s%fluid / rn
    s%coup = s%coup / rn; s%xchg = s%xchg / rn
  end function reduce_times


  !> Sum a per-rank cycle count over the communicator, as a real so that the
  !> 64-bit totals at high core counts stay exact enough to print.
  real(R8) function sum_cycles_over_ranks(c) result(total)
    integer(I8), intent(in) :: c
    real(R8) :: local
#ifdef USE_MPI
    real(R8) :: g
    integer  :: ierr
#endif
    local = real(c, R8)
#ifdef USE_MPI
    call MPI_ALLREDUCE(local, g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    total = g
#else
    total = local
#endif
  end function sum_cycles_over_ranks


  !> Percentage of `whole`, guarded against a zero denominator.
  pure real(R8) function pct(part, whole)
    real(R8), intent(in) :: part, whole
    if (whole > 0.0_R8) then
      pct = part / whole * 100.0_R8
    else
      pct = 0.0_R8
    end if
  end function pct


  !> Report the current window and reset it. Collective.
  subroutine timer_report(step)
    use FUSS_Mod_MPI, only: mpi_is_root
    integer, intent(in) :: step
    type(stats_type) :: s
    integer(I8) :: c_now, c_win
    real(R8)    :: cyc_step, ghz

    s = reduce_times(a_step, a_fuss, a_fluid, a_coup, a_xchg, n_acc)

    ! Total core-cycles spent in this window, summed over every rank and the
    ! threads it owns. Divided by the steps in the window this is the clock-free
    ! cost of a coupling step: compare it across core counts and the ratio is
    ! parallel efficiency with the processor's frequency scaling removed.
    cyc_step = 0.0_R8
    ghz      = 0.0_R8
    if (cyc_ok) then
      block
        integer(c_int64_t) :: craw
        call hydra_cyc_read_all(craw)
        c_now = int(craw, I8)
      end block
      c_win = c_now - c_win_beg
      c_win_beg = c_now
      cyc_step = sum_cycles_over_ranks(c_win) / real(max(n_acc,1), R8)
      if (s%tmax > 0.0_R8) ghz = cyc_step / s%tmax / 1.0e9_R8
    end if

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A,ES11.4,A,ES11.4,A,F5.1,A,ES11.4,A,F5.1,A,ES11.4,A,F5.1,A,ES11.4,A,F5.1,A)') &
        ' HYDRA Timing | Step ', step, ' | ', max(n_acc,1), &
        ' steps | wall/step ', s%tmax, &
        ' s | FUSS ', s%fuss, ' s (', pct(s%fuss, s%tmax), &
        ' %) | '//timer_fluid_name//' ', s%fluid, ' s (', pct(s%fluid, s%tmax), &
        ' %) | coupling ', s%coup, ' s (', pct(s%coup, s%tmax), &
        ' %) | of which MPI ', s%xchg, ' s (', pct(s%xchg, s%coup), ' %)'
      if (cyc_ok) &
        write(*,'(A,ES11.4,A,F9.2,A)') &
          ' HYDRA Cycles | core-cycles/step ', cyc_step, &
          ' | aggregate ', ghz, ' GHz'
    end if

    a_step = 0.0_R8; a_fuss = 0.0_R8; a_fluid = 0.0_R8
    a_coup = 0.0_R8; a_xchg = 0.0_R8
    n_acc  = 0
  end subroutine timer_report


  !> End-of-run summary. Collective, so every rank must call it.
  subroutine timer_summary()
    use FUSS_Mod_MPI, only: mpi_is_root, mpi_size_, mpi_allreduce_max_r8
    type(stats_type) :: s
    real(R8) :: elapsed
    integer  :: iface_mpi, iface_tot

    call mpi_allreduce_max_r8(timer_wtime() - t_run_beg, elapsed)
    s = reduce_times(r_step, r_fuss, r_fluid, r_coup, r_xchg, n_run)
    call sum_i(n_iface_mpi,   iface_mpi)
    call sum_i(n_iface_total, iface_tot)

    if (mpi_is_root) then
      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Timing (HYDRA coupling)'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A,T35,I0)')       '   Coupling steps', n_run
      write(*,'(A,T35,ES12.5,A)') '   Coupling loop', s%tmax * real(max(n_run,1), R8), ' s'
      write(*,'(A,T35,ES12.5,A)') '   Coupling step', s%tmax, ' s'
      write(*,'(A,T35,ES12.5,A)') '   Elapsed', elapsed, ' s'
      write(*,'(A,T35,ES12.5,A,F7.1,A)') '   FUSS per step', s%fuss, ' s', pct(s%fuss, s%tmax), ' %'
      write(*,'(A,T35,ES12.5,A,F7.1,A)') '   '//timer_fluid_name//' per step', s%fluid, ' s', pct(s%fluid, s%tmax), ' %'
      write(*,'(A,T35,ES12.5,A,F7.1,A)') '   Coupling per step', s%coup, ' s', pct(s%coup, s%tmax), ' %'
      if (mpi_size_ > 1) then
        write(*,'(A,T35,ES12.5,A,F7.1,A)') '   Coupling MPI wait', s%xchg, ' s', pct(s%xchg, s%coup), ' %'
        if (iface_tot > 0) &
          write(*,'(A,T35,I0,A,I0,A,F7.1,A)') '   Interface over MPI', iface_mpi, ' of ', &
            iface_tot, ' entries', real(iface_mpi,R8)/real(iface_tot,R8)*100.0_R8, ' %'
      end if
      write(*,'(A)') ' ========================================================================================='
    end if

  end subroutine timer_summary


  !> MPI_ALLREDUCE with MPI_SUM for one integer; identity in serial.
  subroutine sum_i(local_val, global_val)
    integer, intent(in)  :: local_val
    integer, intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    global_val = local_val
#endif
  end subroutine sum_i

end module hydra_Mod_Timers
