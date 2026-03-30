program hydra_MQ2D
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use hydra_Types_m
  use hydraMQ2D_Lib_Metrics
  use hydraMQ2D_Lib_Extrapolation
  use hydraMQ2D_Mod_Explicit
  use hydraMQ2D_Mod_Coupling, only: Build_Coupling_Lookup_Table, Free_Coupling_Lookup_Table
  use Q2D_Procedures_m
  use Q2D_Types_m
  use MOSKA_Procedures_m
  use MOSKA_Types_m
  ! FLINT module
	use FLINT_Lib_Chemistry_data, only : free_chemistry_data
  ! DA CANCELLARE E SOSTITUIRE CON NUOVO FREE_CHEMISTRY_DATA
  use FLINT_Lib_Thermodynamic, only : wm_tab, Ri_tab, cp_tab, s_tab, dcpi_tab, h_tab, mi_tab, k_tab
  
  implicit none
  type(HYDRA_type)            :: HYDRA
  type(MOSKA_type)            :: MOSKA
  type(MOSKA_simulation_type) :: MOSKA_simulation
  type(Q2D_type)              :: Q2D
  type(Q2D_simulation_type)   :: Q2D_simulation
  ! Local
  integer                     :: nthreads, TODO

! Per il parallelo, valuta la possibilità di usare un numero di thread diverso per i due solver
! e.g., se uno è più pesante dell'altro, oppure soluzione pulita con task
#if defined (_OPENMP)
  !$omp parallel
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  Q2D_simulation%config%nthreads = nthreads
  MOSKA_simulation%config%nthreads = nthreads
  write(*,'(A)')   ' OpenMP:'
  write(*,'(A,I4)')' -  Number of threads --> ', nthreads
#else
  Q2D_simulation%config%nthreads = 1
  MOSKA_simulation%config%nthreads = 1
#endif

  ! Setup HYDRA
  call HYDRA%ini_setup()

  ! Setup Solvers
  !call Assign_Phases() ! Da vedere come poterlo utilizzare qui
  call Q2D%setup( Q2D_simulation )
  deallocate(wm_tab)
  deallocate(Ri_tab,cp_tab,s_tab,dcpi_tab,h_tab,mi_tab,k_tab)  
  call free_chemistry_data() ! Deallocate chemistry data loaded by Q2D
  call MOSKA%setup( MOSKA_simulation )

  ! Build lookup table for O(1) coupling boundary access
  call Build_Coupling_Lookup_Table( Q2D_simulation%domain(1), MOSKA_simulation%domain )

  call Compute_Metrics_Coupled( Q2D_simulation%domain(1), MOSKA_simulation%domain )
  call Compute_Extrapolation_Coupled( Q2D_simulation%domain(1), MOSKA_simulation%domain )

  ! Solve & Postprocess
  TODO = 1
  Q2D_simulation%config%TODO = TODO
  MOSKA_simulation%config%TODO = TODO

  do while (Q2D_simulation%config%TODO <= 2 .and.  MOSKA_simulation%config%TODO <= 2)

    call Explicit_Step ( Q2D_simulation%domain(1), MOSKA_simulation%domain, &
                         Q2D_simulation%config, MOSKA_simulation%config )
    if ( Q2D_simulation%config%TODO <= 2 ) call Q2D%postprocess( Q2D_simulation )
    if ( MOSKA_simulation%config%TODO <= 2 ) call MOSKA%postprocess( MOSKA_simulation)
  enddo

  call Q2D%postprocess( Q2D_simulation )
  call MOSKA%postprocess( MOSKA_simulation )

  ! Free coupling lookup table
  call Free_Coupling_Lookup_Table()

end program hydra_MQ2D
