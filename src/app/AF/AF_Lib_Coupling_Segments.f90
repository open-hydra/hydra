!> Type 104 (chimera coupling) as a non-conforming 103, FUSS <-> ARES.
!>
!> Port of hydraMF_Lib_Coupling_Segments. FUSS and ARES meet along the interface
!> without overlapping, so the coupling does not need the overset blend over the
!> ghost cells. The interface is split into segments, the overlaps between one
!> FUSS 104 face and one ARES 104 face, and each segment is coupled as a 103 pair
!> of faces:
!>  - the 103 wall-temperature routines are called with the real partner cell,
!>    the cell of the other solver that owns the other face of the segment;
!>  - the energy flux of the segment is 0.5*own - 0.5*partner estimate on both
!>    sides, so what leaves one solver enters the other exactly;
!>  - as in 103, both sides build the interface temperature with the same fluid
!>    wall conductivity k_l+k_t, the one the ARES routine evaluated on the segment
!>    at the previous call (kS, one per segment).
!> On conforming faces there is one segment per face and the scheme is 103.
!>
!> Still to come:
!>  - partner metric: only the normal distance of the partner cell centre is
!>    used, exact when the mesh is orthogonal at the interface; the partner
!>    tangential terms are dropped. The full metric, oriented as
!>    Setup_Coupled_Metrics does for 103, is the next step.
!>  - MPI: partner cells are read in place, so 104 runs serial
!>    (Check_Chimera_Coupling_Donors).
module hydraAF_Lib_Coupling_Segments
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use ARES_Advanced_Types_m

  implicit none
  private
  public :: Check_Chimera_Coupling_Donors, Setup_Coupling_Segments, Compute_Coupling_Segment_Fluxes
  public :: Segment_Wall_Output_FUSS, Segment_Wall_Output_ARES

  integer, parameter :: max_levels = 16

  !> Segments of the 104 interface of one multigrid level
  type :: segment_set_type
    logical :: ready = .false.
    integer :: ns = 0                            ! number of segments
    integer,  allocatable :: sF(:), sA(:)        ! FUSS / ARES bc index of the two faces of the segment
    real(R8), allocatable :: wF(:), wA(:)        ! segment area over the FUSS / ARES face area
    real(R8), allocatable :: dF(:), dA(:)        ! distance of the FUSS / ARES cell centre from the interface
    integer :: nF = 0, nA = 0                    ! number of FUSS / ARES 104 faces
    integer,  allocatable :: eF(:), eA(:)        ! bc index of each 104 face
    integer,  allocatable :: rowF(:), rowA(:)    ! 104 face number of a bc index, 0 if not a 104 face
    integer,  allocatable :: pF(:), pA(:)        ! segments of face e: l(p(e) : p(e+1)-1)
    integer,  allocatable :: lF(:), lA(:)
    real(R8), allocatable :: qF(:)               ! FUSS estimate of the energy flux through the segment
    real(R8), allocatable :: fA(:,:)             ! ARES estimate of the fluxes through the segment (nprim, ns)
    real(R8), allocatable :: kS(:)               ! fluid wall conductivity k_l+k_t of the segment, lagged by one call
  end type segment_set_type

  type(segment_set_type), target, save :: segs(max_levels)

  !> Geometry of one 104 face, used at setup only
  type :: face_geo_type
    real(R8) :: p(3,4)                           ! nodes, in order around the face
    real(R8) :: c(3)                             ! centre
    real(R8) :: n(3)                             ! unit normal, out of the own cell
    real(R8) :: d                                ! distance of the own cell centre, as the 103 routines see it
    real(R8) :: ortho                            ! tangential / normal metric rows on the face normal
  end type face_geo_type

  type :: int_grid_type
    integer, allocatable :: a(:,:,:)
  end type int_grid_type

contains

  !> Type 104 has no cross-rank exchange yet: the coupling reads the partner cells
  !> of the other solver in place, so every donor of an entry whose block this
  !> rank owns must belong to this rank too. Run once per level at setup, before
  !> Setup_Coupling_Segments; always satisfied in serial.
  subroutine Check_Chimera_Coupling_Donors ( FUSS_domain, ARES_domain, level )
    use FUSS_Mod_MPI, only: FUSS_is_local_block => is_local_block, mpi_abort_all
    use ARES_Mod_MPI, only: ARES_is_local_block => is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(ARES_domain_type), intent(in) :: ARES_domain
    integer, intent(in)                :: level
    ! Local
    integer :: ii, i, c, n_remote
    character(len=256) :: msg

    ! Entries with at least one donor on another rank: FUSS entries read ARES
    ! donors, ARES entries read FUSS donors
    n_remote = 0
    do ii = 1, FUSS_domain % n_local_bc
      i = FUSS_domain % local_bc_idx(ii)
      if (FUSS_domain % bc(i) % type /= 104) cycle
      do c = 1, sum ( FUSS_domain % bc(i) % ni )
        if (.not. ARES_is_local_block ( FUSS_domain % bc(i) % donorID(c,1) )) then
          n_remote = n_remote + 1
          exit
        endif
      enddo
    enddo
    do ii = 1, ARES_domain % n_local_bc
      i = ARES_domain % local_bc_idx(ii)
      if (ARES_domain % bc(i) % type /= 104) cycle
      do c = 1, sum ( ARES_domain % bc(i) % ni )
        if (.not. FUSS_is_local_block ( ARES_domain % bc(i) % donorID(c,1) )) then
          n_remote = n_remote + 1
          exit
        endif
      enddo
    enddo

    if (n_remote > 0) then
      write(msg,'(I0,A,I0,A)') n_remote, ' coupled chimera (104) entries in the level-', level, &
        ' BC files have donors on another MPI rank. BC 104 does not exchange donors between '// &
        'ranks yet: run it serial.'
      call mpi_abort_all(trim(msg))
    endif

  end subroutine Check_Chimera_Coupling_Donors


  !> Build the segments of the 104 interface of one level, once at setup.
  !> The partner faces of a FUSS 104 face are searched among the ARES cells that
  !> ATLAS lists as donors of its first ghost row (the ARES cells next to the
  !> interface overlap that ghost cell); an ARES 104 face of such a cell makes a
  !> segment when it looks at the FUSS face, lies on the same surface and the
  !> two faces overlap. Every 104 face must end up with at least one segment.
  subroutine Setup_Coupling_Segments ( FUSS_domain, ARES_domain, level )
    use ARES_Global_m, only: nprim
    use FUSS_Mod_MPI,  only: mpi_is_root, mpi_abort_all
    use hydraAF_Mod_CouplingExchange, only: K_WALL_UNSET
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(ARES_domain_type), intent(in) :: ARES_domain
    integer, intent(in)                :: level
    ! Local
    type(segment_set_type), pointer  :: set
    type(face_geo_type), allocatable :: gF(:), gA(:)
    type(int_grid_type), allocatable :: head(:)
    integer,  allocatable :: next(:), found(:), pos(:)
    real(R8), allocatable :: covF(:), covA(:)
    integer  :: i, j, e, m, c, sg, nfound, Bs, Is, Js, Ks
    character(len=256) :: msg

    if (level > max_levels) call mpi_abort_all('Setup_Coupling_Segments: more multigrid levels than max_levels')
    set => segs(level)
    call Reset_Set ( set )

    ! 104 faces of each solver
    set%nF = count ( FUSS_domain % bc(1:FUSS_domain%nbound) % type == 104 ) ! Conta le celle di FUSS con bc 104
    set%nA = count ( ARES_domain % bc(1:ARES_domain%nbound) % type == 104 ) ! Conta le celle di ARES con bc 104
    if (set%nF == 0 .and. set%nA == 0) return
    if (set%nF == 0 .or. set%nA == 0) then
      write(msg,'(A,I0,A,I0,A,I0)') 'BC 104 on one solver only at level ', level, ': FUSS faces ', set%nF, &
                                    ', ARES faces ', set%nA
      call mpi_abort_all(trim(msg))
    endif

    allocate ( set%eF(set%nF), set%eA(set%nA) ) ! alloco spazio per i "segmenti" con il numero di celle che ho trovato per ogni faccia 104
    allocate ( set%rowF(FUSS_domain%nbound), set%rowA(ARES_domain%nbound) )
    set%rowF = 0
    set%rowA = 0
    e = 0
    do i = 1, FUSS_domain % nbound
      if (FUSS_domain % bc(i) % type /= 104) cycle
      e = e + 1
      set%eF(e) = i
      set%rowF(i) = e
    enddo
    ! CERCO LE BC 104 tra tutte
    ! i    :  1  2  3  4  5  6 -> facce del blocco
    ! BC type :  1 104 104 7 104 1 -> bc per ogni faccia
    ! eF   = [2, 3, 5]             facce che hanno la bc 104
    ! rowF = [0, 1, 2, 0, 3, 0]    ! salvo

    e = 0
    do j = 1, ARES_domain % nbound
      if (ARES_domain % bc(j) % type /= 104) cycle
      e = e + 1
      set%eA(e) = j
      set%rowA(j) = e
    enddo


    ! Calcolo metria e geometria delle facce con 104
    allocate ( gF(set%nF), gA(set%nA) )
    do e = 1, set%nF
      gF(e) = FUSS_Face_Geometry ( FUSS_domain, set%eF(e) )
    enddo
    do e = 1, set%nA
      gA(e) = ARES_Face_Geometry ( ARES_domain, set%eA(e) )
    enddo

    ! ARES cell -> its 104 faces (up to 6), as a linked list per block
    allocate ( head(ARES_domain%nb), next(set%nA) )
    ! head(m)%a(i,j,k) -> un array 3D per ogni blocco ARES che contiene facce 104 -> contiene il numero e dell'ultima faccia 104 inserita per quella cella; 0 se la cella non ne ha
    ! next(e)	-> dim(nA) (una casella per faccia 104) -> contiene la faccia 104 successiva della stessa cella; 0 se è l'ultima della lista
    do e = 1, set%nA
      j = set%eA(e)
      m = ARES_domain % bc(j) % b
      if (.not. allocated(head(m)%a)) then
        allocate ( head(m)%a(ARES_domain%blk(m)%dim(1), ARES_domain%blk(m)%dim(2), ARES_domain%blk(m)%dim(3)) )
        head(m)%a = 0
      endif
      next(e) = head(m)%a(ARES_domain%bc(j)%i, ARES_domain%bc(j)%j, ARES_domain%bc(j)%k)
      head(m)%a(ARES_domain%bc(j)%i, ARES_domain%bc(j)%j, ARES_domain%bc(j)%k) = e
    enddo

    ! Segments, face by face on the FUSS side
    ! Creazione dei segmenti
    call Grow_Segments ( set, 4*set%nF + 16 )
    allocate ( found(6 * max(1, maxval( [ (FUSS_domain%bc(set%eF(e))%ni(1), e = 1, set%nF) ] ))) )
    do e = 1, set%nF
      i = set%eF(e) ! indice bc della faccia FUSS
      nfound = 0  ! svuota la lista delle facce ARES già provate
      do c = 1, FUSS_domain % bc(i) % ni(1) ! Ciclo sui donatori ATLAS della prima fila di ghost cioè celle ARES
        Bs = FUSS_domain % bc(i) % donorID(c,1)
        Is = FUSS_domain % bc(i) % donorID(c,2)
        Js = FUSS_domain % bc(i) % donorID(c,3)
        Ks = FUSS_domain % bc(i) % donorID(c,4)
        if (Bs < 1 .or. Bs > ARES_domain%nb) cycle ! 	donatore non valido o di un altro dominio
        if (.not. allocated(head(Bs)%a)) cycle     ! quel blocco ARES non ha facce 104, inutile cercare
        if (Is < 1 .or. Js < 1 .or. Ks < 1) cycle   ! indice di ghost o di padding, fuori dal blocco
        if (Is > size(head(Bs)%a,1) .or. Js > size(head(Bs)%a,2) .or. Ks > size(head(Bs)%a,3)) cycle ! oltre le dimensioni del blocco
        m = head(Bs)%a(Is,Js,Ks)          !  prima faccia 104 della cella donatrice (0 = nessuna)
        do while (m > 0)
          if (.not. any(found(1:nfound) == m)) then
            nfound = nfound + 1
            found(nfound) = m
            call Try_Segment ( e, m ) ! ! prova la coppia (faccia FUSS e, faccia ARES m)
          endif
          m = next(m)                   ! faccia successiva della stessa cella
        enddo
      enddo
    enddo

    ! Segments of each face, both sides
    allocate ( set%pF(set%nF+1), set%pA(set%nA+1), set%lF(set%ns), set%lA(set%ns) )
    set%pF = 0
    set%pA = 0
    do sg = 1, set%ns
      set%pF(set%rowF(set%sF(sg))+1) = set%pF(set%rowF(set%sF(sg))+1) + 1
      set%pA(set%rowA(set%sA(sg))+1) = set%pA(set%rowA(set%sA(sg))+1) + 1
    enddo
    set%pF(1) = 1
    set%pA(1) = 1
    do e = 1, set%nF
      set%pF(e+1) = set%pF(e+1) + set%pF(e)
    enddo
    do e = 1, set%nA
      set%pA(e+1) = set%pA(e+1) + set%pA(e)
    enddo
    allocate ( pos(max(set%nF, set%nA)) )
    pos(1:set%nF) = set%pF(1:set%nF)
    do sg = 1, set%ns
      e = set%rowF(set%sF(sg))
      set%lF(pos(e)) = sg
      pos(e) = pos(e) + 1
    enddo
    pos(1:set%nA) = set%pA(1:set%nA)
    do sg = 1, set%ns
      e = set%rowA(set%sA(sg))
      set%lA(pos(e)) = sg
      pos(e) = pos(e) + 1
    enddo

    ! Coverage of each face by its segments: 1 on a closed interface
    allocate ( covF(set%nF), covA(set%nA) )
    covF = 0d0
    covA = 0d0
    do sg = 1, set%ns
      covF(set%rowF(set%sF(sg))) = covF(set%rowF(set%sF(sg))) + set%wF(sg)
      covA(set%rowA(set%sA(sg))) = covA(set%rowA(set%sA(sg))) + set%wA(sg)
    enddo
    if (any(covF <= 0d0) .or. any(covA <= 0d0)) then
      write(msg,'(I0,A,I0,A,I0,A)') count(covF <= 0d0), ' FUSS and ', count(covA <= 0d0), &
        ' ARES 104 faces of level ', level, ' have no partner face on the other solver: check that '// &
        'the two sides of the interface touch and that ATLAS found their donors.'
      call mpi_abort_all(trim(msg))
    endif

    allocate ( set%qF(set%ns), set%fA(nprim, set%ns), set%kS(set%ns) )
    set%qF = 0d0
    set%fA = 0d0
    ! No fluid wall conductivity yet: both wall routines fall back to their local estimate
    set%kS = K_WALL_UNSET
    set%ready = .true.

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' BC 104, level ', level, ': ', set%nF, ' FUSS faces, ', set%nA, &
                                         ' ARES faces, ', set%ns, ' segments'
      write(*,'(A,2F10.6,A,2F10.6)') '   coverage of the faces by their segments: FUSS', minval(covF), maxval(covF), &
                                     '   ARES', minval(covA), maxval(covA)
      if (max(maxval(gF%ortho), maxval(gA%ortho)) > 1d-3) then
        write(*,'(A,ES9.2,A)') '   WARNING: the mesh is not orthogonal at the 104 interface (tangential/normal metric ', &
          max(maxval(gF%ortho), maxval(gA%ortho)), '): this version neglects the partner tangential terms.'
      endif
    endif

  contains

    !> Segment between FUSS face e and ARES face m, if they overlap
    subroutine Try_Segment ( e, m )
      integer, intent(in) :: e, m
      real(R8) :: As, AF, AA, tol

      ! The two faces must look at each other and lie on the same surface
      if (dot_product(gF(e)%n, gA(m)%n) > -0.5d0) return !  Le due facce si guardano?
      tol = 0.05d0 * (gF(e)%d + gA(m)%d)
      if (abs(dot_product(gA(m)%c - gF(e)%c, gF(e)%n)) > tol) return !   Le due facce sono sulla stessa superficie?

      call Overlap_Area ( gA(m)%p, gF(e)%p, gF(e)%c, gF(e)%n, As, AA, AF ) ! Le due facce si sovrappongono
      if (As <= 1d-9 * AF) return

      ! Overlap_Area funziona così:
      !proietta i due quadrilateri sul piano della faccia FUSS, usando una base 2D (e1, e2);
      !li orienta entrambi in senso antiorario;
      !ritaglia il quadrilatero ARES con quello FUSS usando l'algoritmo di Sutherland–Hodgman;
      !calcola con la formula di Gauss (shoelace) l'area dell'intersezione As e le aree proiettate AA e AF.

      if (set%ns == size(set%sF)) call Grow_Segments ( set, 2*size(set%sF) )
      set%ns = set%ns + 1
      set%sF(set%ns) = set%eF(e)
      set%sA(set%ns) = set%eA(m)
      set%wF(set%ns) = As / AF
      set%wA(set%ns) = As / AA
      set%dF(set%ns) = gF(e)%d
      set%dA(set%ns) = gA(m)%d
    end subroutine Try_Segment

  end subroutine Setup_Coupling_Segments


  !> Coupled fluxes of the 104 faces. On every segment the FUSS and ARES estimates
  !> come from the 103 routines with the partner cell; each face then takes,
  !> summed over its segments, 0.5*own - 0.5*partner energy flux, as the 103 pair
  !> of faces does. As in Compute_FUSS_Flux / Compute_ARES_Flux, the FUSS side
  !> reads the fluid wall conductivity of the previous call and the ARES side
  !> reads the same value, then stores the one it has just evaluated.
  subroutine Compute_Coupling_Segment_Fluxes ( FUSS_domain, ARES_domain, level )
    use ARES_Global_m, only: nprim, nh
    use hydraAF_Mod_CouplingExchange, only: K_WALL_UNSET
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    integer, intent(in)                   :: level
    ! Local
    type(segment_set_type), pointer :: set
    integer  :: e, k, sg
    real(R8) :: flux, fluxv(nprim), acc, accv(nprim), accE, wsum, k_new

    if (level > max_levels) return
    set => segs(level)
    if (.not. set%ready) return

    ! FUSS estimate on every segment; a face owns its segments, so one face per thread
    !$omp parallel do schedule (dynamic, 16) private(e, k, sg, flux)
    do e = 1, set%nF
      do k = set%pF(e), set%pF(e+1) - 1
        sg = set%lF(k)
        call Eval_FUSS_Segment ( set, sg, FUSS_domain, ARES_domain, flux )
        set%qF(sg) = flux * set%wF(sg)
      enddo
    enddo
    !$omp end parallel do

    ! ARES estimate on every segment (the ARES routine also sets the RANS ghost of its
    ! face), then the fresh wall conductivity of the segment for the next call
    !$omp parallel do schedule (dynamic, 16) private(e, k, sg, fluxv, k_new)
    do e = 1, set%nA
      do k = set%pA(e), set%pA(e+1) - 1
        sg = set%lA(k)
        k_new = K_WALL_UNSET
        call Eval_ARES_Segment ( set, sg, FUSS_domain, ARES_domain, fluxv, k_new = k_new )
        set%fA(:,sg) = fluxv * set%wA(sg)
        set%kS(sg) = k_new
      enddo
    enddo
    !$omp end parallel do

    ! Exchange as in 103, segment by segment: what one side takes, the other gives
    !$omp parallel do schedule (static) private(e, k, sg, acc)
    do e = 1, set%nF
      acc = 0d0
      do k = set%pF(e), set%pF(e+1) - 1
        sg = set%lF(k)
        acc = acc + 0.5d0 * set%qF(sg) - 0.5d0 * set%fA(nh,sg)
      enddo
      FUSS_domain % bc(set%eF(e)) % ext_flux = acc
    enddo
    !$omp end parallel do

    !$omp parallel do schedule (static) private(e, k, sg, accv, accE, wsum)
    do e = 1, set%nA
      accv = 0d0
      accE = 0d0
      wsum = 0d0
      do k = set%pA(e), set%pA(e+1) - 1
        sg = set%lA(k)
        accv = accv + 0.5d0 * set%fA(:,sg)
        accE = accE + 0.5d0 * set%fA(nh,sg) - 0.5d0 * set%qF(sg)
        wsum = wsum + set%wA(sg)
      enddo
      ! Wall stress and RANS terms: face average, whatever the thermal coverage
      ! (ARES_new_flux carries the factor 2 that the 0.5 cancels, as in 103)
      ARES_domain % bc(set%eA(e)) % ext_flux     = accv / wsum
      ! Energy: sum over the segments; a part of the face with no partner is adiabatic
      ARES_domain % bc(set%eA(e)) % ext_flux(nh) = accE
    enddo
    !$omp end parallel do

  end subroutine Compute_Coupling_Segment_Fluxes


  !> Wall output of a FUSS 104 face: (Tw, qw) of its segments, averaged by area
  subroutine Segment_Wall_Output_FUSS ( FUSS_domain, ARES_domain, level, i, Ovar )
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(in)    :: ARES_domain
    integer, intent(in)                   :: level, i
    real(R8), intent(out)                 :: Ovar(2)
    ! Local
    type(segment_set_type), pointer :: set
    integer  :: e, k, sg
    real(R8) :: flux, O(2), wsum

    Ovar = 0d0
    if (level > max_levels) return
    set => segs(level)
    if (.not. set%ready) return
    e = set%rowF(i)
    if (e == 0) return

    wsum = 0d0
    do k = set%pF(e), set%pF(e+1) - 1
      sg = set%lF(k)
      call Eval_FUSS_Segment ( set, sg, FUSS_domain, ARES_domain, flux, O )
      Ovar = Ovar + set%wF(sg) * O
      wsum = wsum + set%wF(sg)
    enddo
    if (wsum > 0d0) Ovar = Ovar / wsum

  end subroutine Segment_Wall_Output_FUSS


  !> Wall output of an ARES 104 face: the wall variables of its segments, averaged
  !> by area. The stored wall conductivity is read, not advanced, as in AF_IO_Wall.
  subroutine Segment_Wall_Output_ARES ( FUSS_domain, ARES_domain, level, j, Ovar )
    use ARES_Global_m, only: nprim
    implicit none
    type(FUSS_domain_type), intent(in)    :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    integer, intent(in)                   :: level, j
    real(R8), intent(out)                 :: Ovar(8)
    ! Local
    type(segment_set_type), pointer :: set
    integer  :: e, k, sg
    real(R8) :: fluxv(nprim), O(8), wsum

    Ovar = 0d0
    if (level > max_levels) return
    set => segs(level)
    if (.not. set%ready) return
    e = set%rowA(j)
    if (e == 0) return

    wsum = 0d0
    do k = set%pA(e), set%pA(e+1) - 1
      sg = set%lA(k)
      call Eval_ARES_Segment ( set, sg, FUSS_domain, ARES_domain, fluxv, O )
      Ovar = Ovar + set%wA(sg) * O
      wsum = wsum + set%wA(sg)
    enddo
    if (wsum > 0d0) Ovar = Ovar / wsum

  end subroutine Segment_Wall_Output_ARES


  !> FUSS side of segment sg: the 103 routine with the ARES cell of the segment as
  !> partner. new_flux is the flux through the whole FUSS face, as in 103.
  subroutine Eval_FUSS_Segment ( set, sg, FUSS_domain, ARES_domain, new_flux, Ovar )
    use hydraAF_Lib_BC_Coupled_Wall_Temperature, only: FUSS_BC_Coupled_Wall_Temperature
    use ARES_Global_m, only: np, nh
    use FLINT_Lib_Thermodynamic
    implicit none
    type(segment_set_type), intent(in)    :: set
    integer, intent(in)                   :: sg
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(in)    :: ARES_domain
    real(R8), intent(out)                 :: new_flux
    real(R8), optional, intent(out)       :: Ovar(2)
    ! Local
    integer  :: i, j, Bm, Im, Jm, Km, Fm, dir, Bs, Is, Js, Ks
    real(R8) :: P_con(np), Tg(6)
    type(FUSS_tensor_3D_type) :: Mcon

    i  = set%sF(sg)
    j  = set%sA(sg)
    Bm = FUSS_domain % bc(i) % b
    Im = FUSS_domain % bc(i) % i
    Jm = FUSS_domain % bc(i) % j
    Km = FUSS_domain % bc(i) % k
    Fm = FUSS_domain % bc(i) % f
    Bs = ARES_domain % bc(j) % b
    Is = ARES_domain % bc(j) % i
    Js = ARES_domain % bc(j) % j
    Ks = ARES_domain % bc(j) % k

    ! Partner: the ARES cell that owns the other face of the segment
    P_con = ARES_domain % blk(Bs) % P (1:np,Is,Js,Ks)
    Tg = ph2vars( ARES_domain%blk(Bs)%P(np,Is,Js,Ks), ARES_domain%blk(Bs)%P(nh,Is,Js,Ks), T_tab )

    ! Partner metric, normal distance only: the normal row of the own cell rescaled
    ! so that the routine sees the partner cell centre at dA; tangential rows zero
    dir = (Fm+1)/2
    Mcon%c = 0d0
    Mcon%c(dir,:) = FUSS_domain % blk(Bm) % M(Im,Jm,Km) % c(dir,:) * set%dF(sg) / set%dA(sg)

    call FUSS_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                           0, 0, 0, 0, P_con, Mcon, Tg, new_flux, Ovar, &
                                           k_fluid_wall = set%kS(sg) )

  end subroutine Eval_FUSS_Segment


  !> ARES side of segment sg: the 103 routine with the FUSS cell of the segment as
  !> partner. new_flux is the flux through the whole ARES face, as in 103; k_new is
  !> the fluid wall conductivity the routine has just evaluated on the segment.
  subroutine Eval_ARES_Segment ( set, sg, FUSS_domain, ARES_domain, new_flux, Ovar, k_new )
    use hydraAF_Lib_BC_Coupled_Wall_Temperature, only: ARES_BC_Coupled_Wall_Temperature
    use ARES_Global_m, only: nprim
    implicit none
    type(segment_set_type), intent(in)    :: set
    integer, intent(in)                   :: sg
    type(FUSS_domain_type), intent(in)    :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    real(R8), intent(out)                 :: new_flux(nprim)
    real(R8), optional, intent(out)       :: Ovar(8)
    real(R8), optional, intent(out)       :: k_new
    ! Local
    integer  :: i, j, Bm, Im, Jm, Km, Fm, dir, Bs, Is, Js, Ks
    real(R8) :: Pg(1,6), matID_con
    type(ARES_tensor_3D_type) :: Mcon

    j  = set%sA(sg)
    i  = set%sF(sg)
    Bm = ARES_domain % bc(j) % b
    Im = ARES_domain % bc(j) % i
    Jm = ARES_domain % bc(j) % j
    Km = ARES_domain % bc(j) % k
    Fm = ARES_domain % bc(j) % f
    Bs = FUSS_domain % bc(i) % b
    Is = FUSS_domain % bc(i) % i
    Js = FUSS_domain % bc(i) % j
    Ks = FUSS_domain % bc(i) % k

    ! Partner: the FUSS cell that owns the other face of the segment
    Pg = FUSS_domain % blk(Bs) % T (Is,Js,Ks)
    matID_con = FUSS_domain % blk(Bs) % matID (Is,Js,Ks)

    ! Partner metric, normal distance only (see Eval_FUSS_Segment)
    dir = (Fm+1)/2
    Mcon%c = 0d0
    Mcon%c(dir,:) = ARES_domain % blk(Bm) % M(Im,Jm,Km) % c(dir,:) * set%dA(sg) / set%dF(sg)

    call ARES_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, ARES_domain % blk(Bm), &
                                           0, 0, 0, 0, matID_con, Mcon, Pg, new_flux, Ovar, &
                                           k_fluid_wall = set%kS(sg), k_fluid_wall_new = k_new )

  end subroutine Eval_ARES_Segment


  function FUSS_Face_Geometry ( domain, i ) result ( geo )
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in)                :: i
    type(face_geo_type)                :: geo
    ! Local
    integer  :: b, Im, Jm, Km, Fm, a1, a2, a3, fi(3)
    real(R8) :: corner(3,0:1,0:1,0:1)

    b  = domain % bc(i) % b ! blocco
    Im = domain % bc(i) % i ! I
    Jm = domain % bc(i) % j ! J
    Km = domain % bc(i) % k ! K
    Fm = domain % bc(i) % f ! Faccia
    do a3 = 0, 1
      do a2 = 0, 1
        do a1 = 0, 1
          corner(:,a1,a2,a3) = domain % blk(b) % node(Im-1+a1, Jm-1+a2, Km-1+a3) % c ! Sono i vertici della cella
        enddo
      enddo
    enddo
    fi = [Im, Jm, Km]
    fi((Fm+1)/2) = fi((Fm+1)/2) - mod(Fm,2) ! Direzione: dispari verso 1,3,5 pari verso faccia 2,4,6
    geo = Face_Geometry ( corner, Fm, domain % blk(b) % M(Im,Jm,Km) % c, &
                          domain % blk(b) % dir((Fm+1)/2) % f(fi(1),fi(2),fi(3)) % N )

  end function FUSS_Face_Geometry


  function ARES_Face_Geometry ( domain, i ) result ( geo )
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    integer, intent(in)                :: i
    type(face_geo_type)                :: geo
    ! Local
    integer  :: b, Im, Jm, Km, Fm, a1, a2, a3, fi(3)
    real(R8) :: corner(3,0:1,0:1,0:1)

    b  = domain % bc(i) % b
    Im = domain % bc(i) % i
    Jm = domain % bc(i) % j
    Km = domain % bc(i) % k
    Fm = domain % bc(i) % f
    do a3 = 0, 1
      do a2 = 0, 1
        do a1 = 0, 1
          corner(:,a1,a2,a3) = domain % blk(b) % node(Im-1+a1, Jm-1+a2, Km-1+a3) % c
        enddo
      enddo
    enddo
    fi = [Im, Jm, Km]
    fi((Fm+1)/2) = fi((Fm+1)/2) - mod(Fm,2)
    geo = Face_Geometry ( corner, Fm, domain % blk(b) % M(Im,Jm,Km) % c, &
                          domain % blk(b) % dir((Fm+1)/2) % f(fi(1),fi(2),fi(3)) % N )

  end function ARES_Face_Geometry


  !> Geometry of face Fm of a cell from its 8 corners, corner(:,a1,a2,a3) with
  !> a = 0/1 the low/high node in i, j, k, its metric rows M and the stored face
  !> normal N. d is the distance of the cell centre as the 103 routines see it:
  !> they weight the cell by 2*(M(dir,:).N) = 1/d.
  function Face_Geometry ( corner, Fm, M, N ) result ( geo )
    implicit none
    real(R8), intent(in) :: corner(3,0:1,0:1,0:1), M(3,3), N(3)
    integer,  intent(in) :: Fm
    type(face_geo_type)  :: geo
    ! Local
    integer  :: dir, h, t1, t2
    real(R8) :: cc(3), mn

    ! p(3,4)	i 4 nodi della faccia, ordinati lungo il bordo
    ! c(3)	il centro della faccia (media dei nodi)
    ! n(3)	la normale unitaria, uscente dalla cella proprietaria
    ! d	la distanza centro cella–faccia, ricavata dalla metrica come nel 103
    ! ortho	il peso dei termini tangenziali rispetto a quello normale (0 se la mesh è ortogonale)

    dir = (Fm+1)/2        ! Direzione normale alla faccia
    h   = 1 - mod(Fm,2)   ! 0 = faccia bassa/sinistra, 1 = faccia alta/destra
    select case (dir)
      case (1)
        geo%p(:,1) = corner(:,h,0,0)
        geo%p(:,2) = corner(:,h,1,0)
        geo%p(:,3) = corner(:,h,1,1)
        geo%p(:,4) = corner(:,h,0,1)
      case (2)
        geo%p(:,1) = corner(:,0,h,0)
        geo%p(:,2) = corner(:,1,h,0)
        geo%p(:,3) = corner(:,1,h,1)
        geo%p(:,4) = corner(:,0,h,1)
      case (3)
        geo%p(:,1) = corner(:,0,0,h)
        geo%p(:,2) = corner(:,1,0,h)
        geo%p(:,3) = corner(:,1,1,h)
        geo%p(:,4) = corner(:,0,1,h)
    end select
    geo%c = 0.25d0 * sum ( geo%p, dim=2 ) ! centro della faccia
    geo%n = Cross ( geo%p(:,3) - geo%p(:,1), geo%p(:,4) - geo%p(:,2) ) ! Normale unitaria
    geo%n = geo%n / norm2 ( geo%n )
    cc = 0.125d0 * sum ( reshape ( corner, [3,8] ), dim=2 )
    if (dot_product ( geo%n, geo%c - cc ) < 0d0) geo%n = -geo%n

    mn = abs ( dot_product ( M(dir,:), N ) )
    geo%d = 0.5d0 / mn ! Distanza dal centro cella dal centro faccia
    t1 = mod(dir,3) + 1
    t2 = mod(dir+1,3) + 1
    geo%ortho = max ( abs(dot_product(M(t1,:), N)), abs(dot_product(M(t2,:), N)) ) / mn ! indice non ortogonalità

  end function Face_Geometry


  !> Area of the overlap of two convex quads pa and pb, both projected on the plane
  !> through c with normal n (Sutherland-Hodgman clipping of pa by pb), and the
  !> projected areas of the two quads
  subroutine Overlap_Area ( pa, pb, c, n, area, area_a, area_b )
    implicit none
    real(R8), intent(in)  :: pa(3,4), pb(3,4), c(3), n(3)
    real(R8), intent(out) :: area, area_a, area_b
    ! Local
    real(R8) :: e1(3), e2(3), a(2,16), b(2,4), out(2,16), s(2), t(2), p(2), q(2), sp, sq
    integer  :: na, no, k, m

    e1 = pb(:,2) - pb(:,1)
    e1 = e1 - dot_product ( e1, n ) * n
    e1 = e1 / norm2 ( e1 )
    e2 = Cross ( n, e1 )
    do k = 1, 4
      a(:,k) = [ dot_product(pa(:,k) - c, e1), dot_product(pa(:,k) - c, e2) ]
      b(:,k) = [ dot_product(pb(:,k) - c, e1), dot_product(pb(:,k) - c, e2) ]
    enddo
    ! Both counter-clockwise, so that "inside" is on the left of every edge
    if (Shoelace(a, 4) < 0d0) a(:,1:4) = a(:,4:1:-1)
    if (Shoelace(b, 4) < 0d0) b(:,1:4) = b(:,4:1:-1)
    area_a = Shoelace(a, 4)
    area_b = Shoelace(b, 4)

    na = 4
    do m = 1, 4
      s = b(:,m)
      t = b(:,mod(m,4)+1)
      no = 0
      do k = 1, na
        p = a(:,k)
        q = a(:,mod(k,na)+1)
        sp = (t(1)-s(1))*(p(2)-s(2)) - (t(2)-s(2))*(p(1)-s(1))
        sq = (t(1)-s(1))*(q(2)-s(2)) - (t(2)-s(2))*(q(1)-s(1))
        if (sp >= 0d0) then
          no = no + 1
          out(:,no) = p
        endif
        if (sp*sq < 0d0) then
          no = no + 1
          out(:,no) = p + (q - p) * sp / (sp - sq)
        endif
      enddo
      na = no
      if (na < 3) exit
      a(:,1:na) = out(:,1:na)
    enddo

    area = 0d0
    if (na >= 3) area = Shoelace(a, na)

  end subroutine Overlap_Area


  pure function Shoelace ( a, n ) result ( area )
    implicit none
    real(R8), intent(in) :: a(:,:)
    integer, intent(in)  :: n
    real(R8)             :: area
    integer :: k, kn

    area = 0d0
    do k = 1, n
      kn = mod(k,n) + 1
      area = area + a(1,k)*a(2,kn) - a(1,kn)*a(2,k)
    enddo
    area = 0.5d0 * area

  end function Shoelace


  pure function Cross ( u, v ) result ( w )
    implicit none
    real(R8), intent(in) :: u(3), v(3)
    real(R8)             :: w(3)

    w = [ u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1) ]

  end function Cross


  subroutine Grow_Segments ( set, n )
    implicit none
    type(segment_set_type), intent(inout) :: set
    integer, intent(in)                   :: n
    ! Local
    integer,  allocatable :: itmp(:)
    real(R8), allocatable :: rtmp(:)

    if (.not. allocated(set%sF)) then
      allocate ( set%sF(n), set%sA(n), set%wF(n), set%wA(n), set%dF(n), set%dA(n) )
      return
    endif
    allocate ( itmp(n) ); itmp(1:set%ns) = set%sF(1:set%ns); call move_alloc ( itmp, set%sF )
    allocate ( itmp(n) ); itmp(1:set%ns) = set%sA(1:set%ns); call move_alloc ( itmp, set%sA )
    allocate ( rtmp(n) ); rtmp(1:set%ns) = set%wF(1:set%ns); call move_alloc ( rtmp, set%wF )
    allocate ( rtmp(n) ); rtmp(1:set%ns) = set%wA(1:set%ns); call move_alloc ( rtmp, set%wA )
    allocate ( rtmp(n) ); rtmp(1:set%ns) = set%dF(1:set%ns); call move_alloc ( rtmp, set%dF )
    allocate ( rtmp(n) ); rtmp(1:set%ns) = set%dA(1:set%ns); call move_alloc ( rtmp, set%dA )

  end subroutine Grow_Segments


  subroutine Reset_Set ( set )
    implicit none
    type(segment_set_type), intent(inout) :: set

    set%ready = .false.
    set%ns = 0
    set%nF = 0
    set%nA = 0
    if (allocated(set%sF))   deallocate ( set%sF, set%sA, set%wF, set%wA, set%dF, set%dA )
    if (allocated(set%eF))   deallocate ( set%eF )
    if (allocated(set%eA))   deallocate ( set%eA )
    if (allocated(set%rowF)) deallocate ( set%rowF )
    if (allocated(set%rowA)) deallocate ( set%rowA )
    if (allocated(set%pF))   deallocate ( set%pF, set%lF )
    if (allocated(set%pA))   deallocate ( set%pA, set%lA )
    if (allocated(set%qF))   deallocate ( set%qF )
    if (allocated(set%fA))   deallocate ( set%fA )
    if (allocated(set%kS))   deallocate ( set%kS )

  end subroutine Reset_Set

end module hydraAF_Lib_Coupling_Segments
