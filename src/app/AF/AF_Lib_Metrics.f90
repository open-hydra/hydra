!> Metric of the coupled walls (103, 104) of hydra AF.
!>
!> The interface temperature Tint is the root of k_s grad(T_s).N = k_f grad(T_f).N, and each
!> side enters that balance only through the normal component of its one-sided wall gradient:
!>
!>   grad(T).N = modfm2 * A * (Tw - Teff)
!>
!> with A = 1/d = 2 M(dir,:).N (d the distance of the cell centre from the face along N) and
!> Teff the cell temperature corrected by the tangential terms, i.e. the value at the foot of
!> the normal through the face centre:
!>
!>   Teff = T - modfm2 * sum_t dT_t M(t,:).N / A
!>
!> (dT_t the difference of the cell values along the tangential index direction t, per index).
!> Nothing in it depends on how the two meshes are oriented with respect to each other, so each
!> side computes its pair (A, Teff) in its own index frame, with its own metric and its own
!> neighbours, on the rank that owns the cell, and hands it to the other side: Tint is
!>
!>   Tint = ( k_1 A_1 Teff_1 + k_2 A_2 Teff_2 ) / ( k_1 A_1 + k_2 A_2 )
!>
!> on both sides and on conforming (103) and non-conforming (104) interfaces alike. On a 104
!> segment both Teff are moved from their face centre to the centroid x_s of the segment with
!> the gradient of the side in the plane of its face: Teff + g_par . (x_s - x_f). g_par is the
!> part in the plane of the face of the gradient of the cell, from the differences to its
!> neighbours along the wall (inside the block, one-sided at its edges) and to the next cell
!> inwards, with the metric of the same differences of the positions: exact for a linear field
!> on any mesh. The normal difference must enter: where the centres of the neighbours along the
!> wall lie at different distances from it, their differences also carry the normal gradient,
!> large next to a wall, which would otherwise be read as a variation along the wall.
!>
!> The metric M of each side follows diffusive-metric of its own solver (FUSS for the solid side, ARES for the fluid
!> side): each side computes its half alone, so the two keys need not agree:
!>  - mean: the cell metric (the inverse of the mean edges of the cell) and central differences
!>    with the ghost cells along the wall. This is the former coupling metric.
!>  - stencil: the metric of the same differences of the positions as those of the values,
!>    J = [2 modfm2 (x_f - x_c), (x_hi - x_lo)/(hi - lo) for each tangential direction],
!>    M = J^-1, so that D = J^T g and Teff, A are exact for a linear field on any mesh (skewed,
!>    stretched). Along the wall the differences stay inside the block, one-sided at the edges
!>    of the block (the positions of the ghost cells of the other blocks are not known here);
!>    an index direction with one cell, or the third one of a 2-D mesh, keeps the cell edge and
!>    no difference. On the solid side it is the wall metric of FUSS (Wall_Face_Metric, one-sided variant), which
!>    keeps the cell metric when the face has no tangential difference at all; the fluid side replaces the normal
!>    column in that case too, as the wall metric of ARES (Wall_Face_Stencils).
!>    On the fluid side the velocity and the turbulence variables, whose tangential differences
!>    are zero at the wall, take the metric Mu whose tangential columns are the edges of the
!>    cell projected on the plane of the face: their gradient is then normal to the wall,
!>    u_c / d_n along N (d_n the normal distance of the centre from the face).
module hydraAF_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Interface_Metric, Interface_Stencil, Interface_Stencil_ARES
  public :: FUSS_Wall_Half, ARES_Wall_Half
  public :: Compute_Modfm, Face_Index

  !> Reals of the half of the balance a side hands to the other one: A, Teff, g_par(1:3)
  integer, parameter, public :: NHALF = 5

contains

  !> Report the metric of the coupled walls (once, on the root rank).
  subroutine Setup_Interface_Metric ()
    use FUSS_Mod_MPI, only: mpi_is_root
    implicit none

    if ( .not. mpi_is_root ) return
    write(*,'(A,T35,A)') '   Coupled wall metric (103, 104)', &
      'solid '//trim(merge('stencil','mean   ',Interface_Stencil()))//', fluid '// &
      trim(merge('stencil','mean   ',Interface_Stencil_ARES()))

  end subroutine Setup_Interface_Metric


  !> .true. when the solid side of the coupled walls takes the stencil metric: diffusive-metric = stencil in FUSS.
  logical function Interface_Stencil ()
    use FUSS_Config_Types_m, only: FUSS_obj_space_scheme => obj_space_scheme
    implicit none

    Interface_Stencil = FUSS_obj_space_scheme % stencil_metric

  end function Interface_Stencil


  !> .true. when the fluid side of the coupled walls takes the stencil metric: diffusive-metric = stencil in ARES.
  logical function Interface_Stencil_ARES ()
    use ARES_Config_Types_m, only: ARES_obj_space_scheme => obj_space_scheme
    implicit none

    Interface_Stencil_ARES = ARES_obj_space_scheme % stencil_metric

  end function Interface_Stencil_ARES


  !> Half of the interface balance of face Fm of FUSS cell (Im,Jm,Km): A = 1/d, Teff and, if
  !> asked, the gradient of T in the plane of the face (see the module header).
  subroutine FUSS_Wall_Half ( blk, Im, Jm, Km, Fm, A, Teff, gpar )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    use FUSS_Global_m,         only: ndir
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    integer, intent(in)               :: Im, Jm, Km, Fm
    real(R8), intent(out)             :: A, Teff
    real(R8), intent(out), optional   :: gpar(3)
    ! Local
    integer  :: dir, fi, fj, fk, d, c(3), lo(3), hi(3)
    logical  :: stencil, has(3)
    real(R8) :: dT(3), xc(3), xf(3), xt(3,3)

    stencil = Interface_Stencil()
    call Face_Index ( Fm, dir, Im, Jm, Km, fi, fj, fk )
    c   = [Im, Jm, Km]
    has = .false.
    dT  = 0d0
    xc  = 0d0
    xf  = 0d0
    xt  = 0d0
    if ( stencil ) then
      xc = FUSS_Centre ( blk, Im, Jm, Km )
      xf = FUSS_Face_Centre ( blk, Im, Jm, Km, Fm )
    end if

    do d = 1, 3
      if ( d == dir .or. d > ndir ) cycle
      lo = c
      hi = c
      lo(d) = c(d) - 1
      hi(d) = c(d) + 1
      if ( stencil ) then
        lo(d) = max ( lo(d), 1 )
        hi(d) = min ( hi(d), blk % dim(d) )
        if ( hi(d) == lo(d) ) cycle
        xt(:,d) = ( FUSS_Centre ( blk, hi(1), hi(2), hi(3) ) - FUSS_Centre ( blk, lo(1), lo(2), lo(3) ) ) &
                / real ( hi(d) - lo(d), R8 )
      end if
      has(d) = .true.
      dT(d)  = ( blk % T(hi(1),hi(2),hi(3)) - blk % T(lo(1),lo(2),lo(3)) ) / real ( hi(d) - lo(d), R8 )
    enddo

    call Wall_Half_Core ( dir, Fm, stencil, blk % M(Im,Jm,Km) % c, blk % dir(dir) % f(fi,fj,fk) % N, &
                          blk % T(Im,Jm,Km), has, dT, xc, xf, xt, A, Teff, cell_if_none=.true. )

    ! Gradient of T in the plane of the face (104 only). Along the wall: with stencil the
    ! differences just taken (inside the block, one-sided at its edges); with mean those took
    ! the ghost cells, whose positions are not known here, so they are taken again inside the block
    if ( present(gpar) ) then
      if ( .not. stencil ) then
        xc  = FUSS_Centre ( blk, Im, Jm, Km )
        has = .false.
        dT  = 0d0
        do d = 1, 3
          if ( d == dir .or. d > ndir ) cycle
          lo = c
          hi = c
          lo(d) = max ( c(d) - 1, 1 )
          hi(d) = min ( c(d) + 1, blk % dim(d) )
          if ( hi(d) == lo(d) ) cycle
          xt(:,d) = ( FUSS_Centre ( blk, hi(1), hi(2), hi(3) ) - FUSS_Centre ( blk, lo(1), lo(2), lo(3) ) ) &
                  / real ( hi(d) - lo(d), R8 )
          dT(d)   = ( blk % T(hi(1),hi(2),hi(3)) - blk % T(lo(1),lo(2),lo(3)) ) / real ( hi(d) - lo(d), R8 )
          has(d)  = .true.
        enddo
      end if
      ! Across the face: the next cell inwards
      if ( blk % dim(dir) > 1 ) then
        hi = c
        hi(dir) = c(dir) + 2*mod(Fm,2) - 1
        xt(:,dir) = FUSS_Centre ( blk, hi(1), hi(2), hi(3) ) - xc
        dT(dir)   = blk % T(hi(1),hi(2),hi(3)) - blk % T(Im,Jm,Km)
        has(dir)  = .true.
      end if
      call In_Plane_Gradient ( blk % M(Im,Jm,Km) % c, blk % dir(dir) % f(fi,fj,fk) % N, has, dT, xt, gpar )
    end if

  end subroutine FUSS_Wall_Half


  !> Half of the interface balance of face Fm of ARES cell (Im,Jm,Km), as FUSS_Wall_Half, plus,
  !> if asked, the metric Mu of the gradients of the velocity and of the turbulence variables.
  subroutine ARES_Wall_Half ( blk, Im, Jm, Km, Fm, A, Teff, gpar, Mu )
    use ARES_Advanced_Types_m,   only: ARES_block_type
    use ARES_Global_m,           only: ndir, np, nh
    use FLINT_Lib_Thermodynamic, only: ph2vars, T_tab
    implicit none
    type(ARES_block_type), intent(in) :: blk
    integer, intent(in)               :: Im, Jm, Km, Fm
    real(R8), intent(out)             :: A, Teff
    real(R8), intent(out), optional   :: gpar(3), Mu(3,3)
    ! Local
    integer  :: dir, fi, fj, fk, d, c(3), lo(3), hi(3)
    logical  :: stencil, has(3)
    real(R8) :: dT(3), xc(3), xf(3), xt(3,3), Tc, Tlo, Thi

    stencil = Interface_Stencil_ARES()
    call Face_Index ( Fm, dir, Im, Jm, Km, fi, fj, fk )
    c   = [Im, Jm, Km]
    has = .false.
    dT  = 0d0
    xc  = 0d0
    xf  = 0d0
    xt  = 0d0
    if ( stencil ) then
      xc = ARES_Centre ( blk, Im, Jm, Km )
      xf = ARES_Face_Centre ( blk, Im, Jm, Km, Fm )
    end if
    Tc = ph2vars ( blk % P(np,Im,Jm,Km), blk % P(nh,Im,Jm,Km), T_tab )

    do d = 1, 3
      if ( d == dir .or. d > ndir ) cycle
      lo = c
      hi = c
      lo(d) = c(d) - 1
      hi(d) = c(d) + 1
      if ( stencil ) then
        lo(d) = max ( lo(d), 1 )
        hi(d) = min ( hi(d), blk % dim(d) )
        if ( hi(d) == lo(d) ) cycle
        xt(:,d) = ( ARES_Centre ( blk, hi(1), hi(2), hi(3) ) - ARES_Centre ( blk, lo(1), lo(2), lo(3) ) ) &
                / real ( hi(d) - lo(d), R8 )
      end if
      has(d) = .true.
      Thi = ph2vars ( blk % P(np,hi(1),hi(2),hi(3)), blk % P(nh,hi(1),hi(2),hi(3)), T_tab )
      Tlo = ph2vars ( blk % P(np,lo(1),lo(2),lo(3)), blk % P(nh,lo(1),lo(2),lo(3)), T_tab )
      dT(d) = ( Thi - Tlo ) / real ( hi(d) - lo(d), R8 )
    enddo

    call Wall_Half_Core ( dir, Fm, stencil, blk % M(Im,Jm,Km) % c, blk % dir(dir) % f(fi,fj,fk) % N, &
                          Tc, has, dT, xc, xf, xt, A, Teff, Mu )

    ! Gradient of T in the plane of the face (104 only), as FUSS_Wall_Half
    if ( present(gpar) ) then
      if ( .not. stencil ) then
        xc  = ARES_Centre ( blk, Im, Jm, Km )
        has = .false.
        dT  = 0d0
        do d = 1, 3
          if ( d == dir .or. d > ndir ) cycle
          lo = c
          hi = c
          lo(d) = max ( c(d) - 1, 1 )
          hi(d) = min ( c(d) + 1, blk % dim(d) )
          if ( hi(d) == lo(d) ) cycle
          xt(:,d) = ( ARES_Centre ( blk, hi(1), hi(2), hi(3) ) - ARES_Centre ( blk, lo(1), lo(2), lo(3) ) ) &
                  / real ( hi(d) - lo(d), R8 )
          Thi = ph2vars ( blk % P(np,hi(1),hi(2),hi(3)), blk % P(nh,hi(1),hi(2),hi(3)), T_tab )
          Tlo = ph2vars ( blk % P(np,lo(1),lo(2),lo(3)), blk % P(nh,lo(1),lo(2),lo(3)), T_tab )
          dT(d)  = ( Thi - Tlo ) / real ( hi(d) - lo(d), R8 )
          has(d) = .true.
        enddo
      end if
      ! Across the face: the next cell inwards
      if ( blk % dim(dir) > 1 ) then
        hi = c
        hi(dir) = c(dir) + 2*mod(Fm,2) - 1
        xt(:,dir) = ARES_Centre ( blk, hi(1), hi(2), hi(3) ) - xc
        Thi       = ph2vars ( blk % P(np,hi(1),hi(2),hi(3)), blk % P(nh,hi(1),hi(2),hi(3)), T_tab )
        dT(dir)   = Thi - Tc
        has(dir)  = .true.
      end if
      call In_Plane_Gradient ( blk % M(Im,Jm,Km) % c, blk % dir(dir) % f(fi,fj,fk) % N, has, dT, xt, gpar )
    end if

  end subroutine ARES_Wall_Half


  !> The part of FUSS_Wall_Half and ARES_Wall_Half that does not depend on the solver.
  !> M: cell metric (rows grad(xi_a)), N: unit normal of the face (growing index), Tc: cell
  !> temperature, has(d)/dT(d): tangential differences per index, xc/xf/xt: centre of the cell,
  !> of the face and differences of the centres per index (stencil only). cell_if_none: a face with no tangential
  !> difference keeps the cell metric (the wall metric of FUSS); without it the normal column is still the position
  !> difference across the half cell (the wall metric of ARES).
  subroutine Wall_Half_Core ( dir, Fm, stencil, M, N, Tc, has, dT, xc, xf, xt, A, Teff, Mu, cell_if_none )
    implicit none
    integer, intent(in)             :: dir, Fm
    logical, intent(in)             :: stencil, has(3)
    real(R8), intent(in)            :: M(3,3), N(3), Tc, dT(3), xc(3), xf(3), xt(3,3)
    real(R8), intent(out)           :: A, Teff
    real(R8), intent(out), optional :: Mu(3,3)
    logical, intent(in), optional   :: cell_if_none
    ! Local
    integer  :: d, modfm2
    logical  :: ok, keep
    real(R8) :: Ac(3,3), J(3,3), Mw(3,3), Minv(3,3), C

    modfm2 = 1 - 2*mod(Fm,2)

    ! Columns of Ac: the edges of the cell (M = Ac^-1)
    call Inverse3 ( M, Ac, ok )

    ! Metric of the temperature gradient
    Mw = M
    J  = Ac
    keep = .false.
    if ( present(cell_if_none) ) keep = cell_if_none .and. .not. any ( has )
    if ( stencil .and. .not. keep ) then
      J(:,dir) = 2d0 * modfm2 * ( xf - xc )
      do d = 1, 3
        if ( d /= dir .and. has(d) ) J(:,d) = xt(:,d)
      enddo
      call Inverse3 ( J, Minv, ok )
      if ( ok ) Mw = Minv
    end if

    A = 2d0 * dot_product ( Mw(dir,:), N )
    C = 0d0
    do d = 1, 3
      if ( d /= dir .and. has(d) ) C = C + dT(d) * dot_product ( Mw(d,:), N )
    enddo
    Teff = Tc - modfm2 * C / A

    ! Metric of the quantities whose tangential differences are zero at the wall (velocity,
    ! turbulence variables): tangential columns in the plane of the face
    if ( present(Mu) ) then
      Mu = M
      if ( stencil ) then
        J(:,dir) = 2d0 * modfm2 * ( xf - xc )
        do d = 1, 3
          if ( d == dir ) cycle
          J(:,d) = Ac(:,d) - dot_product ( Ac(:,d), N ) * N
        enddo
        call Inverse3 ( J, Minv, ok )
        if ( ok ) Mu = Minv
      end if
    end if

  end subroutine Wall_Half_Core


  !> Part in the plane of the face (unit normal N) of the gradient of a cell: J g = D, with the
  !> differences D(d) of the values and the columns J(:,d) of the same differences of the positions
  !> (has(d)); a direction without a difference takes the edge of the cell (from the cell metric M)
  !> and D = 0. Exact for a linear field.
  subroutine In_Plane_Gradient ( M, N, has, D, dx, gpar )
    implicit none
    real(R8), intent(in)  :: M(3,3), N(3), D(3), dx(3,3)
    logical, intent(in)   :: has(3)
    real(R8), intent(out) :: gpar(3)
    ! Local (the index is q: d would be the argument D, the names are case-insensitive)
    integer  :: q
    logical  :: ok
    real(R8) :: Ac(3,3), J(3,3), Minv(3,3), g(3)

    gpar = 0d0
    call Inverse3 ( M, Ac, ok )
    J = Ac
    do q = 1, 3
      if ( has(q) ) J(:,q) = dx(:,q)
    enddo
    call Inverse3 ( J, Minv, ok )
    if ( .not. ok ) return
    g    = matmul ( D, Minv )                 ! g_j = sum_q D(q) (J^-1)(q,j)
    gpar = g - dot_product ( g, N ) * N

  end subroutine In_Plane_Gradient


  !> Inverse of a 3x3 matrix; ok = .false. (and Minv = 0) when it is singular relative to the
  !> size of its columns.
  pure subroutine Inverse3 ( J, Minv, ok )
    implicit none
    real(R8), intent(in)  :: J(3,3)
    real(R8), intent(out) :: Minv(3,3)
    logical, intent(out)  :: ok
    ! Local
    real(R8) :: det, scale

    Minv = 0d0
    det = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
        - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
        + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
    scale = norm2(J(:,1)) * norm2(J(:,2)) * norm2(J(:,3))
    ok = ( abs(det) > 1d-12 * scale )
    if ( .not. ok ) return

    Minv(1,1) = ( J(2,2)*J(3,3) - J(2,3)*J(3,2) ) / det
    Minv(1,2) = ( J(1,3)*J(3,2) - J(1,2)*J(3,3) ) / det
    Minv(1,3) = ( J(1,2)*J(2,3) - J(1,3)*J(2,2) ) / det
    Minv(2,1) = ( J(2,3)*J(3,1) - J(2,1)*J(3,3) ) / det
    Minv(2,2) = ( J(1,1)*J(3,3) - J(1,3)*J(3,1) ) / det
    Minv(2,3) = ( J(1,3)*J(2,1) - J(1,1)*J(2,3) ) / det
    Minv(3,1) = ( J(2,1)*J(3,2) - J(2,2)*J(3,1) ) / det
    Minv(3,2) = ( J(1,2)*J(3,1) - J(1,1)*J(3,2) ) / det
    Minv(3,3) = ( J(1,1)*J(2,2) - J(1,2)*J(2,1) ) / det

  end subroutine Inverse3


  !> Centre of a FUSS / ARES cell (mean of its eight nodes) and of face Fm of a cell (mean of
  !> its four nodes), as the stencil metrics of the solvers define them.
  pure function FUSS_Centre ( blk, i, j, k ) result ( x )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    integer, intent(in)               :: i, j, k
    real(R8)                          :: x(3)

    x = 0.125d0 * ( blk % node(i-1,j-1,k-1) % c + blk % node(i,j-1,k-1) % c &
                  + blk % node(i-1,j  ,k-1) % c + blk % node(i,j  ,k-1) % c &
                  + blk % node(i-1,j-1,k  ) % c + blk % node(i,j-1,k  ) % c &
                  + blk % node(i-1,j  ,k  ) % c + blk % node(i,j  ,k  ) % c )

  end function FUSS_Centre

  pure function FUSS_Face_Centre ( blk, i, j, k, Fm ) result ( x )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    integer, intent(in)               :: i, j, k, Fm
    real(R8)                          :: x(3)
    ! Local
    integer :: lo(3), a1, a2, n(3), t(2)

    call Face_Nodes ( i, j, k, Fm, lo, t )
    x = 0d0
    do a2 = 0, 1
      do a1 = 0, 1
        n = lo
        n(t(1)) = n(t(1)) + a1
        n(t(2)) = n(t(2)) + a2
        x = x + blk % node(n(1),n(2),n(3)) % c
      enddo
    enddo
    x = 0.25d0 * x

  end function FUSS_Face_Centre

  pure function ARES_Centre ( blk, i, j, k ) result ( x )
    use ARES_Advanced_Types_m, only: ARES_block_type
    implicit none
    type(ARES_block_type), intent(in) :: blk
    integer, intent(in)               :: i, j, k
    real(R8)                          :: x(3)

    x = 0.125d0 * ( blk % node(i-1,j-1,k-1) % c + blk % node(i,j-1,k-1) % c &
                  + blk % node(i-1,j  ,k-1) % c + blk % node(i,j  ,k-1) % c &
                  + blk % node(i-1,j-1,k  ) % c + blk % node(i,j-1,k  ) % c &
                  + blk % node(i-1,j  ,k  ) % c + blk % node(i,j  ,k  ) % c )

  end function ARES_Centre

  pure function ARES_Face_Centre ( blk, i, j, k, Fm ) result ( x )
    use ARES_Advanced_Types_m, only: ARES_block_type
    implicit none
    type(ARES_block_type), intent(in) :: blk
    integer, intent(in)               :: i, j, k, Fm
    real(R8)                          :: x(3)
    ! Local
    integer :: lo(3), a1, a2, n(3), t(2)

    call Face_Nodes ( i, j, k, Fm, lo, t )
    x = 0d0
    do a2 = 0, 1
      do a1 = 0, 1
        n = lo
        n(t(1)) = n(t(1)) + a1
        n(t(2)) = n(t(2)) + a2
        x = x + blk % node(n(1),n(2),n(3)) % c
      enddo
    enddo
    x = 0.25d0 * x

  end function ARES_Face_Centre


  !> Lowest node of face Fm of cell (i,j,k) and its two tangential index directions: the four
  !> nodes are lo, lo + e_t1, lo + e_t2, lo + e_t1 + e_t2.
  pure subroutine Face_Nodes ( i, j, k, Fm, lo, t )
    implicit none
    integer, intent(in)  :: i, j, k, Fm
    integer, intent(out) :: lo(3), t(2)
    ! Local
    integer :: dir

    dir = (Fm + 1) / 2
    lo  = [i-1, j-1, k-1]
    lo(dir) = lo(dir) + 1 - mod(Fm,2)        ! low face: node i-1, high face: node i
    select case ( dir )
      case (1) ; t = [2, 3]
      case (2) ; t = [1, 3]
      case default ; t = [1, 2]
    end select

  end subroutine Face_Nodes


  pure subroutine Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    implicit none
    integer, intent(in)  :: Fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1-modFm
    modfm2 = 1-2*modFm
    modfm3 = 2*( 2*modFm-1 )

  end subroutine Compute_Modfm


  pure subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
    implicit none
    integer, intent(in)  :: face, cell_i, cell_j, cell_k
    integer, intent(out) :: dir, face_i, face_j, face_k

    select case ( face )
      case(1:2)
        dir = 1
        face_i = cell_i - mod ( face, 2 )
        face_j = cell_j
        face_k = cell_k
      case(3:4)
        dir = 2
        face_i = cell_i
        face_j = cell_j - mod ( face, 2 )
        face_k = cell_k
      case default
        dir = 3
        face_i = cell_i
        face_j = cell_j
        face_k = cell_k - mod ( face, 2 )
    end select

  end subroutine Face_Index

end module hydraAF_Lib_Metrics
