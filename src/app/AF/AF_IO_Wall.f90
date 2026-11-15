module hydraAF_IO_Wall
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use hydra_Config_Types_m, only: hydra_obj_io
  use ARES_IO_Wall, only: ARES_obj_io_bc => obj_io_bc
  use FUSS_IO_Wall, only: FUSS_obj_io_bc => obj_io_bc
  use ARES_Global_m
  use FUSS_Global_m
  use Lib_ORION_data
  use FUSS_Mod_MPI, only: mpi_is_root, mpi_reduce_sum_r8_array
  implicit none
  type(ORION_data), allocatable :: IOwall(:)

contains

  subroutine Initialize_Wall_File ( FUSS_domain, ARES_domain, extension, level )
    use ARES_Advanced_Types_m
    use FUSS_Advanced_Types_m
    use IR_Precision
    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(ARES_domain_type), intent(in) :: ARES_domain
    character(len=*), intent(in)       :: extension
    integer, intent(in)                :: level
    ! Local
    integer :: b, f, bb, i, j, k

    IOwall%tec%extension = extension
    hydra_obj_io%Onwall = 7
    hydra_obj_io%Owallnames = '"y+" "tauX" "tauY" "tauZ" "pw" "Tw" "qw"'

    allocate(IOwall(level)%block(1:count(FUSS_obj_io_bc%coupling_flag)+count(ARES_obj_io_bc%coupling_flag)))

    bb = 0
    
    do b = 1, FUSS_domain % nb; do f = 1, 6
      if (FUSS_obj_io_bc%coupling_flag( b , f )) then
        bb = bb + 1
        IOwall(level)%block(bb)%name = trim(FUSS_phase_prefix)//'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:0,0:FUSS_domain%blk(b)%dim(2),0:FUSS_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1,1:FUSS_domain%blk(b)%dim(2),1:FUSS_domain%blk(b)%dim(3)))
          do k = 0, FUSS_domain%blk(b)%dim(3); do j = 0, FUSS_domain%blk(b)%dim(2)
              IOwall(level) % block(bb) % mesh(1,0,j,k) = FUSS_domain % blk(b) % node(0,j,k) % c(1)
              IOwall(level) % block(bb) % mesh(2,0,j,k) = FUSS_domain % blk(b) % node(0,j,k) % c(2)
              IOwall(level) % block(bb) % mesh(3,0,j,k) = FUSS_domain % blk(b) % node(0,j,k) % c(3)
              ! Section to fool compiler...
              if (IOwall(level)%block(bb)%mesh(1,0,j,k)==0d0 .and. &
                  IOwall(level)%block(bb)%mesh(2,0,j,k)==0d0 .and. &
                  IOwall(level)%block(bb)%mesh(3,0,j,k)==0d0) &
              print*, IOwall(level) % block(bb) % mesh(:,0,j,k)
          enddo; enddo
        case(2)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:0,0:FUSS_domain%blk(b)%dim(2),0:FUSS_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1,1:FUSS_domain%blk(b)%dim(2),1:FUSS_domain%blk(b)%dim(3)))
          do k = 0, FUSS_domain%blk(b)%dim(3); do j = 0, FUSS_domain%blk(b)%dim(2)
              IOwall(level) % block(bb) % mesh(:,0,j,k) = FUSS_domain % blk(b) % node(FUSS_domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:0,0:FUSS_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:FUSS_domain%blk(b)%dim(1),1,1:FUSS_domain%blk(b)%dim(3)))
          do k = 0, FUSS_domain%blk(b)%dim(3); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = FUSS_domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:0,0:FUSS_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:FUSS_domain%blk(b)%dim(1),1,1:FUSS_domain%blk(b)%dim(3)))
          do k = 0, FUSS_domain%blk(b)%dim(3); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = FUSS_domain % blk(b) % node(i,FUSS_domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:FUSS_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:FUSS_domain%blk(b)%dim(1),1:FUSS_domain%blk(b)%dim(2),1))
          do j = 0, FUSS_domain%blk(b)%dim(2); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = FUSS_domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:FUSS_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:FUSS_domain%blk(b)%dim(1),1:FUSS_domain%blk(b)%dim(2),1))
          do j = 0, FUSS_domain%blk(b)%dim(2); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = FUSS_domain % blk(b) % node(i,j,FUSS_domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
      endif
    enddo; enddo

    do b = 1, ARES_domain % nb; do f = 1, 6
      if (ARES_obj_io_bc%coupling_flag( b , f )) then
        bb = bb + 1
        IOwall(level)%block(bb)%name = trim(ARES_phase_prefix)//'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:0,0:ARES_domain%blk(b)%dim(2),0:ARES_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1,1:ARES_domain%blk(b)%dim(2),1:ARES_domain%blk(b)%dim(3)))
          do k = 0, ARES_domain%blk(b)%dim(3); do j = 0, ARES_domain%blk(b)%dim(2)
              IOwall(level) % block(bb) % mesh(1,0,j,k) = ARES_domain % blk(b) % node(0,j,k) % c(1)
              IOwall(level) % block(bb) % mesh(2,0,j,k) = ARES_domain % blk(b) % node(0,j,k) % c(2)
              IOwall(level) % block(bb) % mesh(3,0,j,k) = ARES_domain % blk(b) % node(0,j,k) % c(3)
              ! Section to fool compiler...
              if (IOwall(level)%block(bb)%mesh(1,0,j,k)==0d0 .and. &
                  IOwall(level)%block(bb)%mesh(2,0,j,k)==0d0 .and. &
                  IOwall(level)%block(bb)%mesh(3,0,j,k)==0d0) &
              print*, IOwall(level) % block(bb) % mesh(:,0,j,k)
          enddo; enddo
        case(2)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:0,0:ARES_domain%blk(b)%dim(2),0:ARES_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1,1:ARES_domain%blk(b)%dim(2),1:ARES_domain%blk(b)%dim(3)))
          do k = 0, ARES_domain%blk(b)%dim(3); do j = 0, ARES_domain%blk(b)%dim(2)
              IOwall(level) % block(bb) % mesh(:,0,j,k) = ARES_domain % blk(b) % node(ARES_domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:ARES_domain%blk(b)%dim(1),0:0,0:ARES_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:ARES_domain%blk(b)%dim(1),1,1:ARES_domain%blk(b)%dim(3)))
          do k = 0, ARES_domain%blk(b)%dim(3); do i = 0, ARES_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = ARES_domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:ARES_domain%blk(b)%dim(1),0:0,0:ARES_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:ARES_domain%blk(b)%dim(1),1,1:ARES_domain%blk(b)%dim(3)))
          do k = 0, ARES_domain%blk(b)%dim(3); do i = 0, ARES_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = ARES_domain % blk(b) % node(i,ARES_domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:ARES_domain%blk(b)%dim(1),0:ARES_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:ARES_domain%blk(b)%dim(1),1:ARES_domain%blk(b)%dim(2),1))
          do j = 0, ARES_domain%blk(b)%dim(2); do i = 0, ARES_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = ARES_domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:ARES_domain%blk(b)%dim(1),0:ARES_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:hydra_obj_io%Onwall,1:ARES_domain%blk(b)%dim(1),1:ARES_domain%blk(b)%dim(2),1))
          do j = 0, ARES_domain%blk(b)%dim(2); do i = 0, ARES_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = ARES_domain % blk(b) % node(i,j,ARES_domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
      endif
    enddo; enddo

  end subroutine Initialize_Wall_File



  subroutine Write_Wall_Solution( FUSS_domain, ARES_domain, file, level )
    use Lib_VTK
    use Lib_Tecplot
    use hydraAF_Lib_Ghost, only: Fill_Coupled_Ghost_Cell
    use hydraAF_Mod_CouplingExchange, only: exchange_coupling_Tg_P, exchange_coupling_Pg_matID
    use hydra_Parameters_m
    use ARES_Advanced_Types_m
    use FUSS_Advanced_Types_m
    use strings, only: parse
    implicit none
    type(FUSS_domain_type), intent(inout)  :: FUSS_domain
    type(ARES_domain_type), intent(inout)  :: ARES_domain
    character(len=llen_), intent(in)       :: file
    integer, intent(in)                    :: level
    ! Local
    character(len=llen_) :: path
    character(len=llen_) :: localpath_vtk
    integer              :: E_IO
    character(len=clen_) :: format(2)


    ! bc%Tg / bc%Pg still hold the stencil of the last BC_Coupled_Fluxes call, which ran
    ! BEFORE the last solver update: refresh it, so that the wall file is built from the
    ! current state of both solvers and the two sides of an interface report the same Tw.
    ! Both exchanges are no-ops in serial and on a single rank; every rank reaches this
    ! point, so the point-to-point pairing stays matched.
    call exchange_coupling_Tg_P    ( FUSS_domain, ARES_domain, level )
    call exchange_coupling_Pg_matID( FUSS_domain, ARES_domain, level )
    call Fill_Coupled_Ghost_Cell   ( FUSS_domain, ARES_domain, level )

    ! Update IOwall variables
    call Extract_Properties_From_BC( FUSS_domain, ARES_domain, level )


    ! Only rank 0 writes the consolidated wall output file: Extract_Properties_From_BC
    ! has already reduced every rank's contribution onto it.
    if (mpi_is_root) then
      path = 'OUTPUT/'
      call parse(hydra_obj_io%coupling_format,' ', format)
      
      ! Write the IOwall accordingly to the solution format
      E_IO = 0
      select case(trim(format(1)))
        case('vtk')
        IOwall(level)%vtk%format = trim(format(2))
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        !TODO: manca time per time-accurate simulation
        E_IO = vtk_write_structured_multiblock(orion=IOwall(level),vtspath=trim(localpath_vtk)//trim(file), &
                                                                  vtmpath=trim(path)//trim(file),varnames=hydra_obj_io%Owallnames)
      case('tecplot')
        IOwall(level)%tec%format = trim(format(2))
        E_IO = tec_write_structured_multiblock(orion=IOwall(level),varnames=hydra_obj_io%Owallnames,filename=trim(path)//trim(file)//trim(IOwall(level)%tec%extension))
      end select
      ! A failed write (disk full, quota, permission) leaves a cut file: say so, and go on (the next output may
      ! succeed), as the writers of ARES
      if ( E_IO /= 0 ) write(*,'(A,I0,A)') ' [WARNING] '//trim(path)//trim(file)//' was not written (status ', E_IO, &
                                           '): disk full, quota or permission?'
    endif

  end subroutine Write_Wall_Solution



  subroutine Extract_Properties_From_BC ( FUSS_domain, ARES_domain, level )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use hydraAF_Lib_BC_Coupled_Wall_Temperature, only: FUSS_BC_Coupled_Wall_Temperature, ARES_BC_Coupled_Wall_Temperature
    use hydraAF_Lib_Coupling_Segments, only: Segment_Wall_Output_FUSS, Segment_Wall_Output_ARES, &
                                             Update_Coupling_Segment_Partners
    use hydraAF_Mod_CouplingExchange, only: is_FUSS_bc_local_local, is_ARES_bc_local_local, &
                                            get_F_remote_P, get_M_remote_matID, &
                                            get_k_wall_FUSS, get_k_wall_ARES, &
                                            get_partner_half_FUSS, get_partner_half_ARES
#ifdef USE_MPI                                      
    use FUSS_Mod_MPI, only: FUSS_block_owner => block_owner, mpi_rank_, mpi_is_root, mpi_reduce_sum_r8_array
    use ARES_Mod_MPI, only: ARES_block_owner => block_owner
#endif   
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(ARES_domain_type), intent(inout) :: ARES_domain
    integer, intent(in)                   :: level
    ! Local
    real(R8), dimension(8) :: Ovar
    integer  :: f, lower, upper, i, b, bb, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    integer :: n_wall_total, ARES_n_wall_total, wall_i
    real(R8), allocatable :: FUSS_wall_buf(:), ARES_wall_buf(:)
    real(R8) :: FUSS_new_flux, ARES_new_flux(nprim)
    real(R8) :: P_con(np), matID_con, A_con, Teff_con
    real(R8), pointer :: p_remote(:)

    ! 104: current partner cells for the segments, also across ranks (collective)
    call Update_Coupling_Segment_Partners ( FUSS_domain, ARES_domain, level )

    ! FUSS
    n_wall_total = 0
    do b = 1, FUSS_domain%nb
      do f = 1, 6
        if (FUSS_obj_io_bc%coupling_flag(b,f)) n_wall_total = n_wall_total + FUSS_domain%n_bf(b,f)
      end do
    end do

    allocate(FUSS_wall_buf(2 * n_wall_total))
    FUSS_wall_buf = 0d0
   
    ! Each rank computes Ovar only for its local blocks
    upper = 0; wall_i = 0
    do b = 1, FUSS_domain % nb
      do f = 1, 6
        lower = upper + 1                        ! Update lower bound
        upper = upper + FUSS_domain % n_bf(b,f)  ! Upper bound: add number of cells on face f of block b

        do i = lower, upper
          Ovar = 0d0
          select case ( FUSS_domain % bc(i) % type )
            case (103)
              Bm  = FUSS_domain % bc(i) % b
              wall_i = wall_i + 1
#ifdef USE_MPI
              ! Skip entries whose FUSS primary block is not on this rank.
              if (allocated(FUSS_block_owner)) then
                if (FUSS_block_owner(Bm) /= mpi_rank_) cycle
              end if
#endif
              Im  = FUSS_domain % bc(i) % i
              Jm  = FUSS_domain % bc(i) % j
              Km  = FUSS_domain % bc(i) % k
              Fm  = FUSS_domain % bc(i) % f
              Bs  = FUSS_domain % bc(i) % bs
              Is  = FUSS_domain % bc(i) % is
              Js  = FUSS_domain % bc(i) % js
              Ks  = FUSS_domain % bc(i) % ks
              Fs  = FUSS_domain % bc(i) % fs
             
              ! Get P(1:np) of ARES secondary cell: local or remote cache.
#ifdef USE_MPI
              if (allocated(ARES_block_owner) .and. ARES_block_owner(Bs) /= mpi_rank_) then
                p_remote => get_F_remote_P(level, i)
                if (associated(p_remote)) then
                  P_con(1:np) = p_remote(1:np)
                else
                  P_con(1:np) = ARES_domain % blk(Bs) % P(1:np, Is, Js, Ks)
                end if
              else
                P_con(1:np) = ARES_domain % blk(Bs) % P(1:np, Is, Js, Ks)
              end if
#else
              P_con(1:np) = ARES_domain % blk(Bs) % P(1:np, Is, Js, Ks)
#endif
              ! Same mirrored wall conductivity and same partner half of the balance
              ! the solver used, so the wall file reports the temperature and the flux
              ! that were actually applied.
              call get_partner_half_FUSS(level, i, A_con, Teff_con)
              call FUSS_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                          P_con, FUSS_domain%bc(i)%Tg(1), A_con, Teff_con, &
                                          FUSS_new_flux, Ovar(6:7), &
                                          k_fluid_wall = get_k_wall_FUSS(level, i) )
           
              FUSS_wall_buf(2*(wall_i-1)+1 : 2*wall_i) = Ovar(6:7)

            case (104)
              Bm  = FUSS_domain % bc(i) % b
              wall_i = wall_i + 1
#ifdef USE_MPI
              if (allocated(FUSS_block_owner)) then
                if (FUSS_block_owner(Bm) /= mpi_rank_) cycle
              end if
#endif

              ! 104: the segments of the face, averaged by area (hydraAF_Lib_Coupling_Segments)
              call Segment_Wall_Output_FUSS ( FUSS_domain, level, i, Ovar(6:7) )

              FUSS_wall_buf(2*(wall_i-1)+1 : 2*wall_i) = Ovar(6:7)
          end select
        enddo

      enddo ! face
    enddo ! block

    ! Reduce FUSS_wall_buf to root. Each cell is owned by exactly one rank, so SUM gives
    ! correct values.
    call mpi_reduce_sum_r8_array(FUSS_wall_buf, 2*n_wall_total)

    ! Root fills IOwall from the gathered buffer
    if (mpi_is_root) then
      upper = 0; wall_i = 0; bb = 0
      do b = 1, FUSS_domain%nb
        do f = 1, 6
          lower = upper + 1
          upper = upper + FUSS_domain%n_bf(b,f)
          if (.not. FUSS_obj_io_bc%coupling_flag(b,f)) cycle
          bb = bb + 1

          do i = lower, upper
            ! Count exactly the entries the first loop counted: a face flagged as coupled
            ! may still carry non-103/104 entries, and taking them in would shift the buffer
            ! for the whole rest of the block.
            if (FUSS_domain%bc(i)%type /= 103 .and. FUSS_domain%bc(i)%type /= 104) cycle
            wall_i = wall_i + 1
            Ovar(6:7) = FUSS_wall_buf(2*(wall_i-1)+1 : 2*wall_i)

            select case(f)
              case(1,2); im = 1                    ; jm = FUSS_domain%bc(i)%j ; km = FUSS_domain%bc(i)%k
              case(3,4); im = FUSS_domain%bc(i)% i ; jm = 1                   ; km = FUSS_domain%bc(i)%k
              case(5,6); im = FUSS_domain%bc(i)% i ; jm = FUSS_domain%bc(i)%j ; km = 1
            end select
            IOwall(level)%block(bb)%vars(1:5, im, jm, km) = 1786.25d-10
            IOwall(level)%block(bb)%vars(6:7, im, jm, km) = Ovar(6:7)

        enddo

        enddo ! face
      enddo ! block
    endif
    
    deallocate(FUSS_wall_buf)
    
    ! ARES
    ARES_n_wall_total = 0
    do b = 1, ARES_domain%nb
      do f = 1, 6
        if (ARES_obj_io_bc%coupling_flag(b,f)) ARES_n_wall_total = ARES_n_wall_total + ARES_domain%n_bf(b,f)
      end do
    end do

    allocate(ARES_wall_buf(7 * ARES_n_wall_total))
    ARES_wall_buf = 0d0
    ! ARES
    upper = 0; wall_i = 0
    do b = 1, ARES_domain % nb
      do f = 1, 6
        lower = upper + 1                        ! Update lower bound
        upper = upper + ARES_domain % n_bf(b,f)  ! Upper bound: add number of cells on face f of block b

        do i = lower, upper
          Ovar = 0d0
          select case ( ARES_domain % bc(i) % type )
            case (103)
              Bm  = ARES_domain % bc(i) % b
              wall_i = wall_i + 1
#ifdef USE_MPI
              ! Skip entries whose ARES primary block is not on this rank.
              if (allocated(ARES_block_owner)) then
                if (ARES_block_owner(Bm) /= mpi_rank_) cycle
              end if
#endif
              Im  = ARES_domain % bc(i) % i
              Jm  = ARES_domain % bc(i) % j
              Km  = ARES_domain % bc(i) % k
              Fm  = ARES_domain % bc(i) % f
              Bs  = ARES_domain % bc(i) % bs
              Is  = ARES_domain % bc(i) % is
              Js  = ARES_domain % bc(i) % js
              Ks  = ARES_domain % bc(i) % ks
              Fs  = ARES_domain % bc(i) % fs
              ! Get matID of FUSS secondary cell: local or remote cache.
#ifdef USE_MPI
              if (allocated(FUSS_block_owner) .and. FUSS_block_owner(Bs) /= mpi_rank_) then
                matID_con = get_M_remote_matID(level, i)
              else
                matID_con = FUSS_domain % blk(Bs) % matID(Is, Js, Ks)
              end if
#else
              matID_con = FUSS_domain % blk(Bs) % matID(Is, Js, Ks)
#endif
              ! k_fluid_wall only: the store must not be advanced from the I/O path.
              call get_partner_half_ARES(level, i, A_con, Teff_con)
              call ARES_BC_Coupled_Wall_Temperature( Im, Jm, Km, Fm, ARES_domain % blk(Bm), &
                                          matID_con, ARES_domain%bc(i)%Pg(1,1), A_con, Teff_con, &
                                          ARES_new_flux, Ovar, &
                                          k_fluid_wall = get_k_wall_ARES(level, i) )
             
              ARES_wall_buf(7*(wall_i-1)+1 : 7*wall_i) = Ovar(1:7)

            case (104)
              Bm  = ARES_domain % bc(i) % b
              wall_i = wall_i + 1
#ifdef USE_MPI
              if (allocated(ARES_block_owner)) then
                if (ARES_block_owner(Bm) /= mpi_rank_) cycle
              end if
#endif

              ! 104: the segments of the face, averaged by area (hydraAF_Lib_Coupling_Segments)
              call Segment_Wall_Output_ARES ( ARES_domain, level, i, Ovar )

              ARES_wall_buf(7*(wall_i-1)+1 : 7*wall_i) = Ovar(1:7)

          end select
        enddo

      enddo ! face
    enddo ! block

    ! Reduce ARES_wall_buf to root. Each cell is owned by exactly one rank, so SUM gives
    ! correct values.
    call mpi_reduce_sum_r8_array(ARES_wall_buf, 7*ARES_n_wall_total)

    ! Root fills IOwall from the gathered buffer
    if (mpi_is_root) then
      upper = 0; wall_i = 0
      do b = 1, ARES_domain%nb
        do f = 1, 6
          lower = upper + 1
          upper = upper + ARES_domain%n_bf(b,f)
          if (.not.  ARES_obj_io_bc%coupling_flag(b,f)) cycle
          bb = bb + 1

          do i = lower, upper
            if (ARES_domain%bc(i)%type /= 103 .and. ARES_domain%bc(i)%type /= 104) cycle
            wall_i = wall_i + 1
            Ovar(1:7) = ARES_wall_buf(7*(wall_i-1)+1 : 7*wall_i)

            select case(f)
              case(1,2); im = 1                    ; jm = ARES_domain%bc(i)%j ; km = ARES_domain%bc(i)%k
              case(3,4); im = ARES_domain%bc(i)% i ; jm = 1                   ; km = ARES_domain%bc(i)%k
              case(5,6); im = ARES_domain%bc(i)% i ; jm = ARES_domain%bc(i)%j ; km = 1
            end select
            IOwall(level)%block(bb)%vars(1,im, jm, km)    = Ovar(1)
            IOwall(level)%block(bb)%vars(2:5, im, jm, km) = Ovar(2:5)
            IOwall(level)%block(bb)%vars(6:7, im, jm, km) = Ovar(6:7)
        

          enddo

        enddo ! face
      enddo ! block
    endif

    deallocate(ARES_wall_buf)
    
  end subroutine Extract_Properties_From_BC
  
end module hydraAF_IO_Wall