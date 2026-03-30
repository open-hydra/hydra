!>@brief: .
module hydraMF_IO_Wall
  use Lib_ORION_data
  use MOSE_IO_Wall, only: MOSE_coupling_flag => coupling_flag
  use FUSS_IO_Wall, only: FUSS_coupling_flag => coupling_flag
  use MOSE_Global_m
  use FUSS_Global_m

  implicit none
  integer                       :: Onwall
  type(ORION_data), allocatable :: IOwall(:)
  character(len=128)            :: Owallnames
  character(len=128)            :: wall_variables


contains


  subroutine Initialize_Wall_File ( FUSS_domain, MOSE_domain, level )
    use MOSE_Types_m
    use FUSS_Types_m
    use IR_Precision

    implicit none
    type(FUSS_domain_type), intent(in) :: FUSS_domain
    type(MOSE_domain_type), intent(in) :: MOSE_domain
    integer, intent(in)                :: level
    ! Local
    integer :: b, f, bb, i, j, k


    Onwall = 7
    Owallnames = '"y+" "tauX" "tauY" "tauZ" "pw" "Tw" "qw"'

    allocate(IOwall(level)%block(1:count(FUSS_coupling_flag)+count(MOSE_coupling_flag)))

    bb = 0
    
    do b = 1, FUSS_domain % nb; do f = 1, 6
      if (FUSS_coupling_flag( b , f )) then
        bb = bb + 1
        IOwall(level)%block(bb)%name = trim(FUSS_phase_prefix)//'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:0,0:FUSS_domain%blk(b)%dim(2),0:FUSS_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1,1:FUSS_domain%blk(b)%dim(2),1:FUSS_domain%blk(b)%dim(3)))
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
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1,1:FUSS_domain%blk(b)%dim(2),1:FUSS_domain%blk(b)%dim(3)))
          do k = 0, FUSS_domain%blk(b)%dim(3); do j = 0, FUSS_domain%blk(b)%dim(2)
              IOwall(level) % block(bb) % mesh(:,0,j,k) = FUSS_domain % blk(b) % node(FUSS_domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:0,0:FUSS_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:FUSS_domain%blk(b)%dim(1),1,1:FUSS_domain%blk(b)%dim(3)))
          do k = 0, FUSS_domain%blk(b)%dim(3); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = FUSS_domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:0,0:FUSS_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:FUSS_domain%blk(b)%dim(1),1,1:FUSS_domain%blk(b)%dim(3)))
          do k = 0, FUSS_domain%blk(b)%dim(3); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = FUSS_domain % blk(b) % node(i,FUSS_domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:FUSS_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:FUSS_domain%blk(b)%dim(1),1:FUSS_domain%blk(b)%dim(2),1))
          do j = 0, FUSS_domain%blk(b)%dim(2); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = FUSS_domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:FUSS_domain%blk(b)%dim(1),0:FUSS_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:FUSS_domain%blk(b)%dim(1),1:FUSS_domain%blk(b)%dim(2),1))
          do j = 0, FUSS_domain%blk(b)%dim(2); do i = 0, FUSS_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = FUSS_domain % blk(b) % node(i,j,FUSS_domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
      endif
    enddo; enddo

    do b = 1, MOSE_domain % nb; do f = 1, 6
      if (MOSE_coupling_flag( b , f )) then
        bb = bb + 1
        IOwall(level)%block(bb)%name = trim(MOSE_phase_prefix)//'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:0,0:MOSE_domain%blk(b)%dim(2),0:MOSE_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1,1:MOSE_domain%blk(b)%dim(2),1:MOSE_domain%blk(b)%dim(3)))
          do k = 0, MOSE_domain%blk(b)%dim(3); do j = 0, MOSE_domain%blk(b)%dim(2)
              IOwall(level) % block(bb) % mesh(1,0,j,k) = MOSE_domain % blk(b) % node(0,j,k) % c(1)
              IOwall(level) % block(bb) % mesh(2,0,j,k) = MOSE_domain % blk(b) % node(0,j,k) % c(2)
              IOwall(level) % block(bb) % mesh(3,0,j,k) = MOSE_domain % blk(b) % node(0,j,k) % c(3)
              ! Section to fool compiler...
              if (IOwall(level)%block(bb)%mesh(1,0,j,k)==0d0 .and. &
                  IOwall(level)%block(bb)%mesh(2,0,j,k)==0d0 .and. &
                  IOwall(level)%block(bb)%mesh(3,0,j,k)==0d0) &
              print*, IOwall(level) % block(bb) % mesh(:,0,j,k)
          enddo; enddo
        case(2)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:0,0:MOSE_domain%blk(b)%dim(2),0:MOSE_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1,1:MOSE_domain%blk(b)%dim(2),1:MOSE_domain%blk(b)%dim(3)))
          do k = 0, MOSE_domain%blk(b)%dim(3); do j = 0, MOSE_domain%blk(b)%dim(2)
              IOwall(level) % block(bb) % mesh(:,0,j,k) = MOSE_domain % blk(b) % node(MOSE_domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:MOSE_domain%blk(b)%dim(1),0:0,0:MOSE_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:MOSE_domain%blk(b)%dim(1),1,1:MOSE_domain%blk(b)%dim(3)))
          do k = 0, MOSE_domain%blk(b)%dim(3); do i = 0, MOSE_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = MOSE_domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:MOSE_domain%blk(b)%dim(1),0:0,0:MOSE_domain%blk(b)%dim(3)))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:MOSE_domain%blk(b)%dim(1),1,1:MOSE_domain%blk(b)%dim(3)))
          do k = 0, MOSE_domain%blk(b)%dim(3); do i = 0, MOSE_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,0,k) = MOSE_domain % blk(b) % node(i,MOSE_domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:MOSE_domain%blk(b)%dim(1),0:MOSE_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:MOSE_domain%blk(b)%dim(1),1:MOSE_domain%blk(b)%dim(2),1))
          do j = 0, MOSE_domain%blk(b)%dim(2); do i = 0, MOSE_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = MOSE_domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall(level)%block(bb)%mesh(1:3,0:MOSE_domain%blk(b)%dim(1),0:MOSE_domain%blk(b)%dim(2),0:0))
          allocate(IOwall(level)%block(bb)%vars(1:Onwall,1:MOSE_domain%blk(b)%dim(1),1:MOSE_domain%blk(b)%dim(2),1))
          do j = 0, MOSE_domain%blk(b)%dim(2); do i = 0, MOSE_domain%blk(b)%dim(1)
              IOwall(level) % block(bb) % mesh(:,i,j,0) = MOSE_domain % blk(b) % node(i,j,MOSE_domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
      endif
    enddo; enddo

  end subroutine Initialize_Wall_File



  subroutine Write_Wall_Solution( FUSS_domain, MOSE_domain, level, file, format )
    use Lib_VTK
    use Lib_Tecplot
    use MOSE_Types_m
    use FUSS_Types_m

    implicit none
    type(FUSS_domain_type), intent(inout)  :: FUSS_domain
    type(MOSE_domain_type), intent(inout)  :: MOSE_domain
    integer, intent(in)                    :: level
    character(len=300), intent(in)         :: file
    character(32), intent(in), optional    :: format(2)
    ! Local
    character(len=128)           :: path
    character(len=128)           :: localpath_vtk
    integer                      :: E_IO


    path = 'OUTPUT/'

    ! Update IOwall variables
    call Extract_Properties_From_BC( FUSS_domain, MOSE_domain, level )

    ! Write the IOwall accordingly to the solution format
    select case(format(1))
    case('vtk')
      IOwall(level)%vtk%format = format(2)
      localpath_vtk = trim(path)//'vtk/'
      call execute_command_line('mkdir -p '//trim(localpath_vtk))
      !TODO: manca time per time-accurate simulation
      E_IO = vtk_write_structured_multiblock(orion=IOwall(level),vtspath=trim(localpath_vtk)//trim(file), &
                                                                 vtmpath=trim(path)//trim(file),varnames=Owallnames)
    case('tecplot')
      IOwall(level)%tec%format = format(2)
      !TODO: manca time per time-accurate simulation
      E_IO = tec_write_structured_multiblock(orion=IOwall(level),varnames=Owallnames,filename=trim(path)//trim(file)//'.tec')
    end select

  end subroutine Write_Wall_Solution



  subroutine Extract_Properties_From_BC ( FUSS_domain, MOSE_domain, level )
    use MOSE_Types_m
    use MOSE_Global_m
    use FUSS_Types_m
    use FUSS_Global_m
    use hydraMF_Mod_Coupling

    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    integer, intent(in)                   :: level
    ! Local
    real(kind=8), dimension(8) :: Ovar
    integer  :: f, lower, upper, i, b, bb, im, jm, km, bm, fm, bs, fs, is, js, ks
    logical  :: assign
    real(R8) :: FUSS_new_flux, MOSE_new_flux(nprim)


    ! BC fluxes are computed in order of block and face type for both solvers.
    bb = 0

    ! FUSS
    upper = 0
    do b = 1, FUSS_domain % nb
      do f = 1, 6
        lower = upper + 1                        ! Update lower bound
        upper = upper + FUSS_domain % n_bf(b,f)  ! Upper bound: add number of cells on face f of block b

        assign = .false.
        if (FUSS_coupling_flag( b , f )) then
          bb = bb + 1
          assign = .true.
        else
          cycle
        endif

        do i = lower, upper
          Ovar = 0d0
          select case ( FUSS_domain % bc(i) % type )
            case (1000)
              Bm  = FUSS_domain % bc(i) % b
              Im  = FUSS_domain % bc(i) % i
              Jm  = FUSS_domain % bc(i) % j 
              Km  = FUSS_domain % bc(i) % k 
              Fm  = FUSS_domain % bc(i) % f
              Bs  = FUSS_domain % bc(i) % bs
              Is  = FUSS_domain % bc(i) % is
              Js  = FUSS_domain % bc(i) % js
              Ks  = FUSS_domain % bc(i) % ks
              Fs  = FUSS_domain % bc(i) % fs
              call FUSS_Wall_Temperature( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                          Is, Js, Ks, Fs, MOSE_domain % blk(Bs), &
                                          FUSS_domain%bc(i)%Mg(1), FUSS_domain%bc(i)%Tg, &
                                          FUSS_new_flux, Ovar(6:7) )
          end select

          if (assign) then
            select case(f)
              case(1,2); im = 1                    ; jm = FUSS_domain%bc(i)%j ; km = FUSS_domain%bc(i)%k
              case(3,4); im = FUSS_domain%bc(i)% i ; jm = 1                   ; km = FUSS_domain%bc(i)%k
              case(5,6); im = FUSS_domain%bc(i)% i ; jm = FUSS_domain%bc(i)%j ; km = 1
            end select
            IOwall(level)%block(bb)%vars(1:5, im, jm, km) = 1786.25d-10
            IOwall(level)%block(bb)%vars(6:7, im, jm, km) = Ovar(6:7)
          endif

        enddo

      enddo ! face
    enddo ! block

    ! MOSE
    upper = 0
    do b = 1, MOSE_domain % nb
      do f = 1, 6
        lower = upper + 1                        ! Update lower bound
        upper = upper + MOSE_domain % n_bf(b,f)  ! Upper bound: add number of cells on face f of block b

        assign = .false.
        if (MOSE_coupling_flag( b , f )) then
          bb = bb + 1
          assign = .true.
        else
          cycle
        endif

        do i = lower, upper
          Ovar = 0d0
          select case ( MOSE_domain % bc(i) % type )
            case (1000)
              Bm  = MOSE_domain % bc(i) % b
              Im  = MOSE_domain % bc(i) % i
              Jm  = MOSE_domain % bc(i) % j 
              Km  = MOSE_domain % bc(i) % k 
              Fm  = MOSE_domain % bc(i) % f
              Bs  = MOSE_domain % bc(i) % bs
              Is  = MOSE_domain % bc(i) % is
              Js  = MOSE_domain % bc(i) % js
              Ks  = MOSE_domain % bc(i) % ks
              Fs  = MOSE_domain % bc(i) % fs
              call MOSE_Wall_Temperature( Im, Jm, Km, Fm, MOSE_domain % blk(Bm), &
                                          Is, Js, Ks, Fs, FUSS_domain % blk(Bs), &
                                          MOSE_domain%bc(i)%Mg(1), MOSE_domain%bc(i)%Pg, &
                                          MOSE_new_flux, Ovar )
          end select

          if (assign) then
            select case(f)
              case(1,2); im = 1                    ; jm = MOSE_domain%bc(i)%j ; km = MOSE_domain%bc(i)%k
              case(3,4); im = MOSE_domain%bc(i)% i ; jm = 1                   ; km = MOSE_domain%bc(i)%k
              case(5,6); im = MOSE_domain%bc(i)% i ; jm = MOSE_domain%bc(i)%j ; km = 1
            end select
            IOwall(level)%block(bb)%vars(1,im, jm, km) = Ovar(1)
            IOwall(level)%block(bb)%vars(2:5, im, jm, km) = Ovar(2:5)
            IOwall(level)%block(bb)%vars(6:7, im, jm, km) = Ovar(6:7)
          endif

        enddo

      enddo ! face
    enddo ! block


  end subroutine Extract_Properties_From_BC
  

end module hydraMF_IO_Wall
