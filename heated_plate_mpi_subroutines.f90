module heated_plate_mpi_subroutines

    use parameters
    use mpi 

    ! MPI_COMM_WORLD variables
    integer :: ierror, rank_world, size_world
    integer :: status(MPI_STATUS_SIZE)

    ! MPI_CART variables
    integer :: rank_2d

    ! counters
    integer :: i,j,iter

    ! Mathematical formulation
    real(kind=8), allocatable, dimension(:,:) :: T, T_old
    real(kind=8) :: dx, dy, dt

    ! MPI comm2d variables
    integer :: comm2d ! cartesian communicator, opaque handle and therefore an integer in fortran
    integer :: north, south, west, east
    integer :: dims(2), coords(2), coords_src(2), nx_local, ny_local
    logical :: periods(2)
    integer :: root ! identify the master (root processor)
    integer :: istart, iend, jstart, jend
    integer :: src ! mpi_send source for final gathering

contains
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
subroutine set_boundary_conditions

    integer :: gi, gj

    ! Left boundary (x=0)
    if (coords(1) == 0) then
        do j = 0, ny_local+1
            T(0,j) = 0.0d0
            T_old(0,j) = T(0,j)
        end do
    end if

    ! Bottom boundary (y=0)
    if (coords(2) == 0) then
        do i = 0, nx_local+1
            T(i,0) = 0.0d0
            T_old(i,0) = T(i,0)
        end do
    end if

    ! Right boundary (x=Lx)
    if (coords(1) == dims(1)-1) then
        do j = 0, ny_local+1
            gj = coords(2)*ny/dims(2) + j
            T(nx_local+1,j) = dble(gj) / dble(ny)
            T_old(nx_local+1,j) = T(nx_local+1,j)
        end do
    end if

    ! Top boundary (y=Ly)
    if (coords(2) == dims(2)-1) then
        do i = 0, nx_local+1
            gi = coords(1)*nx/dims(1) + i
            T(i,ny_local+1) = dble(gi) / dble(nx)
            T_old(i,ny_local+1) = T(i,ny_local+1)
        end do
    end if

end subroutine set_boundary_conditions
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
subroutine print_results

    real(kind=8), allocatable, dimension(:,:) :: Tglob, buf
    real(kind=8), allocatable, dimension(:) :: buf_bound_x, buf_bound_y

    ! Root takes them all to print the dat file with all the points
    root = 0

    ! Allocate receive buffer for root (global array incl. halos)
    if (rank_world == root) then
        allocate(Tglob(0:nx+1, 0:ny+1))
        Tglob = 0.0d0
    end if

    ! Pack local interior (no halos) into a contiguous buffer
    allocate(buf(1:nx_local, 1:ny_local),buf_bound_x(0:nx_local),buf_bound_y(0:ny_local))
    buf = T(1:nx_local, 1:ny_local)

    ! Send/recv to root
    if (rank_world == root) then
        ! First, copy my own block
        istart = 1 + coords(1)*nx_local
        iend   = istart + nx_local - 1
        jstart = 1 + coords(2)*ny_local
        jend   = jstart + ny_local - 1
        Tglob(istart:iend, jstart:jend) = T(1:nx_local, 1:ny_local)
        ! am i in the corners?
        Tglob(istart-1:iend, 0) = T(0:nx_local, 0)
        Tglob(0, jstart-1:jend) = T(0, 0:ny_local)

        ! Now receive from all other ranks
        do src = 0, size_world-1
            if (src == root) cycle

            call MPI_Recv(coords, 2, MPI_INTEGER, src, 100, comm2d, status, ierror)

            istart = 1 + coords(1)*nx_local
            iend   = istart + nx_local - 1
            jstart = 1 + coords(2)*ny_local
            jend   = jstart + ny_local - 1

            call MPI_Recv(Tglob(istart:iend, jstart:jend), nx_local*ny_local, MPI_DOUBLE_PRECISION, &
                        src, 101, comm2d, status, ierror)
            ! check boundaries: corners
            if((coords(1).eq.0).and.(coords(2).eq.(dims(2)-1))) then ! NORTHWEST
                call MPI_Recv(Tglob(istart-1:iend,ny+1), nx_local+1,MPI_DOUBLE_PRECISION, &
                        src, 102, comm2d, status, ierror)
                call MPI_Recv(Tglob(0,jstart:jend+1), ny_local+1,MPI_DOUBLE_PRECISION, &
                        src, 103, comm2d, status, ierror)
            end if
            if((coords(1).eq.(dims(2)-1)).and.(coords(2).eq.(dims(2)-1))) then ! NORTHEAST
                call MPI_Recv(Tglob(istart:iend+1,ny+1), nx_local+1, MPI_DOUBLE_PRECISION, &
                        src, 104, comm2d, status, ierror)
                call MPI_Recv(Tglob(nx+1,jstart:jend+1), ny_local+1, MPI_DOUBLE_PRECISION, &
                        src, 105, comm2d, status, ierror)
            end if
            if((coords(1).eq.0).and.(coords(2).eq.0)) then ! SOUTHWEST
                call MPI_Recv(Tglob(istart-1:iend,0), nx_local+1,MPI_DOUBLE_PRECISION, &
                        src, 106, comm2d, status, ierror)
                call MPI_Recv(Tglob(0,jstart-1:jend), ny_local+1,MPI_DOUBLE_PRECISION, &
                        src, 107, comm2d, status, ierror)
            end if
            if((coords(1).eq.(dims(2)-1)).and.(coords(2).eq.0)) then ! SOUTHEAST
                call MPI_Recv(Tglob(istart:iend+1,0), nx_local+1, MPI_DOUBLE_PRECISION, &
                        src, 108, comm2d, status, ierror)
                call MPI_Recv(Tglob(nx+1,jstart-1:jend), ny_local+1, MPI_DOUBLE_PRECISION, &
                        src, 109, comm2d, status, ierror)
            end if
            ! check boundaries: edges (no corners)
            if((coords(1).ne.0).and.(coords(1).ne.(dims(2)-1)).and.(coords(2).eq.(dims(2)-1))) then ! NORTH
                call MPI_Recv(Tglob(istart:iend,ny+1), nx_local, MPI_DOUBLE_PRECISION, &
                        src, 110, comm2d, status, ierror)
            end if
            if((coords(1).ne.0).and.(coords(1).ne.(dims(2)-1)).and.(coords(2).eq.0)) then ! SOUTH
                call MPI_Recv(Tglob(istart:iend,ny+1), nx_local, MPI_DOUBLE_PRECISION, &
                        src, 111, comm2d, status, ierror)
            end if
            if((coords(2).ne.0).and.(coords(2).ne.(dims(2)-1)).and.(coords(1).eq.0)) then ! WEST
                buf_bound_y = T(0, 1:ny_local)
                call MPI_Recv(Tglob(0,jstart:jend), ny_local, MPI_DOUBLE_PRECISION, &
                        src, 112, comm2d, status, ierror)
            end if
            if((coords(2).ne.0).and.(coords(2).ne.(dims(2)-1)).and.(coords(1).eq.(dims(1)-1))) then ! EAST
                call MPI_Recv(Tglob(nx+1,jstart:jend), ny_local, MPI_DOUBLE_PRECISION, &
                        src, 113, comm2d, status, ierror)
            end if
        end do
    else
        ! Send my coords, then my data
        call MPI_Send(coords, 2, MPI_INTEGER, root, 100, comm2d, ierror)
        call MPI_Send(buf, nx_local*ny_local, MPI_DOUBLE_PRECISION, root, 101, comm2d, ierror)
        ! check boundaries
        if((coords(1).eq.0).and.(coords(2).eq.(dims(2)-1))) then ! NORTHWEST
            buf_bound_x = T(0:nx_local, ny_local+1)
            call MPI_Send(buf_bound_x, nx_local+1, MPI_DOUBLE_PRECISION, &
                    root, 102, comm2d, ierror)
            buf_bound_y = T(0,1:ny_local+1)
            call MPI_Send(buf_bound_y, ny_local+1, MPI_DOUBLE_PRECISION, &
                    root, 103, comm2d, ierror)
        end if
        if((coords(1).eq.(dims(1)-1)).and.(coords(2).eq.(dims(2)-1))) then ! NORTHEAST
            buf_bound_x = T(1:nx_local+1, ny_local+1)
            call MPI_Send(buf_bound_x, nx_local+1, MPI_DOUBLE_PRECISION, &
                    root, 104, comm2d, ierror)
            buf_bound_y = T(nx_local+1,1:ny_local+1)
            call MPI_Send(buf_bound_y, ny_local+1, MPI_DOUBLE_PRECISION, &
                    root, 105, comm2d, ierror)
        end if
        if((coords(1).eq.0).and.(coords(2).eq.0)) then ! SOUTHWEST
            buf_bound_x = T(0:nx_local, 0)
            call MPI_Send(buf_bound_x, nx_local+1, MPI_DOUBLE_PRECISION, &
                    root, 106, comm2d, ierror)
            buf_bound_y = T(0,0:ny_local)
            call MPI_Send(buf_bound_y, ny_local+1, MPI_DOUBLE_PRECISION, &
                    root, 107, comm2d, ierror)
        end if
        if((coords(1).eq.(dims(1)-1)).and.(coords(2).eq.0)) then ! SOUTHEAST
            buf_bound_x = T(1:nx_local+1, 0)
            call MPI_Send(buf_bound_x, nx_local+1, MPI_DOUBLE_PRECISION, &
                    root, 108, comm2d, ierror)
            buf_bound_y = T(nx_local+1,0:ny_local)
            call MPI_Send(buf_bound_y, ny_local+1, MPI_DOUBLE_PRECISION, &
                    root, 109, comm2d, ierror)
        end if
        if((coords(1).ne.0).and.(coords(1).ne.(dims(2)-1)).and.(coords(2).eq.(dims(2)-1))) then ! NORTH
            buf_bound_x = T(1:nx_local, ny_local+1)
            call MPI_Send(buf_bound_x, nx_local, MPI_DOUBLE_PRECISION, &
                    root, 110, comm2d, ierror)
        end if
        if((coords(1).ne.0).and.(coords(1).ne.(dims(2)-1)).and.(coords(2).eq.0)) then ! SOUTH
            buf_bound_x = T(1:nx_local, 0)
            call MPI_Send(buf_bound_x, nx_local, MPI_DOUBLE_PRECISION, &
                    root, 111, comm2d, ierror)
        end if
        if((coords(2).ne.0).and.(coords(2).ne.(dims(2)-1)).and.(coords(1).eq.0)) then ! WEST
            buf_bound_y = T(0, 1:ny_local)
            call MPI_Send(buf_bound_y, ny_local, MPI_DOUBLE_PRECISION, &
                    root, 112, comm2d, ierror)
        end if
        if((coords(2).ne.0).and.(coords(2).ne.(dims(2)-1)).and.(coords(1).eq.(dims(1)-1))) then ! EAST
            buf_bound_y = T(nx_local+1, 1:ny_local)
            call MPI_Send(buf_bound_y, ny_local, MPI_DOUBLE_PRECISION, &
                    root, 113, comm2d, ierror)
        end if
    end if

    deallocate(buf,buf_bound_x,buf_bound_y)

    ! Root writes the full field (including global halos) to text
    if (rank_world == root) then
        open(49, file='T_global.dat', status='replace', action='write', form='formatted')
        do j = 0, ny+1
            do i = 0, nx+1
                write(49,'(ES23.15)', advance='no') Tglob(i,j)
                if (i < nx+1) write(49,'(A)', advance='no') ' '
            end do
            write(49,*)
        end do
        close(49)
        deallocate(Tglob)
    end if

end subroutine print_results
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
subroutine update_halos

    real(kind=8), allocatable, dimension(:) :: buf_inbound_x, buf_outbound_x
    real(kind=8), allocatable, dimension(:) :: buf_inbound_y, buf_outbound_y
    
    ! every rank sends its edge and gets an halo in return
    ! (unless they lie on boundaries of course)

    allocate(buf_inbound_x(ny_local), buf_outbound_x(ny_local))
    allocate(buf_inbound_y(nx_local), buf_outbound_y(nx_local))

    if(coords(1).ne.0) then ! sends to the west, receives from the west
        buf_outbound_x = T(1,1:ny_local) ! sends part of its domain x = 1
        call MPI_sendrecv(buf_outbound_x, ny_local, MPI_DOUBLE_PRECISION, west, 100, &
                        buf_inbound_x, ny_local, MPI_DOUBLE_PRECISION, west, 101, &
                        comm2d, status, ierror)
        T(0,1:ny_local) = buf_inbound_x ! gets a halo in return x = 0
    end if
    if(coords(1).ne.(dims(1)-1)) then ! sends to the east, receives from the east
        buf_outbound_x = T(nx_local,1:ny_local) ! sends part of its domain x = nx_local
        call MPI_sendrecv(buf_outbound_x, ny_local, MPI_DOUBLE_PRECISION, east, 101, &
                        buf_inbound_x, ny_local, MPI_DOUBLE_PRECISION, east, 100, &
                        comm2d, status, ierror)
        T(nx_local+1,1:ny_local) = buf_inbound_x ! gets a halo in return x = nx_local +1
    end if
    if(coords(2).ne.0) then ! sends to the south, receives from the south
        buf_outbound_y = T(1:nx_local,1) ! sends part of its domain y = 1
        call MPI_sendrecv(buf_outbound_y, nx_local, MPI_DOUBLE_PRECISION, south, 102, &
                        buf_inbound_y, nx_local, MPI_DOUBLE_PRECISION, south, 103, &
                        comm2d, status, ierror)
        T(1:nx_local,0) = buf_inbound_y ! gets a halo in return y = 0
    end if
    if(coords(2).ne.(dims(2)-1)) then ! sends to the north, receives from the north
        buf_outbound_y = T(1:nx_local,ny_local) ! sends part of its domain y = ny_local
        call MPI_sendrecv(buf_outbound_y, nx_local, MPI_DOUBLE_PRECISION, north, 103, &
                        buf_inbound_y, nx_local, MPI_DOUBLE_PRECISION, north, 102, &
                        comm2d, status, ierror)
        T(1:nx_local,ny_local+1) = buf_inbound_y ! gets a halo in return y = ny_local +1
    end if

    deallocate(buf_inbound_x, buf_outbound_x)
    deallocate(buf_inbound_y, buf_outbound_y)

end subroutine update_halos

end module heated_plate_mpi_subroutines