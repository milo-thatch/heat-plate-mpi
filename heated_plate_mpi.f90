program heated_plate_mpi

    use mpi 
    use parameters
    implicit none

    ! MPI_COMM_WORLD variables
    integer :: ierror, rank_world, size_world
    integer :: status(MPI_STATUS_SIZE)

    ! MPI_CART variables
    integer :: rank_2d

    ! counters
    integer :: i,j,iteration

    ! Mathematical formulation
    real(kind=8), allocatable, dimension(:,:) :: T, T_old
    real(kind=8) :: Lx=1.0d0, Ly=1.0d0
    real(kind=8) :: dx,dy,dt

    ! MPI comm2d variables
    integer :: comm2d ! cartesian communicator, opaque handle and therefore an integer in fortran
    integer :: north, south, west, east
    integer :: dims(2), coords(2), coords_src(2), nx_local, ny_local
    logical :: periods(2)
    integer :: root ! identify the master (root processor)
    integer :: istart, iend, jstart, jend
    integer :: src ! mpi_send source for final gathering
    real(kind=8), allocatable, dimension(:,:) :: Tglob, buf
    real(kind=8), allocatable, dimension(:) :: buf_bound_x, buf_bound_y

    !-----------------------------------------------------------------------------------------------------------------
    ! MPI startup (general) ------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------------------------
    call MPI_Init(ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_world, ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, size_world, ierror)

    !-----------------------------------------------------------------------------------------------------------------
    ! Setup 2D Cartesian topology ------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------------------------
    dims = (/ int(sqrt(dble(size_world))), int(sqrt(dble(size_world))) /)
    call MPI_Dims_create(size_world, 2, dims, ierror) ! creates 2D topology
    periods = (/ .false., .false. /)  ! No periodic BCs
    ! creates a 2d communicator --> arranges the processors on a 2d grid
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm2d, ierror) 
    ! Use the cartesian communicator from now on
    call MPI_Comm_rank(comm2d, rank_2d, ierror)
    call MPI_Cart_coords(comm2d, rank_2d, 2, coords, ierror)
    ! Get neighbors ranks (MPI_PROC_NULL if no neighbor)
    ! syntax is MPI_Cart_shift(communicator, direction, displacement, source, dest, ierror)
    call MPI_Cart_shift(comm2d, 0, 1, west, east, ierror) ! positive x west -> east
    call MPI_Cart_shift(comm2d, 1, 1, south, north, ierror) ! positive y south -> north

    ! each processor has access to its square + 2 halos
    ! global domain is (nx+2)*(ny+2)
    nx_local = nx/dims(1)
    ny_local = ny/dims(2)
    allocate(T(0:nx_local+1,0:ny_local+1),T_old(0:nx_local+1,0:ny_local+1)) ! allocate and include halos 
    ! initialize the temperature
    T = dble(rank_2d)
    T_old = dble(rank_2d)

    ! SET BOUNDARY CONDITIONS
    call set_boundary_conditions(T, T_old, nx_local, ny_local, coords, dims, nx, ny)  

    print*, 'Hi, I am processor ',rank_world,' of ',size_world

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
                call MPI_Recv(Tglob(istart-1:iend,ny+1), nx_local+1,MPI_DOUBLE_PRECISION, &
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
        end do
    else
        ! Send my coords, then my data
        call MPI_Send(coords, 2, MPI_INTEGER, root, 100, comm2d, ierror)
        call MPI_Send(buf, nx_local*ny_local, MPI_DOUBLE_PRECISION, root, 101, comm2d, ierror)
        ! check boundaries
        if((coords(1).eq.0).and.(coords(2).eq.(dims(2)-1))) then ! NORTHWEST
            buf_bound_x = T(0:nx_local, ny_local+1)
            call MPI_Send(buf_bound_x, nx_local+1,MPI_DOUBLE_PRECISION, &
                    root, 102, comm2d, ierror)
            buf_bound_y = T(0,1:ny_local+1)
            call MPI_Send(buf_bound_y, ny_local+1,MPI_DOUBLE_PRECISION, &
                    root, 103, comm2d, ierror)
        end if
        if((coords(1).eq.(dims(2)-1)).and.(coords(2).eq.(dims(2)-1))) then ! NORTHEAST
            buf_bound_x = T(1:nx_local+1, ny_local+1)
            call MPI_Send(buf_bound_x, nx_local+1,MPI_DOUBLE_PRECISION, &
                    root, 104, comm2d, ierror)
            buf_bound_y = T(nx_local+1,1:ny_local+1)
            call MPI_Send(buf_bound_y, ny_local+1,MPI_DOUBLE_PRECISION, &
                    root, 105, comm2d, ierror)
        end if
        if((coords(1).eq.0).and.(coords(2).eq.0)) then ! SOUTHWEST
            buf_bound_x = T(0:nx_local, 0)
            call MPI_Send(buf_bound_x, nx_local+1,MPI_DOUBLE_PRECISION, &
                    root, 106, comm2d, ierror)
            buf_bound_y = T(0,0:ny_local)
            call MPI_Send(buf_bound_y, ny_local+1,MPI_DOUBLE_PRECISION, &
                    root, 107, comm2d, ierror)
        end if
        if((coords(1).eq.(dims(2)-1)).and.(coords(2).eq.0)) then ! SOUTHEAST
            buf_bound_x = T(1:nx_local+1, 0)
            call MPI_Send(buf_bound_x, nx_local+1,MPI_DOUBLE_PRECISION, &
                    root, 108, comm2d, ierror)
            buf_bound_y = T(nx_local+1,0:ny_local)
            call MPI_Send(buf_bound_y, ny_local+1,MPI_DOUBLE_PRECISION, &
                    root, 109, comm2d, ierror)
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

    deallocate(T,T_old)
    call MPI_Finalize(ierror)

end program heated_plate_mpi
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
subroutine set_boundary_conditions(T, T_old, nx_local, ny_local, coords, dims, nx, ny)
    implicit none
    integer, intent(in) :: nx_local, ny_local, nx, ny
    integer, dimension(2), intent(in) :: coords, dims
    double precision, dimension(0:nx_local+1, 0:ny_local+1), intent(inout) :: T, T_old
    integer :: i, j, gi, gj

    ! Left boundary (x=0)
    if (coords(1) == 0) then
        do j = 0, ny_local+1
            T(0,j) = -1.0d0
            T_old(0,j) = -1.0d0
        end do
    end if

    ! Bottom boundary (y=0)
    if (coords(2) == 0) then
        do i = 0, nx_local+1
            T(i,0) = 1.0d0
            T_old(i,0) = 1.0d0
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