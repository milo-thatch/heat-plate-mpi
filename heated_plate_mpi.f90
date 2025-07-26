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
    integer :: dims(2), coords(2), nx_local, ny_local
    logical :: periods(2)
    integer :: root ! identify the master (root processor)
    integer :: istart, iend, jstart, jend
    integer :: src ! mpi_send source for final gathering
    real(kind=8), allocatable :: Tglob(:,:)
    real(kind=8), allocatable :: buf(:,:)

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
    T = 0.0d0
    T_old = 0.0d0

    ! impose boundary conditions
    if(coords(1).eq.0) then ! dirichlet condition (zero temperature)
        T(0,:) = 1.0d0
        T_old(0,:) = 1.0d0
    end if 
    if(coords(2).eq.0) then ! dirichlet condition (zero temperature)
        T(:,0) = 1.0d0
        T_old(:,0) = 1.0d0
    end if 
    if(coords(1).eq.(dims(1)-1)) then ! dirichlet condition (linear temperature distribution)
        do j=0,ny_local+1
            T(nx_local+1,j)=dble(ny_local+1-j)/dble(ny_local+1)
        enddo
    end if
    if(coords(2).eq.(dims(2)-1)) then ! dirichlet condition (linear temperature distribution)
        do j=0,ny_local+1
            T(j,ny_local+1)=dble(nx_local+1-j)/dble(nx_local+1)
        enddo
    end if

    print*, 'Hi, I am processor ',rank_world,' of ',size_world

    root = 0

    ! Allocate receive buffer for root (global array incl. halos)
    if (rank_world == root) then
        allocate(Tglob(0:nx+1, 0:ny+1))
        Tglob = 0.0d0
    end if

    ! Pack local interior (no halos) into a contiguous buffer
    allocate(buf(1:nx_local, 1:ny_local))
    buf = T(1:nx_local, 1:ny_local)

    ! Send/recv to root
    if (rank_world == root) then
        ! First, copy my own block
        istart = 1 + coords(1)*nx_local
        iend   = istart + nx_local - 1
        jstart = 1 + coords(2)*ny_local
        jend   = jstart + ny_local - 1
        Tglob(istart:iend, jstart:jend) = buf

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
        end do
    else
        ! Send my coords, then my data
        call MPI_Send(coords, 2, MPI_INTEGER, root, 100, comm2d, ierror)
        call MPI_Send(buf, nx_local*ny_local, MPI_DOUBLE_PRECISION, root, 101, comm2d, ierror)
    end if

    deallocate(buf)

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