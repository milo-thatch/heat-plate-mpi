program heated_plate_mpi

    use parameters
    use heated_plate_mpi_subroutines
    implicit none

    real(kind=8) :: start_wall, end_wall, elapsed_wall
    real(kind=8) :: start_cpu, end_cpu, elapsed_cpu

    !----------------------------------------------------------------------------------------------------
    ! MPI startup (general) -----------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
    call MPI_Init(ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_world, ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, size_world, ierror)
    ! nobody is perfect
    if( nx.ne.ny ) then
        if(rank_world.eq.0) print*, 'nx and ny are different, this is above my pay grade'
        call MPI_Finalize(ierror)
        stop
    end if

    ! start timers
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) ! sync before timing
    call cpu_time(start_cpu) ! CPU time
    start_wall = MPI_WTIME() ! wall-clock time

    !----------------------------------------------------------------------------------------------------
    ! Setup 2D Cartesian topology -----------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
    if( abs(sqrt(dble(size_world)) - int(sqrt(dble(size_world)))) .gt. 1e-10 ) then
        if(rank_world.eq.0) print*, 'the number of processors ',size_world,' is not the square of an integer'
        if(rank_world.eq.0) print*, 'unable to generate an MPI square cartesian topology'
        call MPI_Finalize(ierror)
        stop
    else if( modulo(dble(nx),sqrt(dble(size_world))).ne.0.0d0 ) then 
        if(rank_world.eq.0) print*, 'either nx = ',nx,' or ',ny,&
                            ' are NOT the multiple of of the number of processors ',size_world
        call MPI_Finalize(ierror)
        stop
    end if
    dims = (/ int(sqrt(dble(size_world))), int(sqrt(dble(size_world))) /)
    call MPI_Dims_create(size_world, 2, dims, ierror) ! creates 2D topology
    periods = (/ .false., .false. /)  ! No periodic BCs
    ! creates a 2d communicator --> arranges the processors on a 2d grid
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .false., comm2d, ierror) 
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
    dx = Lx/nx
    dy = Ly/ny
    dt = 0.5d0/(1.0d0/dx**2.0d0 + 1.0d0/dy**2.0d0 )

    ! SET BOUNDARY CONDITIONS
    call set_boundary_conditions
    
    !----------------------------------------------------------------------------------------------------
    ! SOLVE THE HEAT EQUATION ---------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
    ! second-order centred finite difference in x and y, euler explicit in time
    iter = 0

    do while (iter<iterations)
        iter = iter + 1
        T_old = T
        ! update interior points
        do j = 1, ny_local
            do i = 1, nx_local
                T(i,j) = T_old(i,j) + dt * ( &
               (T_old(i+1,j) - 2.0d0*T_old(i,j) + T_old(i-1,j)) / dx**2 + &
                (T_old(i,j+1) - 2.0d0*T_old(i,j) + T_old(i,j-1)) / dy**2 )
            end do
        end do
        ! update the halos
        call update_halos
        if(modulo(iter,1000).eq.0 .and. rank_2d.eq.0) print*,'iter = ',iter
    end do

    !----------------------------------------------------------------------------------------------------
    ! PRINT THE RESULTS, FINALIZE -----------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
    call print_results
    deallocate(T,T_old)
    
    ! end timers
    end_wall = MPI_WTIME()
    call cpu_time(end_cpu)
    elapsed_cpu  = end_cpu  - start_cpu
    elapsed_wall = end_wall - start_wall
    print *, 'Rank', rank_world, ': CPU time =', elapsed_cpu, 's ; Wall time =', elapsed_wall, 's'

    call MPI_Finalize(ierror)

end program heated_plate_mpi
