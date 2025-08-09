program heated_plate_mpi

    use parameters
    use heated_plate_mpi_subroutines
    implicit none

    !----------------------------------------------------------------------------------------------------
    ! MPI startup (general) -----------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
    call MPI_Init(ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_world, ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, size_world, ierror)

    !----------------------------------------------------------------------------------------------------
    ! Setup 2D Cartesian topology -----------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
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
    
    ! SOLVE THE HEAT EQUATION 
    ! second-order centred finite difference in x and y, euler explicit in time
    iteration = 0
    error = 1d3

    do while (iteration<50000)
        iteration = iteration + 1
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
        ! local max error
        local_error = maxval(abs(T - T_old))
        ! global max error
        call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm2d, ierror)
        error = global_error
        if(modulo(iteration,1000).eq.0 .and. rank_2d.eq.0) print*,'iteration = ',iteration
    end do

    ! PRINT THE RESULTS, FINALIZE
    call print_results
    if(rank_2d.eq.0) print*,'iterations = ',iteration

    deallocate(T,T_old)
    call MPI_Finalize(ierror)

end program heated_plate_mpi
