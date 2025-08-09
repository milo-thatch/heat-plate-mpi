module common_variables

    use mpi 

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

end module common_variables