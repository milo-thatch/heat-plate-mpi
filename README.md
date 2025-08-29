# Solve the 2D heat equation with MPI!
This Fortran program solve the 2D homogeneous heat equation 

$$ \left(\partial_t - \partial_{xx} - \partial_{yy}\right) T = 0 $$

with Dirichlet boundary conditions on a cartesian grid. It uses a cartesian parallelization to speed up the computation.

## For the impatient
Follow these steps
1. Install MPI if needed by running `sudo apt install libopenmpi-dev openmpi-bin -y` in the Bash terminal in Ubuntu/Linux.
2. Set up the parameter file `parameters.f90`. Specify the number of iterations, the domain size and the number of grid points. 
3. Open the bash script `deploy.sh` and specify the number of processors `X` in `mpirun -np X ./heated_plate_mpi.x`
4. Change the permissions of the bash script by running `chmod +x deploy.sh` in the Bash terminal in Ubuntu/Linux.
5. Run `./deploy.sh` in the Bash terminal in Ubuntu/Linux to start the computation.

## Technicalities 
The code uses an Euler explicit scheme for time integration and centred finite differences for the spatial derivatives. The integration time step is determined on the grid size to ensure numerical stability

$$ \Delta t = \frac{1}{2} \left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right)^{-1} $$

Works for squared domains and an equal number of grid points along $x$ and $y$. 
