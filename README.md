# FemSolver
A Fortran code to solve an Initial-Boundary-value scalar wave problem in 2D. The solver uses semidiscrete formulation with p-FEM for space and implicit-Euler for time discretizations, respectively. The code is written for Linux, however, the system calls could be modified to be run over other operating systems as well. 
## Features
* Direct import 2D meshes generated in [Gmsh](https://gmsh.info/) - an open source mesh generator.
* Dynamic allocation of all variables and arrays, depending on the imported mesh and order of FEM polynomials.
* Plot the mesh and numerical solutions in [Paraview](https://www.paraview.org/) - an open source alternative to Tecplot.
* Apply Neuman, Dirichlet and Robin boundary conditions on edges marked respectively in Gmsh. 
* LAPACK libraries for matrix solutions.

## Output files
The output files are stored in a folder named *case_default* created at runtime. The following files are created:
* *logfile.txt* - This file logs the information about the problem solved such as total degrees of freedom, mesh coordinates, node and edge mappings, integration points used, etc.
* *Plots* - Paraview plot files that contain the numerical (and analytical if available) values over the mesh for each time step.
* *error_data* - If analytical solutions are given, then the normed errors in numerical solution are stored in these files.

## Usage:
All the files should be in the same folder. Open a terminal in the code folder, and type the following for:
* Compilation: `./CleanNCompile`
* Run: `./femSolver`

 The terminal then outputs the time step currently being processed, with the error in numerical solution if the analytical solution is available. 
 
 ## Description of files:
 * *dat* - This is the "su2" format mesh file supplied by the user, generated from Gmsh. The file should be named as "dat" to be read by the solver.
 * *femSolver.f90* - The main code that coordinates the subroutines and functions. This is where you could change some numerical parameters e.g. 
    * the wavenumber and angular frequency of the problem, 
    * number of integration points to be used, 
    * step size in time for finite differences, 
    * total number of timesteps, 
    * number of plots to be stored, etc.
 * *pellib.f90* - This file allows to specify the boundary sources as well as any sources inside the domain.
 * *pellib_DIC.f90* - Modify this file to specify initial conditions.
 * *ln_norm.f90* and *proslib.f90* - These two files are used to specify the analytical solution (if available) for the computation of normed errors and plotting of analytical values over mesh, respectively.
 
 ## Example files:
 An example *dat* file that has a 2D mesh with 4-th order elements is located in the "Example" folder. The Paraview plots of the mesh and an example numerical solution for a progressive plane wave with Neuman boundaries solved over this mesh, are also available. 
