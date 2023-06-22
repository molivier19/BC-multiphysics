# Heated Cylinder

This simulation represents the case of a hypersonic flow around a cylinder. The
aerodynamic heating of the fluid interacts with the cylinder.

Two types of interface boundary conditions are tested. The Dirichlet-Neumann
boundary conditions combination is used in the case "DirichletNeumann", while
mixed boundary conditions used for the case "mixed".

Requirements:

    - HiSA
    - OpenFOAM-extensions

To run a simulation, launch the "Allrun" script present in each simulation
folder. To clean the simulation folder, use the "Allclean" script. To open the
simulation in paraview, use the "paraFoam" command or load the .foam file.

The user may want to adjust the Courant number for the simulations or the
number of processors for parallel computing.  To do so, the user can access the
file "./system/include/CaseControl" and change the "MaxCourantNo" parameters.
Similarly, the user can access the same file to change the "nProcess"
parameters.

More details on this simulation are provided in section 3.2 of the article
entitled "Modular Framework for the Solution of Boundary-Coupled Multiphysics
Problems". Specifically, the two cases were used to produce Figure 9. Also, raw
data and gnuplot scripts to reproduce Figure 9 are included. Note that the
solution is very sensitive to initial conditions and the initial time-step
size. The simulation cases provided have been adapted to run in a reasonable
time frame by increasing the time-step size and by limiting the number of outer
iterations to 5. Hence, the results will slightly differ from the paper.

Each simulation took about 1 hour to complete on an Intel i7-7700 HQ CPU using
4 cores.
