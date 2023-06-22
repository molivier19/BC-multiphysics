# Cross-flow flexible plate

This simulation represents the case of a flexible plate interacting with a
transverse hypersonic flow. The aerodynamic forces and aerodynamic heating of
the fluid interact with the plate.

Requirement:

    - HiSA
    - OpenFOAM-extensions

To run a simulation, launch the "Allrun" script present in each simulation
folder.  To clean the simulation folder, use the "Allclean" script.  To open
the simulation in paraview, use the "paraFoam" command or open the .foam file
with paraview.

The user may want to adjust the courant number for a simulation or the number
of processors for parallel computing.  To do so, the user can access the file
"./system/include/CaseControl" and change the "MaxCourantNo" parameter.
Similarly, the user can access the same file to change the "nProcess" parameter
to change the number of processors used for the simulation.

More details on this simulation are provided in section 3.3 of the article
entitled "Modular Framework for the Solution of Boundary-Coupled Multiphysics
Problems".

The simulation took under 3 hours to complete on an Intel i7-7700 HQ CPU
using 4 cores.
