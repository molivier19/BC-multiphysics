# Flexible Tube

This simulation represents the case of flow driven by a pressure pulse within a
flexible tube.  Two types of boundary conditions are used for this particular
case. A step total-pressure boundary condition is used for the simulation
contained in the "Step" folder while a smooth total-pressure boundary
condition is used for the simulation contained in the "Smooth" folder.

Requirements:

    - GMSH
    - OpenFOAM-extensions
    - Compile special boundary condition ("smoothPulseTotalPressure" for
      the "Smooth" case and "pulseTotalPressure" for the "Step" case)

Launch the "Allrun" scripts to run the simulations.  To clean a simulation
folder, use the "Allclean" script.  To open the simulations in Paraview, use
the "paraFoam" command or load the .foam file as usual.

The user may want to adjust the mesh size for a simulation or the number of
processors for parallel computing.  To do so, the user can access the file
"./Add-On/meshes/flexibleTubeParameters.geo" and change the "refinementFactor"
parameter.  Similarly, the user can access the file "./system/NumProcessor" to
change the "NumProcess" parameter.

During the simulation, pressure and flow rate data are extracted and can be
accessed in the postProcessing folder.

Section 3.1 of the article "Modular Framework for the Solution of
Boundary-Coupled Multiphysics Problems" provides more details on this
simulation. The pressure and flow rate data can be used to reproduce results
from Figures 6 ("smooth" case) and 7 ("step" case) of the paper. Also, raw data
and gnuplot scripts to create both figures are included. 

Each simulation with the coarse mesh setting (refinement factor of 2) took less
than an hour to complete on an Intel i7-7700 HQ CPU using 4 cores.
