# BC-multiphysics

## About

BC-multiphysics is an add-on library for OpenFOAM that provide a framework to
run boundary-coupled multiphysics problems. The library is presented in further
details in:

St-Onge, G. and Olivier, M. "Modular framework for the solution of boundary-coupled
multiphysics problems", OpenFOAM Journal, 2023

## Installation

Required software:

    - OpenFOAM v2006 (Go to:
      "https://www.openfoam.com/news/main-news/openfoam-v20-06" to download
      binary or source code.)
    - HiSA (Go to: "https://gitlab.com/hisa/hisa" to download the source code. The
      code was tested using the master branch at commit 538736bf)
    - BC-multiphysics
    - Paraview (https://www.paraview.org/)
    - GMSH (For Ubuntu users: "sudo apt install gmsh" to download. For others:
      https://gmsh.info/)


Make sure that all OpenFOAM libraries are in the same folder. Also, the folder
containing the hisa library must be named "hisa". The OpenFOAM file tree should
look like this:

```
~/OpenFOAM/
      |-- hisa/
      |-- BC-multiphysics/
```

The OpenFOAM installation may be local (typically in the ~/OpenFOAM folder) or
system-wide (e.g., installed with a binary package). In the latter case, make
sure you have installed the system requirements to compile code.

Before the compilation steps, the following line must be added to the end of
the  "~/.bashrc" file:

    source ~/OpenFOAM/BC-multiphysics/etc/bashrc

Compilation steps:

1. If not done yet, compile or install OpenFOAM v2006.
2. Compile the hisa library with the provided "Allwmake" script.
3. Compile the BC-multiphysics library:
    - go to the root folder (cd ~/OpenFOAM/BC-multiphysics).
    - launch the compilation process with the "./Allwmake" script.

The BC-multiphysics library should now be compiled and ready to use! 

The folder "./BC-multiphysics/tutorials/" contains all three test cases
that are discussed in the paper.
