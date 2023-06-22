/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 2022 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    multiFieldSolver

Description
    This application was created to manage multiphysics simulations with 
    multiple regions. The number of regions is defined by the user on a
    case-by-case basis as well as the specific field solvers (physicsSolver
    objects).  

\*---------------------------------------------------------------------------*/

//#include "fvCFD.H"
#include "Time.H"
#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"


#include "physicsSolver.H"
#include "multiFieldIterator.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "fvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
//    #include "createTimeControls.H"
    #include "createSolvers.H"
    #include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Residuals display
//    lduMatrix::debug = 0;

    while (runTime.run())
    {

        Info << endl
             << "##############################"
             << "##############################"
             << endl;

        // Pre-solution operations
        for(int i = 0; i < nFields; i++)
        {
            solverPtr[i]->startTimeStep();
        }

        // Solution
        #include "adaptDeltaT.H"
        runTime++;

        Info << "Time: " << runTime.value() << nl << endl;

        #include "readSolverControls.H"

        iterator.init();
        bool isConverged(true);
        do
        {
            isConverged = true;
            for(int i = 0; i < nFields; i++)
            {
                solverPtr[i]->solveFields();
                isConverged = solverPtr[i]->isConverged() && isConverged;
            }

            Info << endl
                 << "------------------------------"
                 << "------------------------------"
                 << endl;

            iterator++;
        }
        while
        (
            (!isConverged || iterator.index() < NOuterCorrMin)
         && iterator.index() < NOuterCorr
        );


        // Post-solution operations
        for(int i = 0; i < nFields; i++)
        {
            solverPtr[i]->endTimeStep();
        }

        runTime.write();

        Info << endl
             << "#################"
             << " End of multi-field loop #"
             << "#################" << endl << endl
             << "No outer-iteration : " << iterator.index() << endl;

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
