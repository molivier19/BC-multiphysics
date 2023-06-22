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
    This file calls the "computeDeltaT" function in every and then call the 
    "returnDeltaT" function. With the evaluated time step from every region, 
    the algorithm can decide which time step will be used for the next 
    iteration of the time loop.    

\*---------------------------------------------------------------------------*/

bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

if (adjustTimeStep && runTime.timeIndex() > 0)
{
    scalar maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    /* maxCo = */
    /*     runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0); */

    /* // Region from which the "phi" field used to compute the new time-step is */
    /* // taken. */
    /* word regionName */
    /* ( */
    /*     runTime.controlDict(). */
    /*         lookupOrDefault<word>("adjustTimeStepRegion", "region0") */
    /* ); */
    
    /* const fvMesh& mesh(runTime.lookupObject<fvMesh>(regionName)); */
        
    /* const surfaceScalarField& phi = */  
    /*     mesh.lookupObject<surfaceScalarField>("phi"); */

    /* #include "CourantNo.H" */

    /* scalar maxDeltaTFact = maxCo/(CoNum + SMALL); */
    /* scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2); */

    /* //return deltaTFact*runTime.deltaTValue(), */

    bool computeDeltaT(false);
    scalar reqDeltaT = VGREAT;
    for(int i = 0; i < nFields; i++)
    {
        if(solverPtr[i]->computeDeltaT())
        {
            computeDeltaT = true;
            scalar dT = solverPtr[i]->returnDeltaT();

            Info<< "Requested time-step in region " << regionList[i]
                << ": " << dT << endl;

            reqDeltaT = min(reqDeltaT, dT);
        }
    }

    if(computeDeltaT)
    {
        scalar maxDeltaTFact = reqDeltaT/(runTime.deltaT().value() + SMALL);

        // Limit growth rate of the time-step
        scalar deltaT = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2)
            *runTime.deltaT().value();

        runTime.setDeltaT
        (
            min
            (
                deltaT,
                maxDeltaT
            )
        );
    }
    else
    {
        WarningIn(args[0])
            << "\n"
            << "    Adaptive time steping requested in controlDict,\n"
            << "    but no solver is configured to compute deltaT.\n" << endl;
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
