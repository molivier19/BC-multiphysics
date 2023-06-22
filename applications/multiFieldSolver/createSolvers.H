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
    This file create all the solver objects for every region given in 
    the "fieldsDictionary" dictionary. Solver objects are created 
    according to the "physicsSolverDict" dictionary needed in every region 
    folder. This file also creates the multi field iterator that is used i
    to control the outer loop.
    
 Usage
    The "fieldsDictionary" file must exist in the constant directory. 
    It is used to specify the order in which each region will be resolved.
    The content of the file must be as follows:

    ---------------------------------------------------------------------
    regionList n     // n: number of regions
    (
     region0         // first region   
     region1        // Second region
     ...
    );

    nOuterCorr 200;     // Maximum number of iterations in the outer loop
    nOuterCorrMin 0;    // Minimum number of iterations in the outer loop
    ---------------------------------------------------------------------

    Also a file named "physicsSolverDict" must exist in every region folder in 
    the constant directory. The content of the file must be as follows:

    ---------------------------------------------------------------------
    physicsSolver         <SpecificSolverName>; // Name of the specific region solver
    ---------------------------------------------------------------------
    
\*---------------------------------------------------------------------------*/

multiFieldIterator iterator
(
    IOobject
    (
        "iterator",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary fieldsDictionary
(
    IOobject
    (
        "fieldsDictionary",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

List<word> regionList(fieldsDictionary.lookup("regionList"));

int nFields = regionList.size();


Info<< "Create solvers for time = "
    << runTime.timeName() << nl << endl;

autoPtr<physicsSolver> solverPtr[nFields];

for(int i = 0; i < nFields; i++)
{
    word fileDict = regionList[i];
    if(regionList[i] == polyMesh::defaultRegion)
    {
        fileDict = "";
    }
    Info << "Region: " << regionList[i] << endl;
    
    solverPtr[i] = physicsSolver::New
    (
        IOobject
        (
            "physicsSolverDict",
            "constant",
            fileDict,
            runTime,
            IOobject::MUST_READ
        ),
        regionList[i]
    );
}
