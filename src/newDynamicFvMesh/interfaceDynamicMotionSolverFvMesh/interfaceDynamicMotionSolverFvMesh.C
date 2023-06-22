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

\*---------------------------------------------------------------------------*/

#include "interfaceDynamicMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "multiFieldIterator.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceDynamicMotionSolverFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicMotionSolverFvMesh,
        interfaceDynamicMotionSolverFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceDynamicMotionSolverFvMesh::interfaceDynamicMotionSolverFvMesh
(
    const IOobject& io
)
:
    dynamicMotionSolverFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceDynamicMotionSolverFvMesh::~interfaceDynamicMotionSolverFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::interfaceDynamicMotionSolverFvMesh::update()
{
    const multiFieldIterator& mfi = 
        time().lookupObject<multiFieldIterator>("iterator");

    label nActiveOuterIterations
    (
        readLabel
        (
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    time().constant(),
                    *this,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            ).lookup("nActiveOuterIterations")
            // ).lookupOrDefault<label>("nActiveOuterIterations", 1)
        )
    );

    if(nActiveOuterIterations < 1)
    {
        Warning<<"nActiveOuterIterations should be larger or equal to 1." 
            << endl << "Setting nActiveOuterIterations to 1." << endl;
        nActiveOuterIterations = 1;
    }

    if(mfi.index() < nActiveOuterIterations)
    {
        return Foam::dynamicMotionSolverFvMesh::update();
    }
    else
    {
        Info<< "Moving mesh boundaries" << endl;
        const displacementMotionSolver& ms
        (
            lookupObject<displacementMotionSolver>("dynamicMeshDict")
        );

        pointVectorField& pointDisplacement = 
            const_cast<pointVectorField&>(ms.pointDisplacement());
//            const_cast<pointVectorField&>
//            (
//                lookupObject<pointVectorField>("pointDisplacement")
//            );

        pointDisplacement.boundaryFieldRef().updateCoeffs();
        pointDisplacement.boundaryFieldRef().evaluate();

//        const pointField& p0 = ms.points0();

        movePoints(ms.curPoints());

        if (foundObject<volVectorField>("U"))
        {
            volVectorField& U =
                const_cast<volVectorField&>(lookupObject<volVectorField>("U"));
            U.correctBoundaryConditions();
        }
        return true;
    }
}


// ************************************************************************* //
