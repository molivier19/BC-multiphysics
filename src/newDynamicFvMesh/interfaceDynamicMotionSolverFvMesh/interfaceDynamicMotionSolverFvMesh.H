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

Class
    Foam::interfaceDynamicMotionSolverFvMesh

Description
    The "interfaceDynamicMotionSolverFvMesh" class is used to activate a dynamic
    mesh solver on the first few outer iterations ("nActiveOuterIterations") of
    a multi-field simulation. Only boundary points are moved on subsequent outer
    iterations. This is particularly useful in FSI simulations using the
    artificial compressibility method.

SourceFiles
    interfaceDynamicMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceDynamicMotionSolverFvMesh_H
#define interfaceDynamicMotionSolverFvMesh_H

#include "dynamicMotionSolverFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
            Class interfaceDynamicMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class interfaceDynamicMotionSolverFvMesh
:
    public dynamicMotionSolverFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceDynamicMotionSolverFvMesh
        (
            const interfaceDynamicMotionSolverFvMesh&
        );

        //- Disallow default bitwise assignment
        void operator=(const interfaceDynamicMotionSolverFvMesh&);


public:

    //- Runtime type information
    TypeName("interfaceDynamicMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        interfaceDynamicMotionSolverFvMesh(const IOobject& io);


    //- Destructor
    ~interfaceDynamicMotionSolverFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
