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
    fluidPhysicsSolver <-- physicsSolver

Description
    Class solver for incompressible flows. The solution method is 
    similar to the PIMPLE formulation.

SourceFiles
    fluidPhysicsSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fluidPhysicsSolver_H
#define fluidPhysicsSolver_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "physicsSolver.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
//class someClass;

/*---------------------------------------------------------------------------*\
                           Class fluidPhysicsSolver Declaration
\*---------------------------------------------------------------------------*/

class fluidPhysicsSolver
:
    public physicsSolver
{

protected:

    // Protected data

        // Mesh
        autoPtr<dynamicFvMesh> meshPtr_;
        dynamicFvMesh& mesh_;

        // Properties
        IOdictionary transportProperties_;
        dimensionedScalar nu_;
        dimensionedScalar rho_;

        // Courant number
        scalar CoNum_;

        // Finite volume fields pointers
        autoPtr<volScalarField> pptr_;
        autoPtr<volVectorField> Uptr_;
        autoPtr<volSymmTensorField> sigmaFluidptr_;
        autoPtr<surfaceScalarField> phiptr_;

        label pRefCell_;
        scalar pRefValue_;
        scalar totalVolume_;
        scalar cumulativeContErr_;

        // Residuals and tolerances
        scalar UInitialResidual_;
        scalar pInitialResidual_;

        scalar UConvergenceTolerance_;
        scalar URelativeTolerance_;
        scalar pConvergenceTolerance_;
        scalar pRelativeTolerance_;

        scalar UOuterConvergenceTolerance_;
        scalar pOuterConvergenceTolerance_;

    // Protected Member Functions

        //- Disallow default bitwise copy construct
        fluidPhysicsSolver(const fluidPhysicsSolver&);

        //- Disallow default bitwise assignment
        void operator=(const fluidPhysicsSolver&);


public:

    //- Runtime type information
        TypeName("fluidPhysicsSolver");


    // Constructors

        //- Construct from components
        fluidPhysicsSolver
        (
            const IOobject& io, 
            const word& regionName
        );


    // Destructor

        virtual ~fluidPhysicsSolver();


    // Member Functions

        virtual void createFields();

        virtual bool isConverged();

        virtual void solveFields();

        virtual void startTimeStep();

        virtual bool computeDeltaT();

        virtual scalar returnDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
