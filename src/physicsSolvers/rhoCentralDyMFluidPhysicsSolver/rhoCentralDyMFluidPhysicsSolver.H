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
    rhoCentralDyMFluidPhysicsSolver <-- physicsSolver

Description
    Object based solvers adapted from the "rhoCentralDyMFluid" application.
    Adapted from the OpenFOAM-v2006 library.

    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor with support for mesh-motion and topology changes.

SourceFiles
    rhoCentralDyMFluidPhysicsSolver.C

\*---------------------------------------------------------------------------*/

#ifndef rhoCentralDyMFluidPhysicsSolver_H
#define rhoCentralDyMFluidPhysicsSolver_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "physicsSolver.H"
#include "runTimeSelectionTables.H"
// #include "postProcess.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
//class someClass;

/*---------------------------------------------------------------------------*\
                           Class rhoCentralDyMFluidPhysicsSolver Declaration
\*---------------------------------------------------------------------------*/

class rhoCentralDyMFluidPhysicsSolver
:
    public physicsSolver
{

protected:

    // Protected data

        // Mesh
        autoPtr<dynamicFvMesh> meshPtr_;
        dynamicFvMesh&  mesh_;

        // properties
        autoPtr<psiThermo> pThermoptr_;
        bool inviscid_;

        // Courant number
        scalar CoNum_;
        scalar meanCoNum_;

        // constant
        const dimensionedScalar v_zero_;

        // Finite volume fields pointers
        autoPtr<volVectorField> Uptr_;
        autoPtr<volScalarField> rhoptr_;
        autoPtr<volVectorField> rhoUptr_;
        autoPtr<volScalarField> rhoEptr_;

        autoPtr<surfaceScalarField> posptr_;
        autoPtr<surfaceScalarField> negptr_;
        autoPtr<surfaceScalarField> phiptr_;

        autoPtr<volScalarField> divPhi_;
        autoPtr<volVectorField> divPhiUp_;
        autoPtr<volScalarField> divEnergy_;

        autoPtr<volSymmTensorField> sigmaFluidptr_;


        // turbulence model
        autoPtr<compressible::turbulenceModel> turbulence_;

        scalar UInitialResidual_;
        scalar hInitialResidual_;

        scalar UConvergenceTolerance_;
        scalar URelativeTolerance_;
        scalar hConvergenceTolerance_;
        scalar hRelativeTolerance_;

        // flux Scheme
        word fluxScheme_;

    // Protected Member Functions

        //- Interpolate field vf according to direction dir
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& dir,
            const word& reconFieldName = word::null
        );

        //- Compute explicit fluxes
        void computeExplicitFluxes();

        //- Disallow default bitwise copy construct
        rhoCentralDyMFluidPhysicsSolver(const rhoCentralDyMFluidPhysicsSolver&);

        //- Disallow default bitwise assignment
        void operator=(const rhoCentralDyMFluidPhysicsSolver&);

        //- Solve the fluid
        // void solveFluidFields();

public:

    //- Runtime type information
        TypeName("rhoCentralDyMFluidPhysicsSolver");


    // Constructors

        //- Construct from components
        rhoCentralDyMFluidPhysicsSolver
        (
            const IOobject& io, 
            const word& regionName
        );


    // Destructor

        virtual ~rhoCentralDyMFluidPhysicsSolver();


    // Member Functions

        virtual void createFields();
    
        virtual bool isConverged();

        virtual void solveFields();

        virtual void startTimeStep();

        virtual void endTimeStep();

        virtual bool computeDeltaT();

        virtual scalar returnDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
