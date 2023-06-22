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
    solidStVKPhysicsSolver <-- physicsSolver

Description
    Class solver for solid elastic materials with a total Lagrangian formulation
    based on the initial configuration. Non-linear geometric deformation is taken
    into account and the Saint-Venant-Kirchoff constitutive law is used.

SourceFiles
    SolidStVKPhysicsSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidStVKPhysicsSolver_H
#define solidStVKPhysicsSolver_H

#include "fvCFD.H"
#include "physicsSolver.H"
#include "pointMesh.H"
#include "pointFields.H"
//#include "gradBasedVolPointInterpolation.H"
//#include "volPointInterpolation.H"
#include "motionInterpolation.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
//class someClass;

/*---------------------------------------------------------------------------*\
                           Class solidStVKPhysicsSolver Declaration
\*---------------------------------------------------------------------------*/

class solidStVKPhysicsSolver
:
    public physicsSolver
{

protected:

    // Protected data

        // Mesh
        fvMesh mesh_;
        pointMesh pMesh_;

        // Interpolator
//        gradBasedVolPointInterpolation cpi_;
//        volPointInterpolation cpi_;
        autoPtr<motionInterpolation> cpi_;

        // Properties
        IOdictionary mechanicalProperties_;
        dimensionedScalar rho_;
        dimensionedScalar E_;
        dimensionedScalar nu_;
        dimensionedScalar mu_;
        dimensionedScalar lambda_;
        dimensionedScalar threeK_;

        // Finite volume fields pointers
        autoPtr<volVectorField> Dptr_;
        autoPtr<pointVectorField> pointDptr_;
        autoPtr<volVectorField> Uptr_;
        autoPtr<volTensorField> gradDptr_;
        autoPtr<volSymmTensorField> sigmaptr_;
        autoPtr<volVectorField> DPrevOuterIterptr_;

        // Residuals and tolerances
        scalar DInitialResidual_;

        scalar DConvergenceTolerance_;
        scalar DRelativeTolerance_;

        scalar DOuterConvergenceTolerance_;

    // Protected Member Functions

        void solveImplicit();

        void solveExplicit();

        //- Disallow default bitwise copy construct
        solidStVKPhysicsSolver(const solidStVKPhysicsSolver&);

        //- Disallow default bitwise assignment
        void operator=(const solidStVKPhysicsSolver&);


public:

    //- Runtime type information
        TypeName("solidStVKPhysicsSolver");


    // Constructors

        //- Construct from components
        solidStVKPhysicsSolver
        (
            const IOobject& io, 
            const word& regionName
        );


    // Destructor

        virtual ~solidStVKPhysicsSolver();


    // Member Functions

        virtual void createFields();

        virtual bool isConverged();

        virtual void solveFields();

        virtual void endTimeStep();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
