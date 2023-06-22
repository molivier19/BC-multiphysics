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
    solidThermalStVKPhysicsSolver <-- physicsSolver

Description
    Object solvers for solid elastic materials with the Total Lagrangian formulation.
    Non-linear geometric deformation is taken into account by the Saint-Venant-Kirchoff
    formulation of the generalized Hooke's law. The solver is formulated in the reference
    configuration as stated by the total Lagrangian formulation. Thermal effect is enabled 
    in the "thermalProperties" directory. The thermal coupling in the energy equation and 
    in the Cauchy tensor constitutive law can be disabled in the "thermalProperties" directory.

    The radiative impact on the solid is still a work in progress...

SourceFiles
    SolidThermalStVKPhysicsSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermalStVKPhysicsSolver_H
#define solidThermalStVKPhysicsSolver_H

#include "fvCFD.H"
#include "physicsSolver.H"
#include "pointMesh.H"
#include "pointFields.H"
//#include "gradBasedVolPointInterpolation.H"
//#include "volPointInterpolation.H"
#include "solidThermo.H"
#include "motionInterpolation.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
//class someClass;

/*---------------------------------------------------------------------------*\
                           Class solidThermalStVKPhysicsSolver Declaration
\*---------------------------------------------------------------------------*/

class solidThermalStVKPhysicsSolver
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

        // dictionnary
        const dictionary solverControls_;

        // Properties
        IOdictionary mechanicalProperties_;
        dimensionedScalar E_;
        dimensionedScalar nu_;
        dimensionedScalar mu_;
        dimensionedScalar lambda_;
        dimensionedScalar threeK_;
        dimensionedScalar alpha_;
        dimensionedScalar T0_;
        dimensionedScalar threeKalpha_;

        // Thermal Properties
        autoPtr<solidThermo> sThermoptr_;

        // bool control
        bool thermalStress_;
        bool thermalStressOnly_;
        bool planeStress_;
        bool thermalCoupling_;

        // Finite volume fields pointers
        autoPtr<volVectorField> Dptr_;
        autoPtr<pointVectorField> pointDptr_;
        autoPtr<volVectorField> Uptr_;
        autoPtr<volScalarField> Tptr_;
        autoPtr<volTensorField> gradDptr_;
        autoPtr<volSymmTensorField> sigmaptr_;
        autoPtr<volSymmTensorField> sigmaCauchyptr_;
        autoPtr<volScalarField> sigmaEqptr_;
        autoPtr<volVectorField> DPrevOuterIterptr_;

        // Residuals and tolerances
        scalar DInitialResidual_;
        scalar TInitialResidual_;

        scalar DConvergenceTolerance_;
        scalar DRelativeTolerance_;

        scalar TConvergenceTolerance_;
        scalar TRelativeTolerance_;

        scalar DOuterConvergenceTolerance_;
        scalar TOuterConvergenceTolerance_;

    // Protected Member Functions

        void solveImplicit();

        void solveExplicit();

        //- Disallow default bitwise copy construct
        solidThermalStVKPhysicsSolver(const solidThermalStVKPhysicsSolver&);

        //- Disallow default bitwise assignment
        void operator=(const solidThermalStVKPhysicsSolver&);


public:

    //- Runtime type information
        TypeName("solidThermalStVKPhysicsSolver");


    // Constructors

        //- Construct from components
        solidThermalStVKPhysicsSolver
        (
            const IOobject& io, 
            const word& regionName
        );


    // Destructor

        virtual ~solidThermalStVKPhysicsSolver();


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
