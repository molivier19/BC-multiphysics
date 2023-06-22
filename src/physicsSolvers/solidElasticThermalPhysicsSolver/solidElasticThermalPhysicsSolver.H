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
    solidElasticThermal <-- physicsSolver

Description
    Object solvers for solid elastic materials with the Lagrangian formulation. 
    Thermal effect is enabled in the "thermalProperties" directory. The thermal 
    coupling in the energy equation and in the Cauchy tensor constitutive law 
    can be disabled in the "thermalProperties" directory. 
    
    The radiative impact on the solid is still a work in progress...

SourceFiles
    solidElasticThermalRadiativePhysicsSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidElasticStrongThermalRadiationPhysicsSolver_H
#define solidElasticStrongThermalRadiationPhysicsSolver_H

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

/*---------------------------------------------------------------------------*\
                           Class solidElasticThermalPhysicsSolver Declaration
\*---------------------------------------------------------------------------*/

class solidElasticThermalPhysicsSolver
:
    public physicsSolver
{

protected:

// Protected data
    
    // Mesh
    fvMesh mesh_;
    pointMesh pMesh_;

    // Interpolator
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
    bool compactNormalStress_;
    bool thermalCoupling_;

    // Finite volume fields pointers
    //autoPtr<volScalarField> Tptr_;
    autoPtr<volVectorField> Dptr_;
    autoPtr<pointVectorField> pointDptr_;
    autoPtr<volVectorField> Uptr_;
    autoPtr<volSymmTensorField> sigmaptr_;
    autoPtr<volScalarField> sigmaEq_;
    autoPtr<volTensorField> gradDptr_;
    autoPtr<volVectorField> divSigmaExpptr_;

    // Residuals and tolerances for the outer loop
    scalar DOuterConvergenceTolerance_;
    scalar TOuterConvergenceTolerance_;

    // Convergence Tolerance for the solver
    scalar InitialResidual_;
    scalar InitialResidualT_;
    scalar DRelativeTolerance_;
    scalar DconvergenceTolerance_;
    scalar TRelativeTolerance_;
    scalar TconvergenceTolerance_;
    int nCorr_;

// Protected Member Functions

    //- Disallow default bitwise copy construct
    solidElasticThermalPhysicsSolver(const solidElasticThermalPhysicsSolver&);

    //- Disallow default bitwise assignment
    void operator=(const solidElasticThermalPhysicsSolver&);

public:

    //- Runtime type information
        TypeName("solidElasticThermalPhysicsSolver");


    // Constructors

        //- Construct from components
        solidElasticThermalPhysicsSolver
        (
            const IOobject& io, 
            const word& regionName
        );


    // Destructor

        virtual ~solidElasticThermalPhysicsSolver();


    // Member Functions

        virtual void createFields();

        virtual bool isConverged();

        virtual void solveFields();

        virtual void startTimeStep();

        virtual void endTimeStep();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
