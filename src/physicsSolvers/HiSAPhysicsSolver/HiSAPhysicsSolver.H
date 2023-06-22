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
    HiSAPhysicsSolver <-- physicsSolver

Description
    Adaptation of the HiSA solver for the "multiFieldSolver" application. 
    Adapted and build with respect to the library HiSA version "hisa 1.4.8" 
    for OpenFOAM-v2006. The source of the HISA solver can be found here:

    https://hisa.gitlab.io
    https://gitlab.com/hisa/hisa

SourceFiles
    HiSAPhysicsSolver.C

\*---------------------------------------------------------------------------*/

#ifndef HiSAPhysicsSolver_H
#define HiSAPhysicsSolver_H

#include "physicsSolver.H"
#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "pimpleControl.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
//#include "fixedRhoFvPatchScalarField.H"
#include "fluxScheme.H"
#include "pseudotimeControl/pseudotimeControl.H"
#include "jacobian.H"
#include "upwind.H"
#include "preconditioner.H"
#include "faceSet.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
//class someClass;

/*---------------------------------------------------------------------------*\
                           Class HiSAPhysicsSolver Declaration
\*---------------------------------------------------------------------------*/

class HiSAPhysicsSolver
:
    public physicsSolver
{

protected:

    // Ordering of residual printed output

        const label residualOrderingArray[3] = {0, 2, 1};

        const labelList residualOrdering =
            labelList(FixedList<label,3>(residualOrderingArray));

    // Private data

        /* Check if reference to time variable affect the program
        //- Reference to the time database
        const Time& time_;
        */

        //Mesh
        autoPtr<dynamicFvMesh> meshPtr_;
        dynamicFvMesh&  mesh_;

        autoPtr< pseudotimeControl > solnControl_;

        autoPtr<psiThermo> pThermo_;

        autoPtr<volVectorField> U_;

        // Stored (if viscous) for use in maxwellSlip boundary condition
        autoPtr<volTensorField> tauMC_;

        autoPtr<surfaceScalarField> phi_;
        autoPtr<surfaceVectorField> phiUp_;
        autoPtr<surfaceScalarField> phiEp_;
        autoPtr<surfaceVectorField> Up_;

        // The value of thermo.he() at T = 0 and p = 0, used as offset to get
        // absolute internal energy
        autoPtr<volScalarField> eZero_;

        autoPtr<volScalarField> rPseudoDeltaT_;
        // boolean to check if the simulation has converged
        bool finished_;
        // 1/dt (for Euler), or equivalent for other schemes
        autoPtr<scalarField> ddtCoeff_;
        // Whether variables had to be bounded this iteration
        bool bounded_;

        //rho and rhoE variables
        PtrList<volScalarField> scalarVars_;
        PtrList<volScalarField> scalarVarsPrevIter_;
        //rhoU variable
        PtrList<volVectorField> vectorVars_;
        PtrList<volVectorField> vectorVarsPrevIter_;

        //rho and rhoE residuals
        PtrList<volScalarField> scalarResiduals_;
        //rhoU residual
        PtrList<volVectorField> vectorResiduals_;

        // fluid stress for fluid structure implementation
        autoPtr<volSymmTensorField> sigmaFluidptr_;

        autoPtr<compressible::turbulenceModel> turbulence_;

        autoPtr<fluxScheme> flux_;

        bool inviscid_;

        //- How frequently mesh is moved during outer corrector loop.
        //  0 means do not move (only at beginning of timestep).
        //  Only applies to transient analysis.
        label moveMeshOuterCorrectors_;

        bool steadyState_;
        bool localTimestepping_;
        bool localTimesteppingBounding_;
        scalar localTimesteppingLowerBound_;
        scalar localTimesteppingUpperBound_;

        //- Store previous residual for CoNum relaxation
        autoPtr<residualIO> initRes_;
        autoPtr<residualIO> prevRes_;

        //- Current pseudo-courant number
        autoPtr< uniformDimensionedScalarField > pseudoCoNum_; // global timestepping
        autoPtr< volScalarField > pseudoCoField_;              // local timestepping
        //- Minimum pseudo-Courant number
        scalar pseudoCoNumMin_;
        //- Maximum pseudo-Courant number
        scalar pseudoCoNumMax_;
        //- Maximum increase as a fraction of current pseudo-Courant number
        scalar pseudoCoNumMaxIncr_;
        //- Minimum decrease as a fraction of current pseudo-Courant number
        scalar pseudoCoNumMinDecr_;

        //- Print debug info for nearest cell
        bool cellDebugging_;
        label debugCell_;

        //- Whether to do parallel redistribute on mesh topo change
        bool rebalance_;
        //- Maximum load imbalance as a fraction before rebalancing done
        scalar maxLoadImbalance_;

    // Protected Member Functions

        /* TODO - Check if related to physicsSolver or rhoCentralFoam
        //- Interpolate field vf according to direction dir
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& dir,
            const word& reconFieldName = word::null
        );*/

        //- Set timestep size based on current Courant no.
        void setPseudoDeltaT();

        //- Adjust Courant number
        void setPseudoCoNum();

        //- Recursively create preconditioners and their Jacobians if applicable
        void createPreconditioners
        (
            PtrList<preconditioner<2,1>>& preconditioners,
            PtrList<jacobian>& jacobians,
            const dictionary& parentDict
        );

        //- Helper function to synchronise fields after mesh redistribute
        template<class FieldType, class Type>
        void parallelSyncFields(const wordList& fields);

        //- Check load-balance and re-distribute mesh if necessary
        void redistributePar();

        //- Search for debug point specified in controlDict
        void findDebugCell();

        //- debug cell 
        void cellDebug();

        //- Disallow default bitwise copy construct
        HiSAPhysicsSolver(const HiSAPhysicsSolver&);

        //- Disallow default bitwise assignment
        void operator=(const HiSAPhysicsSolver&);

public:

    //- Runtime type information
        TypeName("HiSAPhysicsSolver");


    // Constructors

        //- Construct from components
        HiSAPhysicsSolver
        (
            const IOobject& io, 
            const word& regionName
        );


    // Destructor

        virtual ~HiSAPhysicsSolver();


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
