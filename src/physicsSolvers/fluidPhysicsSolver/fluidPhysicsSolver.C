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

#include "fluidPhysicsSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fluidPhysicsSolver, 0);

addToRunTimeSelectionTable(physicsSolver, fluidPhysicsSolver, IOobject);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPhysicsSolver::fluidPhysicsSolver
(
    const IOobject& io, 
    const word& regionName
)
:
    physicsSolver(io, regionName),

    meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region_,
                io.time().timeName(),
                io.time(),
                IOobject::MUST_READ
            )
        )
    ),

    mesh_(meshPtr_()),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            io.time().constant(),
            meshPtr_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    nu_("nu", dimViscosity, transportProperties_),

    rho_("rho", dimDensity, transportProperties_),

    CoNum_(0),

    pptr_(nullptr),
    Uptr_(nullptr),
    sigmaFluidptr_(nullptr),
    phiptr_(nullptr),

    pRefCell_(0),
    pRefValue_(0.0),
    totalVolume_(sum(meshPtr_->V()).value()),

    cumulativeContErr_(0.0),
    UInitialResidual_(0.0),
    pInitialResidual_(0.0),

    UConvergenceTolerance_(0.0),
    URelativeTolerance_(0.0),
    pConvergenceTolerance_(0.0),
    pRelativeTolerance_(0.0),
    UOuterConvergenceTolerance_(0.0),
    pOuterConvergenceTolerance_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhysicsSolver::~fluidPhysicsSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidPhysicsSolver::createFields()
{
    pptr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Uptr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );
    
    sigmaFluidptr_.reset
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigmaFluid",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // No shear stress to avoid very large force on 
            // the solid on the first timestep when the 
            // startup is impulsive.
            //-rho_*(pptr_()*symmTensor(1,0,0,1,0,1))

            rho_*
            (
                (-pptr_()*symmTensor(1,0,0,1,0,1))
              + (nu_*2*symm(fvc::grad(Uptr_())))
            )
        )
    );
    
    phiptr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(Uptr_()) & mesh_.Sf()
        )
    );

    setRefCell
    (
        pptr_(),
        mesh_.solutionDict().subDict("solverControls"), pRefCell_,
        pRefValue_
    );

    // Fluid solver controls
    const dictionary& solverControls = 
        mesh_.solutionDict().subDict("solverControls");
    UConvergenceTolerance_ = readScalar(solverControls.lookup("U"));
    URelativeTolerance_ = readScalar(solverControls.lookup("relU"));
    pConvergenceTolerance_ = readScalar(solverControls.lookup("p"));
    pRelativeTolerance_ = readScalar(solverControls.lookup("relp"));
    UOuterConvergenceTolerance_ = readScalar(solverControls.lookup("UOuter"));
    pOuterConvergenceTolerance_ = readScalar(solverControls.lookup("pOuter"));
}


bool Foam::fluidPhysicsSolver::isConverged()
{
    if
    (
        UInitialResidual_ < UOuterConvergenceTolerance_
     && pInitialResidual_ < pOuterConvergenceTolerance_
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fluidPhysicsSolver::solveFields()
{
    volScalarField& p = pptr_();
    volVectorField& U = Uptr_();
    volSymmTensorField& sigmaFluid = sigmaFluidptr_();
    surfaceScalarField& phi = phiptr_();

    UInitialResidual_ = 0.0;
    pInitialResidual_ = 0.0;

    // These references are created so that it is 
    // possible to use standard OpenFOAM include 
    // files such as CourantNo.H, etc.
    const Time& runTime = mesh_.time();
    dynamicFvMesh& mesh = mesh_;
    scalar& cumulativeContErr = cumulativeContErr_;


    // Fluid solver controls
    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");
    int nPIMPLECorr(readInt(solverControls.lookup("nCorrectors")));
    int nPISOCorr(readInt(solverControls.lookup("nPISOCorrectors")));
    int nNonOrthCorr
    (
        readInt(solverControls.lookup("nNonOrthogonalCorrectors"))
    );
    UConvergenceTolerance_ = readScalar(solverControls.lookup("U"));
    URelativeTolerance_ = readScalar(solverControls.lookup("relU"));
    pConvergenceTolerance_ = readScalar(solverControls.lookup("p"));
    pRelativeTolerance_ = readScalar(solverControls.lookup("relp"));
    UOuterConvergenceTolerance_ = readScalar(solverControls.lookup("UOuter"));
    pOuterConvergenceTolerance_ = readScalar(solverControls.lookup("pOuter"));

    // PIMPLE loop
    mesh_.update();

    scalar UResidual = VGREAT;
    scalar pResidual = VGREAT;
    int corr = 0;

    do
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu_, U)
        );
        UEqn.relax();

        SolverPerformance<vector> USolPerf;
//        solverPerformance USolPerf;
        USolPerf = solve(UEqn == -fvc::grad(p));
        UResidual = cmptMax(USolPerf.initialResidual());

        p.storePrevIter();

        // --- PISO loop
        for (int PISOCorr=0; PISOCorr<nPISOCorr; PISOCorr++)
        {
            U = UEqn.H()/UEqn.A();

            phi = fvc::interpolate(U) & mesh_.Sf();

            if (p.needReference())
            {
                fvc::makeRelative(phi, U);
                adjustPhi(phi, U, p);
                fvc::makeAbsolute(phi, U);
            }

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/UEqn.A(), p) == fvc::div(phi)
                );
                pEqn.setReference(pRefCell_, pRefValue_);

                if(PISOCorr == 0 && nonOrth == 0)
                {
                    solverPerformance pSolPerf;
                    pSolPerf = pEqn.solve();
                    pResidual = pSolPerf.initialResidual();
                }
                else
                {
                    pEqn.solve();
                }

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
                if(corr == 0 && PISOCorr == 0 && nonOrth == 0)
                {
                    UInitialResidual_ = UResidual;
                    pInitialResidual_ = pResidual;
                }
            }

            p.relax();
            U -= fvc::grad(p)/UEqn.A();
            U.correctBoundaryConditions();
        }

        // Make the fluxes relative
        //phi -= fvc::meshPhi(U);
        fvc::makeRelative(phi, U);

    } 
    while 
    (
        (++corr < nPIMPLECorr)
        && 
        (
            (
                UResidual/(UInitialResidual_ + SMALL)
              > URelativeTolerance_
            ) 
            ||
            (
                pResidual/(pInitialResidual_ + SMALL)
              > pRelativeTolerance_
            )
        )
        && 
        (
            (UResidual > UConvergenceTolerance_) 
         || (pResidual > pConvergenceTolerance_)
        )
    );
#   include "movingMeshContinuityErrs.H"

    sigmaFluid = rho_*
    (
        (-p*symmTensor(1,0,0,1,0,1))
      + (nu_*twoSymm(fvc::grad(U)))
    );
    //sigmaFluid = -rho_*(p*symmTensor(1,0,0,1,0,1));

    Info << endl << "Fluid Solver Monitoring Summary" << endl
         << "U initial residual : " << UInitialResidual_ << endl
         << "U final residual : " << UResidual << endl
         << "p initial residual : " << pInitialResidual_ << endl
         << "p final residual : " << pResidual << endl
         << "No iteration : " << corr << endl;
}


void Foam::fluidPhysicsSolver::startTimeStep()
{
    surfaceScalarField& phi = phiptr_();
    const Time& runTime = mesh_.time();
    dynamicFvMesh& mesh = mesh_;

    #include "CourantNo.H"
    CoNum_ = CoNum;
}


bool Foam::fluidPhysicsSolver::computeDeltaT()
{
    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    bool computeDeltaT = 
        solverControls.lookupOrDefault<bool>("computeDeltaT", false);

    return computeDeltaT;
}


scalar Foam::fluidPhysicsSolver::returnDeltaT()
{
    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    scalar maxCo(readScalar(solverControls.lookup("maxCo")));

    return time().deltaT().value()*maxCo/(CoNum_ + SMALL);
}
// ************************************************************************* //
