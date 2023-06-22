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

#include "solidStVKPhysicsSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "steadyStateDdtScheme.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "multiFieldIterator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidStVKPhysicsSolver, 0);

addToRunTimeSelectionTable(physicsSolver, solidStVKPhysicsSolver, IOobject);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solidStVKPhysicsSolver::solveImplicit()
{
    volVectorField& D = Dptr_();
    //pointVectorField& pointD = pointDptr_();
    volVectorField& U = Uptr_();
    volTensorField& gradD = gradDptr_();
    volSymmTensorField& sigma = sigmaptr_();
    //volVectorField& DPrevOuterIter = DPrevOuterIterptr_();

    // Solver controls
    const dictionary& solverControls = mesh_.solutionDict().subDict("solverControls");
    int nCorr = readInt(solverControls.lookup("nCorrectors"));
    DConvergenceTolerance_ = readScalar(solverControls.lookup("D"));
    DRelativeTolerance_ = readScalar(solverControls.lookup("relD"));

#   include "solidStVKDdtCoeffs.H"

    scalar DResidual = VGREAT;
    int corr = 0;
    DInitialResidual_ = 0;

    do
    {
        D.storePrevIter();
           
        fvVectorMatrix DEqn
        (
            rho_*(Cn*fvm::ddt(D) - Co*U.oldTime() + Coo*U.oldTime().oldTime())
         ==
            fvm::laplacian(2*mu_ + lambda_, D, "laplacian(DD,D)") 
          - fvc::laplacian(2*mu_ + lambda_, D, "laplacianExp(DD,D)") 
          + fvc::div(sigma & (I + gradD))
        );

        SolverPerformance<vector> DSolPerf;
//        solverPerformance DSolPerf;
        DSolPerf = DEqn.solve();
        DResidual = cmptMax(DSolPerf.initialResidual());
        D.relax();

        if(corr == 0)
        {
            DInitialResidual_ = DResidual;
        }   

        gradD = fvc::grad(D);
        sigma = symm((mu_*(gradD + gradD.T())) 
               + ((lambda_*I)*tr(gradD)) 
               + (mu_*(gradD & gradD.T())) 
               + (0.5*lambda_*I*tr(gradD & gradD.T())));

    } while
    (
        (++corr < nCorr)
     && (DResidual/(DInitialResidual_ + SMALL) > DRelativeTolerance_)
     && (DResidual > DConvergenceTolerance_)
    );

    Info << endl << "Solid Solver Monitoring Summary" << endl
         << "D initial residual : " << DInitialResidual_ << endl
         << "D final residual : " << DResidual << endl
         << "No iteration : " << corr << endl;

    U = fvc::ddt(D);    
}


void Foam::solidStVKPhysicsSolver::solveExplicit()
{
    volVectorField& D = Dptr_();
    volVectorField& U = Uptr_();
    volTensorField& gradD = gradDptr_();
    volSymmTensorField& sigma = sigmaptr_();

    // Required to ajust the load coming from fluid in FSI.
    //D.correctBoundaryConditions();
    gradD = fvc::grad(D);
    sigma = symm((mu_*(gradD + gradD.T())) 
           + ((lambda_*I)*tr(gradD)) 
           + (mu_*(gradD & gradD.T())) 
           + (0.5*lambda_*I*tr(gradD & gradD.T())));

    dimensionedScalar Dn("Cn", dimless/(dimTime*dimTime), 0.0);
    dimensionedScalar Do("Co", dimless/(dimTime*dimTime), 0.0);
    dimensionedScalar Doo("Coo", dimless/(dimTime*dimTime), 0.0);
    dimensionedScalar Uo("Co", dimless/dimTime, 0.0);

    if(mesh_.time().timeIndex() == 1)
    {
        // Explicit central (2nd order)
        // Start-up procedure
        scalar deltaT = mesh_.time().deltaT().value();
        Dn.value() = 2.0/(deltaT*deltaT);
        Do.value() = -2.0/(deltaT*deltaT);
        Uo.value() = -2.0/(deltaT);
    }
    else
    {
        // Explicit central (2nd order)
        // Note: does not support variable timesteps
        scalar deltaT = mesh_.time().deltaT().value();
        Dn.value() = 1.0/(deltaT*deltaT);
        Do.value() = -2.0/(deltaT*deltaT);
        Doo.value() = 1.0/(deltaT*deltaT);
    }
    
    fvVectorMatrix DEqn
    (
        rho_*
        (
            fvm::Sp(Dn,D)
          + Do*D.oldTime()
          + Doo*D.oldTime().oldTime()
          + Uo*U.oldTime()
        )
      - fvc::div(sigma & (I + gradD))
    );
    DEqn.solve();

    gradD = fvc::grad(D);
    sigma = symm((mu_*(gradD + gradD.T())) 
           + ((lambda_*I)*tr(gradD)) 
           + (mu_*(gradD & gradD.T())) 
           + (0.5*lambda_*I*tr(gradD & gradD.T())));

    //U_ = fvc::ddt(D);    
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidStVKPhysicsSolver::solidStVKPhysicsSolver
(
    const IOobject& io, 
    const word& regionName
)
:
    physicsSolver(io, regionName),

    mesh_
    (
        IOobject
        (
            region_,
            io.time().timeName(),
            io.time(),
            IOobject::MUST_READ
        )
    ),

    pMesh_(mesh_),

//    cpi_(mesh_),
    cpi_
    (
        mesh_.solutionDict().found("cellToPointInterpolation")
      ? motionInterpolation::New
        (
            mesh_, mesh_.solutionDict().lookup("cellToPointInterpolation")
        )
      : motionInterpolation::New(mesh_)
    ),

    mechanicalProperties_
    (
        IOobject
        (
            "mechanicalProperties",
            io.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    rho_("rho", dimDensity, mechanicalProperties_),
    E_("E", dimPressure, mechanicalProperties_),
    nu_("nu", dimless, mechanicalProperties_),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    threeK_(E_/(1.0 - 2.0*nu_)),

    Dptr_(nullptr),
    pointDptr_(nullptr),
    Uptr_(nullptr),
    gradDptr_(nullptr),
    sigmaptr_(nullptr),
    DPrevOuterIterptr_(nullptr),

    DInitialResidual_(1.0),

    DConvergenceTolerance_(0.0),
    DRelativeTolerance_(0.0),
    DOuterConvergenceTolerance_(0.0)
{
    // Correct Lamé coefficients for if plane stress
    Switch planeStress(mechanicalProperties_.lookup("planeStress"));

    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        //- change lambda and threeK for plane stress
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));
        threeK_ = E_/(1.0 - nu_);
    }

    Info<< "mu = " << mu_.value() << " Pa\n";
    Info<< "lambda = " << lambda_.value() << " Pa\n";
    Info<< "threeK = " << threeK_.value() << " Pa\n";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidStVKPhysicsSolver::~solidStVKPhysicsSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidStVKPhysicsSolver::createFields()
{
    Dptr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "D",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    gradDptr_.reset
    (
        new volTensorField
        (
            IOobject
            (
                "gradD",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(Dptr_())
        )
    );

    pointDptr_.reset
    (
        new pointVectorField
        (
            IOobject
            (
                "pointD",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
//            cpi_.interpolate(Dptr_(),gradDptr_())
//            cpi_.interpolate(Dptr_())
            pointMesh::New(mesh_),
            dimensionedVector("", vector::zero)
        )
    );
    cpi_->interpolate(Dptr_(), pointDptr_());
    pointDptr_->storePrevIter();

    Uptr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::ddt(Dptr_())
            //mesh_,
            //dimensionedVector("", dimVelocity, vector::zero)
        )
    );

    sigmaptr_.reset
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            symm((mu_*(gradDptr_() + gradDptr_().T())) 
            + ((lambda_*I)*tr(gradDptr_())) 
            + (mu_*(gradDptr_() & gradDptr_().T())) 
            + (0.5*lambda_*I*tr(gradDptr_() & gradDptr_().T())))
        )
    );

    DPrevOuterIterptr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "DPrevOuterIter",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Dptr_()
        )
    );

    // Solver controls
    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");
    //int nCorr = readInt(solverControls.lookup("nCorrectors"));
    DConvergenceTolerance_ = readScalar(solverControls.lookup("D"));
    DRelativeTolerance_ = readScalar(solverControls.lookup("relD"));
    //DOuterConvergenceTolerance_ = readScalar(solverControls.lookup("DOuter"));
}


bool Foam::solidStVKPhysicsSolver::isConverged()
{
    const dictionary& solverControls = 
        mesh_.solutionDict().subDict("solverControls");
    DOuterConvergenceTolerance_ = readScalar(solverControls.lookup("DOuter"));
    if
    (
        DInitialResidual_ < DOuterConvergenceTolerance_
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::solidStVKPhysicsSolver::solveFields()
{
    pointDptr_->storePrevIter();

    volVectorField& D = Dptr_();
    pointVectorField& pointD = pointDptr_();
//    volTensorField& gradD = gradDptr_();
    volVectorField& DPrevOuterIter = DPrevOuterIterptr_();

    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    Switch explicitPredictor =
        solverControls.lookupOrDefault("explicitFirstOuterIteration",false);

    D.correctBoundaryConditions();

    DPrevOuterIter = D;

    const multiFieldIterator& mfi
    (
        db().time().lookupObject<multiFieldIterator>
        (
            "iterator"
        )
    );
    if(explicitPredictor && mfi.index() == 0)
    {
        solveExplicit();
    }
    else
    {
        solveImplicit();
    }
    // vol to point interpolation.
    {
        /*wordList wl;
        wl.clear();
        wl.setSize(D.boundaryField().size());
        forAll(wl, patchI)
        {
            if (D.boundaryField()[patchI].type() == "fixedValue")
            {
                wl[patchI] = D.boundaryField()[patchI].type();
            }
            else
            {
                wl[patchI] = "calculated";
            }
            Info << wl[patchI] << endl;
        }
        pointD = cpi_.interpolate(D, wl);*/
//        pointD = cpi_.interpolate(D, gradD);
//        pointD = cpi_.interpolate(D);
//        D.correctBoundaryConditions();
        cpi_->interpolate(D, pointD);
        pointD.correctBoundaryConditions();
    }

    // Outer loop relaxation
    if(solverControls.found("outerRelaxFactorD"))
    {
        scalar outerRelaxFactorD
        (
            readScalar(solverControls.lookup("outerRelaxFactorD"))
        );

//        D = outerRelaxFactorD*D + (1 - outerRelaxFactorD)*DPrevOuterIter;
//        gradD = fvc::grad(D);
//        sigmaptr_() = symm((mu_*(gradD + gradD.T())) 
//               + ((lambda_*I)*tr(gradD)) 
//               + (mu_*(gradD & gradD.T())) 
//               + (0.5*lambda_*I*tr(gradD & gradD.T())));
//        D.correctBoundaryConditions();

        pointDptr_() = pointDptr_->prevIter() + outerRelaxFactorD*
            (pointDptr_() - pointDptr_->prevIter());

    }
}


void Foam::solidStVKPhysicsSolver::endTimeStep()
{
    ////////////////////////////////////////
    // Nécessaire pour les restart avec FSI
    if (time().outputTime() && mesh_.moving())
    {
        //mesh_.movePoints(mesh_.points());
        mesh_.V0().write();
        mesh_.phi().write();
    }
    ////////////////////////////////////////
}

// ************************************************************************* //
