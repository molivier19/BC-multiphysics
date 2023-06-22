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

#include "solidElasticThermalPhysicsSolver.H"

#include "addToRunTimeSelectionTable.H"
#include "multiFieldIterator.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidElasticThermalPhysicsSolver, 0);

addToRunTimeSelectionTable(physicsSolver, solidElasticThermalPhysicsSolver, IOobject);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidElasticThermalPhysicsSolver::solidElasticThermalPhysicsSolver
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

    cpi_
    (
        mesh_.solutionDict().found("cellToPointInterpolation")
      ? motionInterpolation::New
        (
            mesh_, mesh_.solutionDict().lookup("cellToPointInterpolation")
        )
      : motionInterpolation::New(mesh_)
    ),

    solverControls_(mesh_.solutionDict().subDict("solverControls")),

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

    E_("E", dimPressure, mechanicalProperties_),
    nu_("nu", dimless, mechanicalProperties_),
    mu_("mu",dimPressure, Zero),
    lambda_("lambda",dimPressure, Zero),
    threeK_("threeK",dimPressure, Zero),

    alpha_("alpha", inv(dimTemperature), mechanicalProperties_),
    T0_("T0",dimTemperature, mechanicalProperties_),
    threeKalpha_("threeKalpha",dimPressure*inv(dimTemperature), Zero),

    sThermoptr_(nullptr),

    thermalStress_(mechanicalProperties_.lookupOrDefault<bool>("thermalStress", false)),
    thermalStressOnly_(mechanicalProperties_.lookupOrDefault<bool>("thermalStressOnly", false)),
    planeStress_(mechanicalProperties_.lookupOrDefault<bool>("planeStress", true)),
    compactNormalStress_(solverControls_.lookupOrDefault<bool>("compactNormalStress", true)),
    thermalCoupling_(solverControls_.lookupOrDefault<bool>("thermalCoupling", false)),

    Dptr_(nullptr),
    pointDptr_(nullptr),
    Uptr_(nullptr),
    sigmaptr_(nullptr),
    sigmaEq_(nullptr),
    gradDptr_(nullptr),
    divSigmaExpptr_(nullptr),

    DOuterConvergenceTolerance_(0.0),
    TOuterConvergenceTolerance_(0.0),
    InitialResidual_(1.0),
    InitialResidualT_(1.0),
    DRelativeTolerance_(solverControls_.get<scalar>("relD")),
    DconvergenceTolerance_(solverControls_.get<scalar>("D")),
    TRelativeTolerance_(solverControls_.get<scalar>("relT")),
    TconvergenceTolerance_(solverControls_.get<scalar>("T")),
    nCorr_(solverControls_.getOrDefault<int>("nCorrectors", 1))
{
    if (!thermalStressOnly_)
    {
        // thermal fields and properties
        Info<< "Reading mechanical properties\n" << endl;

        Info<< "Calculating Lame's coefficients\n" << endl;

        mu_ = E_/(2.0*(1.0 + nu_));
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
        threeK_ = E_/((1.0 - 2.0*nu_));

        // Correct Lamé coefficients for if plane stress
        if (planeStress_)
        {
            Info<< "Plane Stress\n" << endl;
            //- change lambda and threeK for plane stress
            lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));
            threeK_ = E_/((1.0 - nu_));
        }
    }
    else
    {
        thermalStress_ = true;
        thermalCoupling_ = false;
    }

    // Correct thermal coefficients if thermal stress
    if (thermalStress_)
    {
        Info<< "Calculating thermal coefficients\n" << endl;

        threeKalpha_ = threeK_*alpha_;
    }

    if (thermalCoupling_)
    {
        Info << "thermal coupling option activated!\n"
             <<    "Energy equation, momentum equation and\n"
             <<    "hook law are strongly coupled." << endl;
    }

    argList::addNote
    (
        "Thermo elastic solver with optional strong stress thermal coupling.\n"
        "This solver use the solid thermo model to load thermal properties.\n"
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidElasticThermalPhysicsSolver::~solidElasticThermalPhysicsSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidElasticThermalPhysicsSolver::createFields()
{
    // Create fields for solidElasticThermalPhysicsSolver

    // thermal fields and properties
    Info<< "Reading thermophysical properties and field\n" << endl;

    sThermoptr_ =
    (
        solidThermo::New(mesh_)
    );
    solidThermo& thermo = sThermoptr_();

    if (!thermo.isotropic())
    {
        FatalErrorInFunction
            << "The usage of non-isotropic\n"
            << " kappa vector field is not implemented yet."
            << abort(FatalError);
    }
    
    Info<< "Reading field D\n" << endl;
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

    Info<< "Reading field gradD\n" << endl;
    gradDptr_.reset
    (
        new volTensorField
        (
            IOobject
            (
                "gradD",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
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
        )
    );

    Info<< "Calculating stress field sigmaD\n" << endl;
    sigmaptr_.reset
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mu_*twoSymm(fvc::grad(Dptr_())) + lambda_*(I*tr(fvc::grad(Dptr_())))
        )
    );

    Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
    divSigmaExpptr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "divSigmaExp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(sigmaptr_())
        )
    );

    Info<< "Calculating VonMises Stress\n" << endl;
    sigmaEq_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "sigmaEq",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigmaptr_())))
        )
    );
    
    if (compactNormalStress_)
    {
        divSigmaExpptr_() -= fvc::laplacian(2*mu_ + lambda_, Dptr_(), "laplacian(DD,D)");
    }
    else
    {
        divSigmaExpptr_() -= fvc::div((2*mu_ + lambda_)*fvc::grad(Dptr_()), "div(sigma)");
    }

    mesh_.setFluxRequired(Dptr_().name());
    

    Info << "End Create" << endl;
}

bool Foam::solidElasticThermalPhysicsSolver::isConverged()
{
    bool DConverged = true;
    bool TConverged = true;

    if (!thermalStressOnly_)
    {
        DOuterConvergenceTolerance_ = solverControls_.get<scalar>("DOuter");
        DConverged = (InitialResidual_ < DOuterConvergenceTolerance_);
    }
    if (thermalStress_)
    {
        TOuterConvergenceTolerance_ = solverControls_.get<scalar>("TOuter");
        TConverged = (InitialResidualT_ < TOuterConvergenceTolerance_);
    }
    
    if
    (
        DConverged && TConverged
    )
    {
        if(debug)
        {
            Info << "Solid solver has converged!" << endl;
        }
        return true;
    }
    else
    {
        if(debug)
        {
            Info << "Solid solver has not converged yet..." << endl;
        }
        return false;
    }
}

void Foam::solidElasticThermalPhysicsSolver::solveFields()
{
    if (!thermalStressOnly_)
    {
        pointDptr_->storePrevIter();
        Dptr_().correctBoundaryConditions();
    }

    // ============================= Start: Solver Solid Displacement ===============
    Info << "Solving Displacement solid" << endl;

    // Declare reference fields
    solidThermo& thermo = sThermoptr_();

    // Declare mechanical properties
    volScalarField& rho = thermo.rho();

    // Declare fields
    volVectorField& D = Dptr_();
    pointVectorField& pointD = pointDptr_();
    volVectorField& U = Uptr_();
    volSymmTensorField& sigma = sigmaptr_();
    volScalarField& sigmaEq = sigmaEq_();
    volVectorField& divSigmaExp = divSigmaExpptr_();

    // Declare mechanical properties
    dimensionedScalar& mu = mu_;
    dimensionedScalar& lambda = lambda_;
    dimensionedScalar& threeKalpha = threeKalpha_;
    dimensionedScalar& T0 = T0_;

    // Declare Thermal properties
    volScalarField C("C", thermo.Cp());
    volScalarField k("k", thermo.kappa());

    // Dictionnary
    const dictionary& solverControls = solverControls_;

    // Solver control
    bool& thermalStress = thermalStress_;
    bool& compactNormalStress = compactNormalStress_;
    bool& thermalCoupling = thermalCoupling_;

    // Convergence Tolerance for the solver
    scalar& DconvergenceTolerance = DconvergenceTolerance_;
    scalar& InitialResidual = InitialResidual_;
    scalar& DRelativeTolerance = DRelativeTolerance_;
    InitialResidual = 0;
    scalar Residual = 0;
    
    // Convergence Tolerance for the solver
    scalar& TconvergenceTolerance = TconvergenceTolerance_;
    scalar& InitialResidualT = InitialResidualT_;
    scalar& TRelativeTolerance = TRelativeTolerance_;
    InitialResidualT = 0;
    scalar ResidualT = 0;

    int nCorr = nCorr_;
    int iCorr = 0;

    // These references are created so that it is
    // possible to use standard OpenFOAM include
    // files such as CourantNo.H, etc.
    // const Time& runTime = mesh_.time();
    // fvMesh& mesh = mesh_;

    // ===================================== read Solid Displacement Controls ==============================================
    nCorr = solverControls.getOrDefault<int>("nCorrectors", 1);

    if (!thermalStressOnly_)
    {
        DconvergenceTolerance = solverControls.get<scalar>("D");
        DRelativeTolerance = solverControls.get<scalar>("relD");
    }
    
    if (thermalStress_)
    {
        TconvergenceTolerance = solverControls.get<scalar>("T");
        TRelativeTolerance = solverControls.get<scalar>("relT");
    }
    
    // ===================================== read Solid Displacement Controls (end) ==============================================

    bool convergence = true;

    do
    {
        if (thermalStress)
        {
            volScalarField& T = thermo.T();

            // Solve energy equation
            if (thermalCoupling)
            {
                fvScalarMatrix TEqn
                (
                    rho*C*fvm::ddt(T)
                    ==
                    fvm::laplacian(k, T, "laplacian(DT,T)")
                    +(sigma && fvc::grad(U))
                );

                SolverPerformance<scalar> TSolPerf;
                TSolPerf = TEqn.solve();
                ResidualT = cmptMax(TSolPerf.initialResidual());
                T.relax();
            }
            else
            {
                fvScalarMatrix TEqn
                (
                    rho*C*fvm::ddt(T)
                    ==
                    fvm::laplacian(k, T, "laplacian(DT,T)")
                );

                SolverPerformance<scalar> TSolPerf;
                TSolPerf = TEqn.solve();
                ResidualT = cmptMax(TSolPerf.initialResidual());
                T.relax();
            }
        }

        if (!thermalStressOnly_)
        {
            fvVectorMatrix DEqn
            (
                rho*fvm::d2dt2(D)
             ==
                fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
                + divSigmaExp
            );

            SolverPerformance<vector> DSolPerf;
            DSolPerf = DEqn.solve();
            Residual = cmptMax(DSolPerf.initialResidual());
            D.relax();

            if (!compactNormalStress)
            {
                divSigmaExp = fvc::div(DEqn.flux());
            }
        
            volTensorField& gradD = gradDptr_();
            gradD = fvc::grad(D);
            sigma = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

            if (thermalStress && thermalCoupling)
            {
                const volScalarField& T = thermo.T();
                
                sigma = sigma - I*(threeKalpha*(T-T0));
            }

            if (compactNormalStress)
            {
                divSigmaExp = fvc::div
                (
                    sigma - (2*mu + lambda)*gradD,
                    "div(sigma)"
                );
            }
            else
            {
                divSigmaExp += fvc::div(sigma);
            }
        }

        if(iCorr == 0)
        {
            InitialResidual = Residual;
            InitialResidualT = ResidualT;
        }

        if (!thermalStressOnly_)
        {
            convergence = (Residual > DconvergenceTolerance_
                            && Residual/(InitialResidual + SMALL) > DRelativeTolerance_
                            && ++iCorr < nCorr);
        }
        else
        {
            convergence = false;
        }

        if (thermalStress_)
        {
            convergence = convergence || (ResidualT > TconvergenceTolerance_
                                            && ResidualT/(InitialResidualT + SMALL) > TRelativeTolerance_
                                            && ++iCorr < nCorr);
        }
        else
        {
            convergence = convergence || false;
        }

    } while (convergence);


    if (!thermalStressOnly_)
    {
        sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma)));
        Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;
    }
        

    Info << endl << "Solid Solver Monitoring Summary" << endl;
    if (!thermalStressOnly_)
    {
        Info << "D initial residual : " << InitialResidual << endl
            << "D final residual : " << Residual << endl;
    }
        
    if (thermalStress_)
    {
        Info << "T initial residual : " << InitialResidualT << endl
            << "T final residual : " << ResidualT << endl;
    }
    
    (iCorr < nCorr)
    ? Info << "Converged after " << ( ++iCorr ) << " iterations." << endl << endl
    : Info << "Not-converged after " << nCorr << " iterations." << endl << endl;

    // ============================= End: Solver Solid Displacement ===============
    /*if (thermalStress_)
    {
        T.correctBoundaryConditions();
    }*/

    if (!thermalStressOnly_)
    {
        U = fvc::ddt(D);

        D.correctBoundaryConditions();

        // vol to point interpolation.
        {
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

            pointD = pointDptr_->prevIter() + outerRelaxFactorD*
                (pointD - pointDptr_->prevIter());
        }
    }
}


void Foam::solidElasticThermalPhysicsSolver::startTimeStep()
{

}


void Foam::solidElasticThermalPhysicsSolver::endTimeStep()
{
    if (!thermalStressOnly_)
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
}

// ************************************************************************* //

