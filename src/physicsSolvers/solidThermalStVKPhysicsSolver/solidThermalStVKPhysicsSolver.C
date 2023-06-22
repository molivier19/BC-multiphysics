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


#include "solidThermalStVKPhysicsSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "steadyStateDdtScheme.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "multiFieldIterator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidThermalStVKPhysicsSolver, 0);

addToRunTimeSelectionTable(physicsSolver, solidThermalStVKPhysicsSolver, IOobject);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solidThermalStVKPhysicsSolver::solveImplicit()
{
    // Declare reference fields
    solidThermo& thermo = sThermoptr_();

    // Declare mechanical properties
    volScalarField& rho = thermo.rho();

    // Declare Thermal properties
    volScalarField C("C", thermo.Cp());
    volScalarField k("k", thermo.kappa());

    volVectorField& D = Dptr_();
    volVectorField& U = Uptr_();
    volTensorField& gradD = gradDptr_();
    volSymmTensorField& sigma = sigmaptr_();
    volSymmTensorField& sigmaCauchy = sigmaCauchyptr_();
    volScalarField& sigmaEq = sigmaEqptr_();

    bool& thermalStress = thermalStress_;
    bool& thermalCoupling = thermalCoupling_;

    // Solver controls
    int nCorr = readInt(solverControls_.lookup("nCorrectors"));
    int corr = 0;

    if (!thermalStressOnly_)
    {
        DConvergenceTolerance_ = readScalar(solverControls_.lookup("D"));
        DRelativeTolerance_ = readScalar(solverControls_.lookup("relD"));
    }

    if (thermalStress_)
    {
        TConvergenceTolerance_ = readScalar(solverControls_.lookup("T"));
        TRelativeTolerance_ = readScalar(solverControls_.lookup("relT"));
        thermalCoupling = solverControls_.get<bool>("thermalCoupling");
    }

    scalar DResidual = VGREAT;
    DInitialResidual_ = 0;
    scalar TResidual = VGREAT;
    TInitialResidual_ = 0;

#   include "solidThermalStVKDdtCoeffs.H"

    bool convergence = true;

    do
    {
        // Solve momentum equation
        D.storePrevIter();

        volTensorField F = (I + gradD);
        volScalarField J = det(F);
        volTensorField PK1Tensor = sigma & F;        
        
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
                    +(PK1Tensor && fvc::ddt(F))
                );

                SolverPerformance<scalar> TSolPerf;
                TSolPerf = TEqn.solve();
                TResidual = cmptMax(TSolPerf.initialResidual());
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
                TResidual = cmptMax(TSolPerf.initialResidual());
                T.relax();
            }
        }

        if (!thermalStressOnly_)
        {
            fvVectorMatrix DEqn
            (
                rho*
                (
                    Cn*fvm::ddt(D) 
                    - Co*U.oldTime() 
                    + Coo*U.oldTime().oldTime()
                )
            ==
                fvm::laplacian(2*mu_ + lambda_, D, "laplacian(DD,D)") 
            - fvc::laplacian(2*mu_ + lambda_, D, "laplacianExp(DD,D)") 
            + fvc::div(PK1Tensor)
            );

            SolverPerformance<vector> DSolPerf;
            DSolPerf = DEqn.solve();
            DResidual = cmptMax(DSolPerf.initialResidual());
            D.relax();
          

            gradD = fvc::grad(D);
            sigma = symm((mu_*(gradD + gradD.T())) 
                + ((lambda_*I)*tr(gradD)) 
                + (mu_*(gradD & gradD.T())) 
                + (0.5*lambda_*I*tr(gradD & gradD.T()))
                );

            if (thermalStress && thermalCoupling)
            {
                const volScalarField& T = thermo.T();
                sigma = sigma - I*(threeKalpha_*(T-T0_));
            }
        }

        if(corr == 0)
        {
            DInitialResidual_ = DResidual;
            TInitialResidual_ = TResidual;
        }

        if (!thermalStressOnly_)
        {
            convergence = (DResidual > DConvergenceTolerance_
                            && DResidual/(DInitialResidual_ + SMALL) > DRelativeTolerance_
                            && ++corr < nCorr);
        }
        else
        {
            convergence = false;
        }

        if (thermalStress_)
        {
            convergence = convergence || (TResidual > TConvergenceTolerance_
                                            && TResidual/(TInitialResidual_ + SMALL) > TRelativeTolerance_
                                            && ++corr < nCorr);
        }
        else
        {
            convergence = convergence || false;
        }

    } while (convergence);

    if (!thermalStressOnly_)
    {
        volTensorField F = (I + gradD);
        volScalarField J = det(F);
        sigmaCauchy = symm(J*(F & sigma & F.T()));
        sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigmaCauchy)));
        Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;
    }

    Info << endl << "Solid Solver Monitoring Summary" << endl;
    if (!thermalStressOnly_)
    {
        Info << "D initial residual : " << DInitialResidual_ << endl
            << "D final residual : " << DResidual << endl;
    }
        
    if (thermalStress_)
    {
        Info << "T initial residual : " << TInitialResidual_ << endl
            << "T final residual : " << TResidual << endl;
    }
    
    (corr < nCorr)
    ? Info << "Converged after " << ( ++corr ) << " iterations." << endl << endl
    : Info << "Not-converged after " << nCorr << " iterations." << endl << endl;

    U = fvc::ddt(D);    
}

void Foam::solidThermalStVKPhysicsSolver::solveExplicit()
{
    // Declare reference fields
    solidThermo& thermo = sThermoptr_();

    // Declare mechanical properties
    volScalarField& rho = thermo.rho();

    // Declare Thermal properties
    volScalarField C("C", thermo.Cp());
    volScalarField k("k", thermo.kappa());

    volVectorField& D = Dptr_();
    //pointVectorField& pointD = pointDptr_();
    volVectorField& U = Uptr_();
    volTensorField& gradD = gradDptr_();
    volSymmTensorField& sigma = sigmaptr_();
    volSymmTensorField& sigmaCauchy = sigmaCauchyptr_();
    volScalarField& sigmaEq = sigmaEqptr_();
    //volVectorField& DPrevOuterIter = DPrevOuterIterptr_();

    bool& thermalStress = thermalStress_;
    bool& thermalCoupling = thermalCoupling_;

    // Required to ajust the load coming from fluid in FSI.
    //D.correctBoundaryConditions();
    gradD = fvc::grad(D);
    sigma = symm((mu_*(gradD + gradD.T())) 
           + ((lambda_*I)*tr(gradD)) 
           + (mu_*(gradD & gradD.T())) 
           + (0.5*lambda_*I*tr(gradD & gradD.T())));

    if (thermalStress && thermalCoupling)
    {
        const volScalarField& T = thermo.T();
        sigma = sigma - I*(threeKalpha_*(T-T0_));
    }

    dimensionedScalar Dn("Cn", dimless/(dimTime*dimTime), 0.0);
    dimensionedScalar Do("Co", dimless/(dimTime*dimTime), 0.0);
    dimensionedScalar Doo("Coo", dimless/(dimTime*dimTime), 0.0);
    dimensionedScalar Uo("Co", dimless/dimTime, 0.0);

    dimensionedScalar Tn("Tn", dimless/(dimTime), 0.0);
    dimensionedScalar To("To", dimless/(dimTime), 0.0);
    dimensionedScalar Too("Too", dimless/(dimTime), 0.0);

    if(mesh_.time().timeIndex() == 1)
    {
        // Explicit central (2nd order)
        // Start-up procedure
        scalar deltaT = mesh_.time().deltaT().value();
        Dn.value() = 2.0/(deltaT*deltaT);
        Do.value() = -2.0/(deltaT*deltaT);
        Uo.value() = -2.0/(deltaT);

        // Explicit forward difference (1st order)
        Tn.value() = 1/(2*deltaT);
        To.value() = -1/(2*deltaT);
    }
    else
    {
        // Explicit central (2nd order)
        // Note: does not support variable timesteps
        scalar deltaT = mesh_.time().deltaT().value();
        Dn.value() = 1.0/(deltaT*deltaT);
        Do.value() = -2.0/(deltaT*deltaT);
        Doo.value() = 1.0/(deltaT*deltaT);

        // Explicit central (2nd order)
        Tn.value() = 1/(deltaT);
        Too.value() = -1/(deltaT);
    }

    volTensorField F = (I + gradD);
    volScalarField J = det(F);
    volTensorField PK1Tensor = sigma & F;        
    
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
                +(PK1Tensor && fvc::ddt(F))
            );
            TEqn.solve();
        }
        else
        {
            fvScalarMatrix TEqn
            (
                rho*C*fvm::ddt(T)
                ==
                fvm::laplacian(k, T, "laplacian(DT,T)")
            );
            TEqn.solve();
        }
        
    }

    if (!thermalStressOnly_)
    {
        fvVectorMatrix DEqn
        (
            rho*
            (
                fvm::Sp(Dn,D)
                + Do*D.oldTime()
                + Doo*D.oldTime().oldTime()
                + Uo*U.oldTime()
            )
        - fvc::div(PK1Tensor)
        );
        DEqn.solve();

        U = fvc::ddt(D); 
        gradD = fvc::grad(D);
        sigma = symm((mu_*(gradD + gradD.T())) 
            + ((lambda_*I)*tr(gradD)) 
            + (mu_*(gradD & gradD.T())) 
            + (0.5*lambda_*I*tr(gradD & gradD.T())));

        if (thermalStress && thermalCoupling)
        {
            const volScalarField& T = thermo.T();
            sigma = sigma - I*(threeKalpha_*(T-T0_));
        }

        sigmaCauchy = symm(J*(F & sigma & F.T()));
        sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigmaCauchy)));
        Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermalStVKPhysicsSolver::solidThermalStVKPhysicsSolver
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
    thermalCoupling_(solverControls_.lookupOrDefault<bool>("thermalCoupling", false)),

    Dptr_(nullptr),
    pointDptr_(nullptr),
    Uptr_(nullptr),
    gradDptr_(nullptr),
    sigmaptr_(nullptr),
    sigmaCauchyptr_(nullptr),
    sigmaEqptr_(nullptr),
    DPrevOuterIterptr_(nullptr),

    DInitialResidual_(1.0),
    TInitialResidual_(1.0),

    DConvergenceTolerance_(0.0),
    DRelativeTolerance_(0.0),
    TConvergenceTolerance_(0.0),
    TRelativeTolerance_(0.0),

    DOuterConvergenceTolerance_(0.0),
    TOuterConvergenceTolerance_(0.0)
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
        "The thermal stress is integrated according to the large deformation\n formulation and the St-Venant-Kirchhoff model.\n"
        "This solver use the solid thermo model to load thermal properties.\n"
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidThermalStVKPhysicsSolver::~solidThermalStVKPhysicsSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidThermalStVKPhysicsSolver::createFields()
{
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
            << " kappa vector field is not yet implemented."
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

    Info<< "Calculating stress field Second Piola-Kirchhoff tensor\n" << endl;
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

    sigmaCauchyptr_.reset
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigmaCauchy",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            symm(det(I + gradDptr_()) * ((I + gradDptr_()) 
                    & sigmaptr_() & (I + gradDptr_().T())))
            
        )
    );

    Info<< "Calculating VonMises Stress\n" << endl;
    sigmaEqptr_.reset
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
            sqrt((3.0/2.0)*magSqr(dev(sigmaCauchyptr_())))
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

    //int nCorr = readInt(solverControls.lookup("nCorrectors"));
    DConvergenceTolerance_ = readScalar(solverControls_.lookup("D"));
    DRelativeTolerance_ = readScalar(solverControls_.lookup("relD"));
    //DOuterConvergenceTolerance_ = readScalar(solverControls.lookup("DOuter"));
}


bool Foam::solidThermalStVKPhysicsSolver::isConverged()
{
    bool DConverged = true;
    bool TConverged = true;

    if (!thermalStressOnly_)
    {
        DOuterConvergenceTolerance_ = solverControls_.get<scalar>("DOuter");
        DConverged = (DInitialResidual_ < DOuterConvergenceTolerance_);
    }
    if (thermalStress_)
    {
        TOuterConvergenceTolerance_ = solverControls_.get<scalar>("TOuter");
        TConverged = (TInitialResidual_ < TOuterConvergenceTolerance_);
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

void Foam::solidThermalStVKPhysicsSolver::solveFields()
{
    pointDptr_->storePrevIter();

    //solidThermo& thermo = sThermoptr_();
    //volScalarField& T = thermo.T();
    volVectorField& D = Dptr_();
    pointVectorField& pointD = pointDptr_();
//    volTensorField& gradD = gradDptr_();
    volVectorField& DPrevOuterIter = DPrevOuterIterptr_();

    const dictionary& solverControls = solverControls_;

    Switch explicitPredictor =
        solverControls.lookupOrDefault("explicitFirstOuterIteration",false);

    D.correctBoundaryConditions();
    //T.correctBoundaryConditions();

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

        pointDptr_() = pointDptr_->prevIter() + outerRelaxFactorD*
            (pointDptr_() - pointDptr_->prevIter());
    }
}


void Foam::solidThermalStVKPhysicsSolver::endTimeStep()
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
