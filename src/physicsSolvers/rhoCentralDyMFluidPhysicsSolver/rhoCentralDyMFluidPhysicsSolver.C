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

#include "rhoCentralDyMFluidPhysicsSolver.H"

#include "addToRunTimeSelectionTable.H"
#include "pointFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rhoCentralDyMFluidPhysicsSolver, 0);

addToRunTimeSelectionTable(physicsSolver, rhoCentralDyMFluidPhysicsSolver, IOobject);

// * * * * * * * * * * * * Protected member functions  * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::rhoCentralDyMFluidPhysicsSolver::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& dir,
    const word& reconFieldName
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        fvc::interpolate
        (
            vf,
            dir,
            "reconstruct("
          + (reconFieldName != word::null ? reconFieldName : vf.name())
          + ')'
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    sf.rename(vf.name() + '_' + dir.name());

    return tsf;
}

void Foam::rhoCentralDyMFluidPhysicsSolver::computeExplicitFluxes()
{
    // References
    volVectorField& U = Uptr_();
    volScalarField& rho = rhoptr_();
    volVectorField& rhoU = rhoUptr_();
    surfaceScalarField& pos = posptr_();
    surfaceScalarField& neg = negptr_();
    surfaceScalarField& phi = phiptr_();
    psiThermo& thermo = pThermoptr_();
    volScalarField& e = thermo.he();
    volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();
    const Time& runTime = mesh_.time();
    dynamicFvMesh& mesh = mesh_;
    // scalar& CoNum = CoNum_;
    // scalar& meanCoNum = meanCoNum_;
    const dimensionedScalar& v_zero = v_zero_;
    word& fluxScheme = fluxScheme_;

    // --- Directed interpolation of primitive fields onto faces

    surfaceScalarField rho_pos(interpolate(rho, pos));
    surfaceScalarField rho_neg(interpolate(rho, neg));

    surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
    surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

    volScalarField rPsi("rPsi", 1.0/psi);
    surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
    surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

    surfaceScalarField e_pos(interpolate(e, pos, T.name()));
    surfaceScalarField e_neg(interpolate(e, neg, T.name()));

    surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
    surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

    surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
    surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

    surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
    surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

    // Make fluxes relative to mesh-motion
    if (mesh.moving())
    {
        phiv_pos -= mesh.phi();
        phiv_neg -= mesh.phi();
    }
    // Note: extracted out the orientation so becomes unoriented
    phiv_pos.setOriented(false);
    phiv_neg.setOriented(false);

    volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));

    surfaceScalarField cSf_pos
    (
        "cSf_pos",
        interpolate(c, pos, T.name())*mesh.magSf()
    );
    surfaceScalarField cSf_neg
    (
        "cSf_neg",
        interpolate(c, neg, T.name())*mesh.magSf()
    );

    surfaceScalarField ap
    (
        "ap",
        max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
    );
    surfaceScalarField am
    (
        "am",
        min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
    );

    surfaceScalarField a_pos("a_pos", ap/(ap - am));

    surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

    surfaceScalarField aSf("aSf", am*a_pos);

    if (fluxScheme == "Tadmor")
    {
        aSf = -0.5*amaxSf;
        a_pos = 0.5;
    }

    surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

    phiv_pos *= a_pos;
    phiv_neg *= a_neg;

    surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
    surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

    // Reuse amaxSf for the maximum positive and negative fluxes
    // estimated by the central scheme
    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

    // Courant number
    if (mesh.nInternalFaces())
    {
        scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

        CoNum_ = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();

        meanCoNum_ =
            0.5*(gSum(sumAmaxSf)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }

    Info<< "Mean and max Courant Numbers = "
        << meanCoNum_ << " " << CoNum_ << endl;


    phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

    surfaceVectorField phiU(aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg);

    // Note: reassembled orientation from the pos and neg parts so becomes
    // oriented
    phiU.setOriented(true);

    surfaceVectorField phiUp(phiU + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf());

    surfaceScalarField phiEp
    (
        "phiEp",
        aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
      + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
      + aSf*p_pos - aSf*p_neg
    );

    // Make flux for pressure-work absolute
    if (mesh.moving())
    {
        surfaceScalarField phia(a_pos*p_pos + a_neg*p_neg);
        phia.setOriented(true);

        phiEp += mesh.phi()*phia;
    }

    volScalarField muEff("muEff", turbulence_->muEff());
    volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

    surfaceScalarField sigmaDotU
    (
        "sigmaDotU",
        (
            fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
          + fvc::dotInterpolate(mesh.Sf(), tauMC)
        )
      & (a_pos*U_pos + a_neg*U_neg)
    );


    divPhi_() = fvc::div(phi);
    divPhiUp_() = fvc::div(phiUp);
    divEnergy_() = fvc::div(phiEp) - fvc::div(sigmaDotU);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoCentralDyMFluidPhysicsSolver::rhoCentralDyMFluidPhysicsSolver
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

    pThermoptr_(nullptr),
    inviscid_(true),

    CoNum_(0.0),
    meanCoNum_(0.0),

    v_zero_(dimVolume/dimTime, Zero),

    Uptr_(nullptr),
    rhoptr_(nullptr),
    rhoUptr_(nullptr),
    rhoEptr_(nullptr),

    posptr_(nullptr),
    negptr_(nullptr),
    phiptr_(nullptr),

    divPhi_(nullptr),
    divPhiUp_(nullptr),
    divEnergy_(nullptr),

    sigmaFluidptr_(nullptr),

    turbulence_(nullptr),

    // cumulativeContErr_(0),

    UInitialResidual_(0.0),
    hInitialResidual_(0.0),

    UConvergenceTolerance_(0.0),
    URelativeTolerance_(0.0),
    hConvergenceTolerance_(0.0),
    hRelativeTolerance_(0.0),

    fluxScheme_("Kurganov")
{
    argList::addNote
    (
        "Density-based compressible flow solver based on central-upwind"
        " schemes of Kurganov and Tadmor.\n"
        "With support for mesh-motion and topology changes."
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoCentralDyMFluidPhysicsSolver::~rhoCentralDyMFluidPhysicsSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rhoCentralDyMFluidPhysicsSolver::createFields()
{

    // Create fields
    Info<< "Reading thermophysical properties\n" << endl;

    pThermoptr_ =
    (
        psiThermo::New(mesh_)
    );

    // Create reference field
    psiThermo& thermo = pThermoptr_();
    volScalarField& e = thermo.he();
    const volScalarField& mu = thermo.mu();

    // Info<< "For info heat capacity ratio (gamma) value : " << pThermoptr_().Cp()/pThermoptr_().Cv() << " \n" << endl;

    Info<< "Reading field U\n" << endl;
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

    rhoptr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermo.rho()
        )
    );

    rhoUptr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "rhoU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Uptr_()*rhoptr_()
        )
    );

    rhoEptr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "rhoE",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rhoptr_()*(e + 0.5*magSqr(Uptr_()))
        )
    );

    Info<< "Reading field pos\n" << endl;
    posptr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "pos",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("pos", dimless, 1.0)
        )
    );


    negptr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "neg",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("neg", dimless, -1.0)
        )
    );
    
    phiptr_.reset
    (
        new surfaceScalarField
        (
            "phi",
            fvc::flux(rhoUptr_())
        )
    );


    divPhi_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "divPhi",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("", dimDensity/dimTime, 0.0)
        )
    );

    divPhiUp_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "divPhiUp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("", dimDensity*dimVelocity/dimTime, vector(0,0,0))
        )
    );

    divEnergy_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "divEnergy",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("", dimEnergy/(dimVolume*dimTime), 0.0)
        )
    );


    Info<< "Creating turbulence model\n" << endl;
    turbulence_ =
    (
        compressible::turbulenceModel::New
        (
            rhoptr_(),
            Uptr_(),
            phiptr_(),
            thermo
        )
    );

    Info<< "Creating Sigma Fluid Field\n" << endl;
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
          - thermo.rho()*(
          thermo.p()*symmTensor(1,0,0,1,0,1)
          + turbulence_->devRhoReff())

        )
    );

    if (max(mu.primitiveField()) > 0.0)
    {
        inviscid_ = false;
    }

    // Validate for turbulence
    turbulence_ -> validate();

    // Read flux scheme
    if (mesh_.schemesDict().readIfPresent("fluxScheme", fluxScheme_))
    {
        if ((fluxScheme_ == "Tadmor") || (fluxScheme_ == "Kurganov"))
        {
            Info<< "\nfluxScheme: " << fluxScheme_ << endl;
        }
        else
        {
            FatalErrorInFunction
                << "fluxScheme: " << fluxScheme_
                << " is not a valid choice. "
                << "Options are: Tadmor, Kurganov"
                << abort(FatalError);
        }
    }

    computeExplicitFluxes();
}


bool Foam::rhoCentralDyMFluidPhysicsSolver::isConverged()
{
    // TODO: check...
    // Central is an explicit solver.
    // There is no convergence criteria.
    return true;
}


void Foam::rhoCentralDyMFluidPhysicsSolver::solveFields()
{
    // Reference
    volVectorField& U = Uptr_();
    volScalarField& rho = rhoptr_();
    volVectorField& rhoU = rhoUptr_();
    volScalarField& rhoE = rhoEptr_();
    volSymmTensorField& sigmaFluid = sigmaFluidptr_();
    psiThermo& thermo = pThermoptr_();
    volScalarField& e = thermo.he();
    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    dynamicFvMesh& mesh = mesh_;
    bool& inviscid = inviscid_;

    // Do any mesh changes
    mesh.update();


    // --- Solve density
    solve(fvm::ddt(rho) + divPhi_());

    // --- Solve momentum
    solve(fvm::ddt(rhoU) + divPhiUp_());

    U.ref() =
        rhoU()
       /rho();
    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

    if (!inviscid)
    {
        volScalarField muEff("muEff", turbulence_->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        solve
        (
            fvm::ddt(rho, U) - fvc::ddt(rho, U)
          - fvm::laplacian(muEff, U)
          - fvc::div(tauMC)
        );
        rhoU = rho*U;
    }

    // --- Solve energy
    solve(fvm::ddt(rhoE) + divEnergy_());

    e = rhoE/rho - 0.5*magSqr(U);
    e.correctBoundaryConditions();
    thermo.correct();
    rhoE.boundaryFieldRef() ==
        rho.boundaryField()*
        (
            e.boundaryField() + 0.5*magSqr(U.boundaryField())
        );

    if (!inviscid)
    {
        solve
        (
            fvm::ddt(rho, e) - fvc::ddt(rho, e)
          - fvm::laplacian(turbulence_->alphaEff(), e)
        );
        thermo.correct();
        rhoE = rho*(e + 0.5*magSqr(U));
    }

    p.ref() =
        rho()
       /psi();
    p.correctBoundaryConditions();
    rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

    turbulence_->correct();


    // TODO: check
    sigmaFluid = -rho*
        (
            p*symmTensor(1,0,0,1,0,1)
          + turbulence_->devRhoReff()
        );
}


void Foam::rhoCentralDyMFluidPhysicsSolver::startTimeStep()
{
    // Info<< "Start Time Step\n" << endl;
}

void Foam::rhoCentralDyMFluidPhysicsSolver::endTimeStep()
{
    // Info<< "End Time Step\n" << endl;

    computeExplicitFluxes();
}

bool Foam::rhoCentralDyMFluidPhysicsSolver::computeDeltaT()
{
    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    bool computeDeltaT =
        solverControls.lookupOrDefault<bool>("computeDeltaT", false);

    return computeDeltaT;
}

scalar Foam::rhoCentralDyMFluidPhysicsSolver::returnDeltaT()
{
    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    scalar maxCo(readScalar(solverControls.lookup("maxCo")));

    return time().deltaT().value()*maxCo/(CoNum_ + SMALL);
}

// ************************************************************************* //
