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

#include "physicsSolverTemplate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(physicsSolverTemplate, 0);

addToRunTimeSelectionTable(physicsSolver, physicsSolverTemplate, IOobject);

// * * * * * * * * * * * * * * * Custom Functions  * * * * * * * * * * * * * //

void Foam::physicsSolverTemplate::computeExplicitFluxes()
{
  // Field references
  volVectorField& U = Uptr_();
  volScalarField& rho = rhoptr_();

  // Compute value
  volVectorField phiU = ...
  volScalarField phiRho = ...

  // Assign value
  divPhiU_() = fvc::div(phiU);
  divPhiRho_() = fvc::div(phiRho);
}

bool Foam::physicsSolverTemplate::validateIODictionary()
{
  // Reference
  IOdictionary& properties = properties_;

  // Local variable
  bool validData = false;

  // Validate Data
  if(...)
  {
    validData = true;
  }

  return validData;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicsSolverTemplate::physicsSolverTemplate
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

    properties_
    (
        IOobject
        (
            "properties",
            io.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    Uptr_(nullptr),
    rhoptr_(nullptr),

    divPhiU_(nullptr),
    divPhiRho_(nullptr)
{
    // =======================================
    //
    // Actions to do upon the construction 
    // of the solver object
    //
    // =======================================
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::physicsSolverTemplate::~physicsSolverTemplate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::physicsSolverTemplate::createFields()
{
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

    divPhiU_.reset
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

    divPhiRho_.reset
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
}


bool Foam::physicsSolverTemplate::isConverged()
{
    bool converged = false;

    // Check convergence
    if (...)
    {
        converged = true;
    }

    return false;   
}


void Foam::physicsSolverTemplate::solveFields()
{
    // Mesh and fields references
    dynamicFvMesh& mesh = mesh_;
    volVectorField& U = Uptr_();
    volScalarField& rho = rhoptr_();
    volVectorField& divPhiU = divPhiUptr_();
    volScalarField& divPhiRho = divPhiRhoptr_();

    // =======================================
    //
    //          Solve physical fields
    //
    // =======================================
}


void Foam::physicsSolverTemplate::startTimeStep()
{
    // =======================================
    //
    // Actions to do before solving the fields
    //
    // =======================================
}

void Foam::physicsSolverTemplate::endTimeStep()
{
    // =======================================
    //
    // Actions to do after solving the fields
    //
    // =======================================
}

bool Foam::physicsSolverTemplate::computeDeltaT()
{
    bool computeDeltaT = false;

    // Check if adaptative time step is required
    if (...)
    {
        computeDeltaT = true;
    }

    return computeDeltaT;
}


scalar Foam::physicsSolverTemplate::returnDeltaT()
{
    // Compute adaptative time step

    return deltaTime;
}
// ************************************************************************* //
