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

#include "pseudotimeControl.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudotimeControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

#if OPENFOAM >= 1906
bool Foam::pseudotimeControl::read()
#else
void Foam::pseudotimeControl::read()
#endif
{
    solutionControl::read(false);


    // Read solution controls
    const dictionary& pseudoDict = dict();
    nCorrOuter_ = readInt(fieldsDictionary_.lookup("nOuterCorr")); // nCorrOuter will be adjusted by the multiField solver
    nCorrOuterMin_ = fieldsDictionary_.lookupOrDefault("nOuterCorrMin", 1); // nCorrOuterMin will be adjusted by the multiField solver
    if (pseudoDict.found("pseudoTol"))
    {
        pseudoDict.lookup("pseudoTol") >> residualTols_;
    }
    if (pseudoDict.found("pseudoTolRel"))
    {
        pseudoDict.lookup("pseudoTolRel") >> residualTolsRel_;
    }
    turbOnFinalIterOnly_ =
        pseudoDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", false);

    if (state_.found("initResiduals"))
    {
        state_.lookup("initResiduals") >> initResiduals_;
    }
    else
    {
        initResiduals_ = residualIO(residualTols_, 0.0);
    }
#if OPENFOAM >= 1906
    return true;
#endif
}


bool Foam::pseudotimeControl::criteriaSatisfied()
{
    // no checks on first iteration after restart or first iteration of 
    // timestep - nothing has been calculated yet
    /*if (firstIteration_ || corr_ == 1 || finalIter())
    {
        firstIteration_ = false;
        return false;
    }*/

    // no check if final interation
    if (finalIter())
    {
        firstIteration_ = false;
        return false;
    }

    bool storeIni = this->storeInitialResiduals();

    if (storeIni)
    {
        initResiduals_ = residuals_;
        state_.add("initResiduals", initResiduals_, true); // Store for restarts
    }

    bool absCheck = (max(residuals_/residualTols_) < 1.0);

    bool relCheck = false;

    if (!storeIni)
    {
        relCheck = (max(residuals_/(initResiduals_+ROOTVSMALL)/(residualTolsRel_+ROOTVSMALL)) < 1.0);
    }

    return (corr_ > nCorrOuterMin_) && (absCheck || relCheck);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pseudotimeControl::pseudotimeControl(fvMesh& mesh, const bool steadyState, const label nScalars, const label nVectors, const residualIO& defaultTol, const residualIO& defaultTolRel)
:
    solutionControl(mesh, "pseudoTime"),
    steadyState_(steadyState),
    nCorrOuter_(0),
    nCorrOuterMin_(0),
    fieldsDictionary_
    (
        IOobject
        (
            "fieldsDictionary",
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    turbOnFinalIterOnly_(false),
    converged_(false),
    firstIteration_(true),
    residualTols_(defaultTol),
    residualTolsRel_(defaultTolRel),
    residuals_(defaultTol, 0.0),
    initResiduals_(defaultTol, 0.0),
    state_
    (
        IOobject
        (
            "pseudotimeState",
            mesh.time().timeName(),
            "uniform",
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    )
{
    read();

    Info<< algorithmName_ << ": ";
    if (!steadyState)
    {
        Info<< "max outer-iterations = " << nCorrOuter_ << ", ";
    }
    Info<< "tolerance = " << residualTols_ << ", "
        << "relTol = " << residualTolsRel_
        << nl;
    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudotimeControl::~pseudotimeControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pseudotimeControl::loop()
{
    read();

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    // Loop detection of the convergence.
    bool completed = false;
    if (!steadyState_ && corr_ == nCorrOuter_-1)
    {
        if (!steadyState_)
        {
                if (debug)
                {   
                    Info<< algorithmName_ << ": Remove Final Iteration Keyword. " << endl;
                }
                mesh_.data::remove("finalIteration");
        }

        if (nCorrOuter_ != 0)
        {
            Info<< algorithmName_ << ": not converged within "
                << corr_ << " outer-iterations" << endl << endl;
        }
    }
    else
    {
        // Criteria satisfied
        if (criteriaSatisfied() && !converged_) 
        {
            converged_ = true;
        }
        // Converged and all criteria are satisfied
        else if (converged_) 
        {   
            Info<< algorithmName_ << ": converged for outer-iteration : " << corr_ << endl << endl;

            storePrevIterFields();

            completed = true;
        }

        // Thing to do every iteration not completed
        if (!completed)
        {
            // Thing to do for final iteration
            if (finalIter())
            {
                if (!steadyState_)
                {
                    if (debug)
                    {   
                        Info<< algorithmName_ << ": Add Final Iteration Keyword. " << endl;
                    }
                    mesh_.data::add("finalIteration", true);
                }
            }

            if (nCorrOuter_ != 0)
            {
                Info<< algorithmName_ << ": outer-iteration " << corr_ << endl << endl;

                storePrevIterFields();
            }
        }
    }

    return completed;
}


// ************************************************************************* //
