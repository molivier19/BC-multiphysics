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
    pseudotimeControl

Description
    Adaptation of the pseudo-time control library for the "multiFieldSolver"
    application. Adapted and build with respect to the library HiSA version
    "hisa 1.4.8" for OpenFOAM-v2006. The source of the HISA solver can be
    found here:

    https://hisa.gitlab.io
    https://gitlab.com/hisa/hisa

\*---------------------------------------------------------------------------*/ 

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::pseudotimeControl::nCorrOuter() const
{
    return nCorrOuter_;
}


inline bool Foam::pseudotimeControl::storeInitialResiduals() const
{
    return firstIter();
}


inline bool Foam::pseudotimeControl::firstIter() const
{
    return corr_ == 0;
}


inline bool Foam::pseudotimeControl::finalIter() const
{
    return converged_ || (corr_ == nCorrOuter_-2);
}


inline bool Foam::pseudotimeControl::turbCorr() const
{
    return !turbOnFinalIterOnly_ || finalIter();
}


inline void Foam::pseudotimeControl::setResidual(const residualIO& res)
{
    residuals_ = res;
}

inline void Foam::pseudotimeControl::init()
{
    if (!steadyState_)
    {
        if (debug)
        {   
            Info<< algorithmName_ << ": Remove Final Iteration Keyword. " << endl;
        }
        mesh_.data::remove("finalIteration");
    }

    converged_ = false;
}


// ************************************************************************* //
