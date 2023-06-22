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

#include "radialExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<radialExtrapolation> radialExtrapolation::New
(
    const pointField& fromPoints,
    const pointField& toPoints,
    const dictionary& dict
)
{
    word type(dict.subDict("radialExtrapolation").lookup("method"));

#if OPENFOAM >= 2112
    auto cstrIter = patchesConstructorTablePtr_->find(type);
#else
    patchesConstructorTable::iterator cstrIter =
        patchesConstructorTablePtr_->find(type);
#endif

    if (cstrIter == patchesConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "radialExtrapolation::New\n\
    (\n\
        const pointField& fromPoints,\n\
        const pointField& toPoints,\n\
        const dictionary& dict\n\
    )"
        )   << "Unknown radialExtrapolation type " << type
            << endl << endl
            << "Valid radialExtrapolation types are :" << endl
            << patchesConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<radialExtrapolation>
    (
        cstrIter()(fromPoints, toPoints, dict)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
