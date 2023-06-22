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

#include "radialExtrapolationNNI.H"
#include "addToRunTimeSelectionTable.H"
#include "pointSelectionNearestPoint.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(radialExtrapolationNNI, 0);
addToRunTimeSelectionTable
(
    radialExtrapolation,
    radialExtrapolationNNI,
    patches
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
radialExtrapolationNNI::radialExtrapolationNNI
(
    const pointField& fromPoints,
    const pointField& toPoints,
    const dictionary& dict
)
:
    radialExtrapolation(fromPoints, toPoints, dict)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

radialExtrapolationNNI::~radialExtrapolationNNI(){}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void radialExtrapolationNNI::generatePointWeights()
{
    forAll(fromPointsWeights_, pI)
    {
        fromPointsWeights_[pI][0] = 1.0;
    }
}


void radialExtrapolationNNI::generatePointIDs()
{
    IDsNeedUpdate_ = false;

#if OPENFOAM >= 2112
    pointSelectionNearestPoint pointSelector(*this);
#else
    pointSelectionNearestPoint pointSelector(this);
#endif

    pointSelector.generatePointIDs
    (
        fromPointsIDs_,
        fromPointsProcIDs_,
        fromPoints_,
        toPoints_
    );
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
