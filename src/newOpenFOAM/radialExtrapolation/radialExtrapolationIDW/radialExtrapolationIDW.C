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

#include "radialExtrapolationIDW.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(radialExtrapolationIDW, 0);
addToRunTimeSelectionTable
(
    radialExtrapolation,
    radialExtrapolationIDW,
    patches
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
radialExtrapolationIDW::radialExtrapolationIDW
(
    const pointField& fromPoints,
    const pointField& toPoints,
    const dictionary& dict
)
:
    radialExtrapolation(fromPoints, toPoints, dict),
    IDWdict_(subDict("IDWparameters")),
    exponent_(readScalar(IDWdict_.lookup("exponent")))
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

radialExtrapolationIDW::~radialExtrapolationIDW(){}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void radialExtrapolationIDW::generatePointWeights()
{
    Info<< "Generating weights using IDW..." << endl;


    List< pointField > fromPointsAll(Pstream::nProcs());
    fromPointsAll[Pstream::myProcNo()] = fromPoints_;
    Pstream::gatherList(fromPointsAll);
    Pstream::scatterList(fromPointsAll);


    for(label pI = 0; pI < toPoints_.size(); pI++)
    {
        labelList fromPointsProcIDsPI = fromPointsProcIDs_[pI];
        labelList fromPointsIDsPI = fromPointsIDs_[pI];
        scalarList& fromPointsWeightsPI = fromPointsWeights_[pI];

        scalar wSum(0.0);
        for(label i = 0; i < fromPointsIDsPI.size(); i++)
        {
            scalar r = mag
            (
                fromPointsAll[fromPointsProcIDsPI[i]][fromPointsIDsPI[i]]
              - toPoints_[pI]
            );

            scalar w = 1/(pow(r,exponent_));

            fromPointsWeights_[pI][i] = w;
            wSum += w;
        }
        for(label i = 0; i < fromPointsWeightsPI.size(); i++)
        {
            fromPointsWeightsPI[i] /= wSum;
        }
    }
    Info<< "done." << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
