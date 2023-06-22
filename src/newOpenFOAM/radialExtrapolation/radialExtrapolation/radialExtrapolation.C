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
#include "polyMesh.H"
#include "globalMeshData.H"
#include "pointSelection.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(radialExtrapolation, 0);

defineRunTimeSelectionTable(radialExtrapolation, patches);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
radialExtrapolation::radialExtrapolation
(
    const pointField& fromPoints,
    const pointField& toPoints,
    const dictionary& dict
)
:
    dictionary(dict.subDict("radialExtrapolation")),
    fromPointsIDs_(0),
    fromPointsProcIDs_(0),
    fromPointsWeights_(0),
    IDsNeedUpdate_(true),
    weightsNeedUpdate_(true),
    fromPointsSubsetIsUnique_(false),
    fromPoints_(fromPoints),
    toPoints_(toPoints),
    size_(fromPoints.size())
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

radialExtrapolation::~radialExtrapolation()
{
    //clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*void radialExtrapolation::initTableSizes(label nPoints)
{
    fromPointsIDs_.clear();
    fromPointsProcIDs_.clear();

    fromPointsIDs_.setSize(toPoints_.size());
    fromPointsProcIDs_.setSize(toPoints_.size());

    fromPointsWeights_.clear();
    fromPointsWeights_.setSize(toPoints().size());

    for(int n = 0; n < toPoints_.size(); n++)
    {
        fromPointsIDs_[n].clear();
        fromPointsProcIDs_[n].clear();
        fromPointsWeights_[n].clear();

        fromPointsIDs_[n].setSize(nPoints);
        fromPointsProcIDs_[n].setSize(nPoints);
        fromPointsWeights_[n].setSize(nPoints);

        for(label nn = 0; nn < nPoints; nn++)
        {
            fromPointsIDs_[n][nn] = 0;
            fromPointsProcIDs_[n][nn] = 0;
        }
    }
}*/


void radialExtrapolation::generatePointIDs()
{
    IDsNeedUpdate_ = false;

    autoPtr<pointSelection> pointSelector
    (
        pointSelection::New
        (
#if OPENFOAM >= 2112
            *this
#else
            this
#endif
        )
    );

    pointSelector->generatePointIDs
    (
        fromPointsIDs_,
        fromPointsProcIDs_,
        fromPoints_,
        toPoints_
    );

    fromPointsSubsetIsUnique_ = pointSelector->fromPointsSubsetIsUnique();
}






// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
