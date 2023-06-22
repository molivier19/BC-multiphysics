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

#include "pointSelectionNearestGroup.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointSelectionNearestGroup, 0);
addToRunTimeSelectionTable
(
    pointSelection,
    pointSelectionNearestGroup,
    patches
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pointSelectionNearestGroup::pointSelectionNearestGroup
(
    const dictionary& dict
)
:
    pointSelection(dict)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

pointSelectionNearestGroup::~pointSelectionNearestGroup(){}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointSelectionNearestGroup::generatePointIDs
(
    labelListList& fromPointsIDs,
    labelListList& fromPointsProcIDs,
    const pointField& fromPoints,
    const pointField& toPoints
)
{
    scalar nClosestPointsUsed =
        readLabel(pointSelectionDict().lookup("nClosestPointsUsed"));

    Info<< "Finding " << nClosestPointsUsed 
        << " nearest points..." <<endl;


    fromPointsIDs.clear();
    fromPointsProcIDs.clear();

    fromPointsIDs.setSize(toPoints.size());
    fromPointsProcIDs.setSize(toPoints.size());

    for(int n = 0; n < toPoints.size(); n++)
    {
        fromPointsIDs[n].clear();
        fromPointsProcIDs[n].clear();

        fromPointsIDs[n].setSize(nClosestPointsUsed);
        fromPointsProcIDs[n].setSize(nClosestPointsUsed);

        for(label nn = 0; nn < nClosestPointsUsed; nn++)
        {
            fromPointsIDs[n][nn] = 0;
            fromPointsProcIDs[n][nn] = 0;
        }
    }



    List< pointField > fromPointsAll(Pstream::nProcs());
    fromPointsAll[Pstream::myProcNo()] = fromPoints;
    Pstream::gatherList(fromPointsAll);
    Pstream::scatterList(fromPointsAll);

    for(int toI = 0; toI < toPoints.size(); toI++)
    {
        List<scalar> nMinDist;
        nMinDist.clear();
        nMinDist.setSize(nClosestPointsUsed);
        List<vector> pointPosition;
        pointPosition.clear();
        pointPosition.setSize(nClosestPointsUsed);
        for(int n = 0; n < nClosestPointsUsed; n++)
        {
            nMinDist[n] = VGREAT;
            pointPosition[n] = vector(VGREAT, VGREAT, VGREAT);
        }

        for(int procI = 0; procI < Pstream::nProcs(); procI++)
        {
            for(int fromI = 0; fromI < fromPointsAll[procI].size(); fromI++)
            {
                vector fromPoint = fromPointsAll[procI][fromI];
                scalar dist = mag
                (
                    toPoints[toI]
                  - fromPoint
                );

                if(dist < nMinDist[nClosestPointsUsed-1])
                {
                    label position(nClosestPointsUsed-1);

                    // Find "position"
                    bool check(true);
                    if(position == 0)
                    {
                        check = false;
                    }
                    while(check)
                    {
                        if(dist > nMinDist[position-1])
                        {
                            check = false;
                        }
                        else
                        {
                            position--;
                        }
                        if(position == 0)
                        {
                            check = false;
                        }
                    }

                    // Downgrade indices if node is new 
                    // (processor boundary nodes are duplicated !)
                    // TODO: is there a better way to detect processor 
                    // boundary nodes ?
                    scalar TOL = readScalar
                    (
                        pointSelectionDict().lookup("matchTolerance")
                    );

                    bool newPoint(true);
                    forAll(pointPosition,pP)
                    {
                        if
                        (
                            mag(pointPosition[pP].x() - fromPoint.x()) < TOL
                         && mag(pointPosition[pP].y() - fromPoint.y()) < TOL
                         && mag(pointPosition[pP].z() - fromPoint.z()) < TOL
                        )
                        {
                            newPoint = false;
                        }
                    }

                    if(newPoint)
                    {
                        for(label i = nClosestPointsUsed-1; i > position; i--)
                        {
                            nMinDist[i] = nMinDist[i-1];
                            pointPosition[i] = pointPosition[i-1];
                            fromPointsIDs[toI][i] = fromPointsIDs[toI][i-1];
                            fromPointsProcIDs[toI][i] = fromPointsProcIDs[toI][i-1];
                        }

                        // Assign new index
                        nMinDist[position] = dist;
                        pointPosition[position] = fromPoint;
                        fromPointsIDs[toI][position] = fromI;
                        fromPointsProcIDs[toI][position] = procI;
                    }
                }
            }
        }
    }
    Info<< "done." <<endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
