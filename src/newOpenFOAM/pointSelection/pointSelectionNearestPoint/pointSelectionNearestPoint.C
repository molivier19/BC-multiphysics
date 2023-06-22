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

#include "pointSelectionNearestPoint.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointSelectionNearestPoint, 0);
addToRunTimeSelectionTable
(
    pointSelection,
    pointSelectionNearestPoint,
    patches
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pointSelectionNearestPoint::pointSelectionNearestPoint
(
    const dictionary& dict
)
:
    pointSelection(dict)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

pointSelectionNearestPoint::~pointSelectionNearestPoint(){}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointSelectionNearestPoint::generatePointIDs
(
    labelListList& fromPointsIDs,
    labelListList& fromPointsProcIDs,
    const pointField& fromPoints,
    const pointField& toPoints
)
{
    Info<< "Finding " << " nearest points..." <<endl;

    fromPointsIDs.clear();
    fromPointsProcIDs.clear();

    fromPointsIDs.setSize(toPoints.size());
    fromPointsProcIDs.setSize(toPoints.size());

    for(int n = 0; n < toPoints.size(); n++)
    {
        fromPointsIDs[n].clear();
        fromPointsProcIDs[n].clear();

        fromPointsIDs[n].setSize(1);
        fromPointsProcIDs[n].setSize(1);

        fromPointsIDs[n][0] = 0;
        fromPointsProcIDs[n][0] = 0;
    }



    List< pointField > fromPointsAll(Pstream::nProcs());
    fromPointsAll[Pstream::myProcNo()] = fromPoints;
    Pstream::gatherList(fromPointsAll);
    // Pstream::scatterList(fromPointsAll, Pstream::blocking);
    Pstream::scatterList(fromPointsAll);

    for(int toI = 0; toI < toPoints.size(); toI++)
//    forAllReverse(toPoints, toI)
    {
        scalar minDist(VGREAT);


        for(int procI = 0; procI < Pstream::nProcs(); procI++)
//        forAllReverse(fromPointsAll, procI)
        {
            for(int fromI = 0; fromI < fromPointsAll[procI].size(); fromI++)
//            forAllReverse(fromPointsAll[procI], fromI)
            {
                vector fromPoint = fromPointsAll[procI][fromI];
                scalar dist = mag
                (
                    toPoints[toI]
                  - fromPoint
                );

                if(dist < minDist)
                {
                    minDist = dist;
                    fromPointsIDs[toI][0] = fromI;
                    fromPointsProcIDs[toI][0] = procI;
                }
            }

        }

        /* label procI = fromPointsProcIDs[toI][0]; */
        /* label fromI = fromPointsIDs[toI][0]; */
        /* if */
        /* ( */
        /*     mag */
        /*     ( */
        /*         toPoints[toI] */ 
        /*       - fromPointsAll[procI][fromI] */
        /*     ) > 1e-10 */
        /* ) */
        /* { */
        /*     WarningInFunction */
        /*         << "Nearest neighbor found is beyond tolerance (1e-10)" */
        /*         << endl */
        /*         << "Matching point distance: " << mag */
        /*         ( */
        /*             toPoints[toI] */ 
        /*           - fromPointsAll[procI][fromI] */
        /*         ) */
        /*         << endl */
        /*         << "toPoint: " << toPoints[toI] << endl */
        /*         << "fromPoint: " << fromPointsAll[procI][fromI]<< endl */
        /*         << endl; */ 
        /* } */
    }
 
    Info<< "done." <<endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
