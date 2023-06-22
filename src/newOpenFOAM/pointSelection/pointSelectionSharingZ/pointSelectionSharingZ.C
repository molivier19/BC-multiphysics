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

#include "pointSelectionSharingZ.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointSelectionSharingZ, 0);
addToRunTimeSelectionTable
(
    pointSelection,
    pointSelectionSharingZ,
    patches
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pointSelectionSharingZ::pointSelectionSharingZ
(
    const dictionary& dict
)
:
    pointSelection(dict)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

pointSelectionSharingZ::~pointSelectionSharingZ(){}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointSelectionSharingZ::generatePointIDs
(
    labelListList& fromPointsIDs,
    labelListList& fromPointsProcIDs,
    const pointField& fromPoints,
    const pointField& toPoints
)
{
//    scalar nClosestPointsUsed =
//        readLabel(pointSelectionDict().lookup("nClosestPointsUsed"));

    Info<< "Finding " << " all points..." <<endl;


    fromPointsIDs.clear();
    fromPointsProcIDs.clear();

    fromPointsIDs.setSize(toPoints.size());
    fromPointsProcIDs.setSize(toPoints.size());

/*    for(int n = 0; n < toPoints.size(); n++)
    {
        fromPointsIDs[n].clear();
        fromPointsProcIDs[n].clear();

        fromPointsIDs[n].setSize(nClosestPointsUsed);
        fromPointsProcIDs[n].setSize(nClosestPointsUsed);

//        for(label nn = 0; nn < nClosestPointsUsed; nn++)
//        {
//            fromPointsIDs[n][nn] = 0;
//            fromPointsProcIDs[n][nn] = 0;
//        }
    }
*/


    List< pointField > fromPointsAll(Pstream::nProcs());
    fromPointsAll[Pstream::myProcNo()] = fromPoints;
    Pstream::gatherList(fromPointsAll);
    Pstream::scatterList(fromPointsAll);

    // Processor boundary nodes are duplicated...
    List< List<bool> > isExcluded(fromPointsAll.size());
    forAll(isExcluded, procI)
    {
        isExcluded[procI].clear();
        isExcluded[procI].setSize(fromPointsAll[procI].size());
        forAll(isExcluded[procI], pI)
        {
            isExcluded[procI][pI] = false;
        }
    }

    scalar TOL = readScalar(pointSelectionDict().lookup("matchTolerance"));
    scalar zValue = readScalar(pointSelectionDict().lookup("zValue"));

    // Find and exclude out-of-plane and duplicated nodes
    label nFromPoints(0);
    forAll(fromPointsAll, procI)
    {
        forAll(fromPointsAll[procI], pI)
        {

            if
            (
                mag
                (
                    fromPointsAll[procI][pI].z()
                  - zValue
                ) > TOL
            )
            {
                isExcluded[procI][pI] = true;
            }
            else
            {
                // Compare with nodes in all previously swept 
                // proc (up to procI-1)
                for(label procJ = 0; procJ < procI; procJ++)
                {
                    forAll(fromPointsAll[procJ], pJ)
                    {
                        if
                        (
                            mag
                            (
                                fromPointsAll[procJ][pJ]
                              - fromPointsAll[procI][pI]
                            ) < TOL
                        )
                        {
                            isExcluded[procI][pI] = true;
                        }
                    }
                }
                for(label pJ = 0; pJ < pI; pJ++)
                {
                    if
                    (
                        mag
                        (
                            fromPointsAll[procI][pJ]
                          - fromPointsAll[procI][pI]
                        ) < TOL
                    )
                    {
                        isExcluded[procI][pI] = true;
                    }
                }
            }

            if(!isExcluded[procI][pI])
            {
                nFromPoints++;
            }
        }
    }

    Info<< "Detected " << nFromPoints 
        << " different points in the \"from\" cloud." << endl;

    for(int n = 0; n < toPoints.size(); n++)
    {
        fromPointsIDs[n].clear();
        fromPointsProcIDs[n].clear();

        fromPointsIDs[n].setSize(nFromPoints);
        fromPointsProcIDs[n].setSize(nFromPoints);

//        for(label nn = 0; nn < nFromPoints; nn++)
//        {
//            fromPointsIDs[n][nn] = 0;
//            fromPointsProcIDs[n][nn] = 0;
//        }
    }

    for(int toI = 0; toI < toPoints.size(); toI++)
    {
        label fromI(0);
        forAll(fromPointsAll, procI)
        {
            forAll(fromPointsAll[procI], pI)
            {
                if(!isExcluded[procI][pI])
                {
                    fromPointsIDs[toI][fromI] = pI;
                    fromPointsProcIDs[toI][fromI] = procI;
                    fromI++;
                }
            }
        }
        
    }
    Info<< "done." <<endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
