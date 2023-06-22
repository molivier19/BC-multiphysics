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
    pointSelectionSharingZ <-- pointSelection

Description
    This point-to-point selection method select the points that share the 
    same position in the z direction.

SourceFiles
    pointSelectionSharingZ.C

\*---------------------------------------------------------------------------*/

#ifndef pointSelectionSharingZ_H
#define pointSelectionSharingZ_H

#include "pointSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class pointSelectionSharingZ Declaration
\*---------------------------------------------------------------------------*/

class pointSelectionSharingZ
:
    public pointSelection
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointSelectionSharingZ(const pointSelectionSharingZ&);

        //- Disallow default bitwise assignment
        void operator=(const pointSelectionSharingZ&);

public:

    //- Runtime type information
    TypeName("pointSelectionSharingZ");

    // Constructors

        //- Construct from components
        pointSelectionSharingZ
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~pointSelectionSharingZ();


    // Member Functions

        //- Generate tables

        virtual void generatePointIDs
        (
            labelListList& fromPointsIDs,
            labelListList& fromPointsProcIDs,
            const pointField& fromPoints,
            const pointField& toPoints
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
