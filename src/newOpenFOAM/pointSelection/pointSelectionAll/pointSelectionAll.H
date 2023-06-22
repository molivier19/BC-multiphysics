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
    pointSelectionAll <-- pointSelection 

Description
    This point-to-point selection method selects a point set made of all
    available points.

SourceFiles
    pointSelectionAll.C

\*---------------------------------------------------------------------------*/

#ifndef pointSelectionAll_H
#define pointSelectionAll_H

#include "pointSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class pointSelectionAll Declaration
\*---------------------------------------------------------------------------*/

class pointSelectionAll
:
    public pointSelection
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointSelectionAll(const pointSelectionAll&);

        //- Disallow default bitwise assignment
        void operator=(const pointSelectionAll&);

public:

    //- Runtime type information
    TypeName("pointSelectionAll");

    // Constructors

        //- Construct from components
        pointSelectionAll
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~pointSelectionAll();


    // Member Functions

        //- Generate tables

        virtual void generatePointIDs
        (
            labelListList& fromPointsIDs,
            labelListList& fromPointsProcIDs,
            const pointField& fromPoints,
            const pointField& toPoints
        );

        virtual bool fromPointsSubsetIsUnique()
        {
            return true;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
