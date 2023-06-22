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
    radialExtrapolationNNI <-- radialExtrapolation

Description
    This radial extrapolation class uses nearest neighbour interpolation. 
    The value on a "toPoints" point/face corresponds to the closest point/face 
    value on the "fromPoints".

SourceFiles
    radialExtrapolationNNI.C

\*---------------------------------------------------------------------------*/

#ifndef radialExtrapolationNNI_H
#define radialExtrapolationNNI_H

#include "radialExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class radialExtrapolationNNI Declaration
\*---------------------------------------------------------------------------*/

class radialExtrapolationNNI
:
    public radialExtrapolation
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        radialExtrapolationNNI(const radialExtrapolationNNI&);

        //- Disallow default bitwise assignment
        void operator=(const radialExtrapolationNNI&);

public:

    //- Runtime type information
    TypeName("radialExtrapolationNNI");

    // Constructors

        //- Construct from components
        radialExtrapolationNNI
        (
            const pointField&fromPoints,
            const pointField& toPoints,
            const dictionary& dict
        );


    // Destructor

        virtual ~radialExtrapolationNNI();


    // Member Functions

        //- Generate tables

        virtual void generatePointWeights();

        virtual void generatePointIDs();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
