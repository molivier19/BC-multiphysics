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
    radialExtrapolationIDWXY <-- radialExtrapolation

Description
    This radial extrapolation class uses inverse-distance weighting
    interpolation in the xy plane.

SourceFiles
    radialExtrapolationIDWXY.C

\*---------------------------------------------------------------------------*/

#ifndef radialExtrapolationIDWXY_H
#define radialExtrapolationIDWXY_H

#include "radialExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class radialExtrapolationIDWXY Declaration
\*---------------------------------------------------------------------------*/

class radialExtrapolationIDWXY
:
    public radialExtrapolation
{
    // Private Member Functions

        dictionary IDWdict_;

        scalar exponent_;

        //- Disallow default bitwise copy construct
        radialExtrapolationIDWXY(const radialExtrapolationIDWXY&);

        //- Disallow default bitwise assignment
        void operator=(const radialExtrapolationIDWXY&);

public:

    //- Runtime type information
    TypeName("radialExtrapolationIDWXY");

    // Constructors

        //- Construct from components
        radialExtrapolationIDWXY
        (
            const pointField&fromPoints,
            const pointField& toPoints,
            const dictionary& dict
        );


    // Destructor

        virtual ~radialExtrapolationIDWXY();


    // Member Functions

        //- Generate tables

        virtual void generatePointWeights();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
