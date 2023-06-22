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
    radialExtrapolationRBF <-- radialExtrapolation

Description
    This radial extrapolation class uses radial basis function (RBF)
    interpolation.

SourceFiles
    radialExtrapolationRBF.C

\*---------------------------------------------------------------------------*/

#ifndef radialExtrapolationRBF_H
#define radialExtrapolationRBF_H

#include "radialExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class radialExtrapolationRBF Declaration
\*---------------------------------------------------------------------------*/

class radialExtrapolationRBF
:
    public radialExtrapolation
{
    // Private data

        dictionary RBFdict_;

        scalar radius_;

        bool xTerm_;
        bool yTerm_;
        bool zTerm_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        radialExtrapolationRBF(const radialExtrapolationRBF&);

        //- Disallow default bitwise assignment
        void operator=(const radialExtrapolationRBF&);

public:

    //- Runtime type information
    TypeName("radialExtrapolationRBF");

    // Constructors

        //- Construct from components
        radialExtrapolationRBF
        (
            const pointField& fromPoints,
            const pointField& toPoints,
            const dictionary& dict
        );


    // Destructor

        virtual ~radialExtrapolationRBF();


    // Member Functions

        // TODO: selection mechanism to choose the RBF. 
        // TODO: Using RBF from OF-1.6-ext ?
        scalar RBF(scalar);

        const bool& xTerm() const
        {
            return xTerm_;
        }
        const bool& yTerm() const
        {
            return yTerm_;
        }
        const bool& zTerm() const
        {
            return zTerm_;
        }

        virtual void generatePointWeights();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
