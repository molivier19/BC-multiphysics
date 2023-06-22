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
    radialExtrapolation

Description
    Base class for radial extrapolation (i.e. distance-based extrapolation).
    Can be used, for example, for FSI interface interpolation or mesh motion.

SourceFiles
    radialExtrapolation.C
    newRadialExtrapolation.C

\*---------------------------------------------------------------------------*/

#ifndef radialExtrapolation_H
#define radialExtrapolation_H

#include "runTimeSelectionTables.H"
#include "PrimitivePatch.H"
#include "face.H"
#include "SubList.H"
#include "pointField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class radialExtrapolation Declaration
\*---------------------------------------------------------------------------*/

class radialExtrapolation
:
    public dictionary
{
protected:
    // Protected data

        // Point tables
        List< List<label> > fromPointsIDs_;
        List< List<label> > fromPointsProcIDs_;
        List< List<scalar> > fromPointsWeights_;

        bool IDsNeedUpdate_;
        bool weightsNeedUpdate_;
        bool fromPointsSubsetIsUnique_;

        //- Reference to the source cloud
        pointField fromPoints_;

        //- Reference to the target cloud
        pointField toPoints_;

        label size_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        radialExtrapolation(const radialExtrapolation&);

        //- Disallow default bitwise assignment
        void operator=(const radialExtrapolation&);


public:

    //- Runtime type information
    TypeName("radialExtrapolation");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            radialExtrapolation,
            patches,
            (
                const pointField& fromPoints,
                const pointField& toPoints,
                const dictionary& dict
            ),
            (fromPoints, toPoints, dict)
        );

    // Selectors

        //- Return a pointer to the selected extrapolation method
        static autoPtr<radialExtrapolation> New
        (
            const pointField& fromPoints,
            const pointField& toPoints,
            const dictionary& dict
        );

    // Constructors

        //- Construct from components
        radialExtrapolation
        (
            const pointField& fromPoints,
            const pointField& toPoints,
            const dictionary& dict
        );


    // Destructor

        virtual ~radialExtrapolation();


    // Member Functions

        //- Interpolate point field
        template<class Type>
        Field<Type> pointInterpolate(const Field<Type>& pf);

        template<class Type>
        Field<Type> pointInterpolate(const tmp<Field<Type> >& tpf);

        //- Set table sizes
        //virtual void initTableSizes(label nPoints);

        //- Generate tables
        virtual void generatePointIDs();

        // Calculate point weights
        virtual void generatePointWeights() = 0;

        //- Access
        const List< List<scalar> >& fromPointsWeights() const
        {
            return fromPointsWeights_;
        }
        List< List<scalar> >& fromPointsWeights()
        {
            return fromPointsWeights_;
        }

        const pointField& fromPoints() const
        {
            return fromPoints_;
        }

        const pointField& toPoints() const
        {
            return toPoints_;
        }

        /*pointField& fromPoints()
        {
            return fromPoints_;
        }

        pointField& toPoints()
        {
            return toPoints_;
        }*/

        const List< List<label> > fromPointsIDs() const
        {
            return fromPointsIDs_;
        }
        List< List<label> > fromPointsIDs()
        {
            return fromPointsIDs_;
        }

        const List< List<label> > fromPointsProcIDs() const
        {
            return fromPointsProcIDs_;
        }
        List< List<label> > fromPointsProcIDs()
        {
            return fromPointsProcIDs_;
        }


        const label& size() const
        {
            return size_;
        }


        // Write

        void write(Ostream& os, bool subDict=true) const
        {
            os.writeKeyword("radialExtrapolation");
            dictionary::write(os, subDict);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radialExtrapolationInterpolate.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
