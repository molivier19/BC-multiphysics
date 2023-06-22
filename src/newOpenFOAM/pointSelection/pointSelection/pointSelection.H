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
    pointSelection

Description
    Base class for point set selection classes. The point sets are used for
    interpolation or extrapolation purposes (e.g. FSI interfaces, mesh motion).

SourceFiles
    pointSelection.C
    newPointSelection.C

\*---------------------------------------------------------------------------*/

#ifndef pointSelection_H
#define pointSelection_H

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
                Class pointSelection Declaration
\*---------------------------------------------------------------------------*/

class pointSelection
:
    public dictionary
{
protected:

        dictionary pointSelectionDict_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointSelection(const pointSelection&);

        //- Disallow default bitwise assignment
        void operator=(const pointSelection&);


public:

    //- Runtime type information
    TypeName("pointSelection");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointSelection,
            patches,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Selectors

        //- Return a pointer to the selected point selection method
        static autoPtr<pointSelection> New
        (
            const dictionary& dict
        );

    // Constructors

        //- Construct from components
        pointSelection
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~pointSelection();


    // Member Functions

        //- Generate tables
        virtual void generatePointIDs
        (
            labelListList& fromPointsIDs,
            labelListList& fromPointsProcIDs,
            const pointField& fromPoints,
            const pointField& toPoints
        ) = 0;

        // Write

        void write(Ostream& os, bool subDict=true) const
        {
            os.writeKeyword("pointSelection");
            dictionary::write(os, subDict);
        }

        virtual bool fromPointsSubsetIsUnique()
        {
            return false;
        }

        const dictionary& pointSelectionDict()
        {
            return pointSelectionDict_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
