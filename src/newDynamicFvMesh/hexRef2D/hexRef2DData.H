
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
    Foam::HexRef2DData

Description (by Luca)
    Various for reading/decomposing/reconstructing/distributing refinement
    data.

Original Author
    Luca Cornolti - written for foam-extend

Modified for OpenFOAM-v1806 by
    Ajit Kumar, Shiv Nadar University, India

Modified for OpenFOAM-v2006 by

SourceFiles
    hexRef2DData.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef2DData_H
#define hexRef2DData_H

#include "labelIOList.H"
#include "uniformDimensionedFields.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class mapDistributePolyMesh;
class refinementHistory;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class hexRef2DData Declaration
\*---------------------------------------------------------------------------*/

class hexRef2DData
{
    // Private Data

        autoPtr<labelIOList> cellLevelPtr_;

        autoPtr<labelIOList> pointLevelPtr_;

        autoPtr<uniformDimensionedScalarField> level0EdgePtr_;

        autoPtr<refinementHistory> refHistoryPtr_;


    // Private Member Functions

        //- No copy construct
        hexRef2DData(const hexRef2DData&) = delete;

        //- No copy assignment
        void operator=(const hexRef2DData&) = delete;


public:

    // Constructors

        //- Construct read. Has special provision for only some processors
        //- having the files so can be used in redistribution.
        explicit hexRef2DData(const IOobject& io);

        //- Construct as subset
        hexRef2DData
        (
            const IOobject& io,
            const hexRef2DData&,
            const labelList& cellMap,
            const labelList& pointMap
        );

        //- Construct from multiple hexRef2DData
        hexRef2DData
        (
            const IOobject& io,
            const UPtrList<const labelList>& cellMaps,
            const UPtrList<const labelList>& pointMaps,
            const UPtrList<const hexRef2DData>&
        );


    //- Destructor
    ~hexRef2DData();


    // Member Functions

        //- Parallel synchronise. This enforces valid objects on all processors
        //  (even if they don't have a mesh). Used by redistributePar.
        void sync(const IOobject& io);

        //- Update local numbering for changed mesh.
        void updateMesh(const mapPolyMesh&);

        //- In-place distribute
        void distribute(const mapDistributePolyMesh&);

        //- Write
        bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
