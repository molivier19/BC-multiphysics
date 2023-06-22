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
    physicsSolverTemplate <-- physicsSolver

Description
    Class solver template for new physical solver.

SourceFiles
    physicsSolverTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef physicsSolverTemplate_H
#define physicsSolverTemplate_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "physicsSolver.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
//class someClass;

/*---------------------------------------------------------------------------*\
                           Class physicsSolverTemplate Declaration
\*---------------------------------------------------------------------------*/

class physicsSolverTemplate
:
    public physicsSolver
{

protected:

    // Protected data
        // ================================================
        // All physical parameters required for the physical 
        // solver should be implemented here.

        // Mesh
        autoPtr<dynamicFvMesh> meshPtr_;
        dynamicFvMesh& mesh_;

        // Dictionary
        IOdictionary properties_;

        // Field
        autoPtr<volVectorField> Uptr_;
        autoPtr<volScalarField> rhoptr_;

        autoPtr<volVectorField> divU_;
        autoPtr<volScalarField> divRho_;
        
        // ================================================

    // Protected Member Functions
        // ================================================
        // Optional custom functions

        // Functions
        void computeExplicitFluxes();
        bool validateIODictionary();
        
        // ================================================

        //- Disallow default bitwise copy construct
        physicsSolverTemplate(const physicsSolverTemplate&);

        //- Disallow default bitwise assignment
        void operator=(const physicsSolverTemplate&);


public:

    //- Runtime type information
        TypeName("physicsSolverTemplate");


    // Constructors

        //- Construct from components
        physicsSolverTemplate
        (
            const IOobject& io, 
            const word& regionName
        );


    // Destructor

        virtual ~physicsSolverTemplate();


    // Member Functions

        virtual void createFields(); // Required

        virtual bool isConverged(); // Required

        virtual void solveFields(); // Required

        virtual void startTimeStep(); // Optional

        virtual void endTimeStep(); // Optional

        virtual bool computeDeltaT(); // Optional

        virtual scalar returnDeltaT(); // Optional
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
