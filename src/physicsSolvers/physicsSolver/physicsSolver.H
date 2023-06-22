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
    physicsSolver

Description
    This class is a base class to create normalized physical solver classestea.
    
SourceFiles
    physicsSolver.C
    newPhysicsSolver.C
    
\*---------------------------------------------------------------------------*/

#ifndef physicsSolver_H
#define physicsSolver_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class physicsSolver Declaration
\*---------------------------------------------------------------------------*/

class physicsSolver
:
    public IOdictionary
{

protected:
    
        // Region
        word region_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        physicsSolver(const physicsSolver&);

        //- Disallow default bitwise assignment
        void operator=(const physicsSolver&);


public:

    //- Runtime type information
    TypeName("physicsSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            physicsSolver,
            IOobject,
            (
                const IOobject& io,
                const word& regionName
            ),
            (io, regionName)
        );


    // Constructors

        //- Construct from objectRegistry
        physicsSolver
        (
            const IOobject& io,
            const word& regionName
        );


    // Selectors

        //- Select null constructed
        static autoPtr<physicsSolver> New
        (
            const IOobject& io, const word& regionName
        );


    // Destructor

        virtual ~physicsSolver();


    // Member Functions

        virtual void createFields() = 0;
        
        virtual bool isConverged() = 0;

        virtual void solveFields() = 0;

        virtual void startTimeStep(){};

        virtual void endTimeStep(){};

        virtual bool computeDeltaT()
        {
            return false;
        }

        virtual scalar returnDeltaT()
        {
            return 0.0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
