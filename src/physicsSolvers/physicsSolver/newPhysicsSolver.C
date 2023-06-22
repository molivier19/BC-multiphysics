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

\*---------------------------------------------------------------------------*/

#include "physicsSolver.H"
#include "polyMesh.H"
#include "Time.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::physicsSolver> Foam::physicsSolver::New
(
    const IOobject& io,
    const word& regionName
)
{
    word rn = regionName;
    if(regionName == polyMesh::defaultRegion)
    {
        rn = "";
    }
    
    IOdictionary physicsSolverDict
    (
        IOobject
        (
            "physicsSolverDict",
            io.time().constant(),
            rn,
            io.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    const word physicsSolverTypeName(physicsSolverDict.lookup("physicsSolver"));

    Info<< "Selecting physicsSolver " << physicsSolverTypeName << endl;

    const_cast<Time&>(io.time()).libs().open
    (
        physicsSolverDict,
        "physicsSolverLibs",
        IOobjectConstructorTablePtr_
    );

    if (!IOobjectConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "physicsSolver::New(const IOobject&)"
        )   << "physicsSolver table is empty"
            << exit(FatalError);
    }

#if OPENFOAM >= 2112
    auto cstrIter = IOobjectConstructorTablePtr_->find(physicsSolverTypeName);
#else
    IOobjectConstructorTable::iterator cstrIter =
        IOobjectConstructorTablePtr_->find(physicsSolverTypeName);
#endif


    if (cstrIter == IOobjectConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "physicsSolver::New(const IOobject&)"
        )   << "Unknown physicsSolver type " << physicsSolverTypeName
            << nl << nl
            << "Valid physicsSolver types are :" << endl
            << IOobjectConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<physicsSolver>(cstrIter()(io, regionName));
}


// ************************************************************************* //
