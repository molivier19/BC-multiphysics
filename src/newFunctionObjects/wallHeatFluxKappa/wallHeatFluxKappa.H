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
    Foam::functionObjects::wallHeatFluxKappa

Group
    grpForcesFunctionObjects

Description
    Computes the wall-heat flux at selected wall patches.

    \table
      Operand        | Type              | Location
      input          | -                 | -
      output file    | dat               <!--
                 --> | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<field\>
      output field   | volScalarField (only boundaryField) <!--
                 --> | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wallHeatFluxKappa1
    {
        // Mandatory entries (unmodifiable)
        type        wallHeatFluxKappa;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        patches     (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");
        qr          qr;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: wallHeatFluxKappa               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patches    | Names of operand patches   | wordList | no | all wall patches
      qr         | Name of radiative heat flux field       | word | no | qr
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func wallHeatFluxKappa
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::wallHeatFluxKappa

SourceFiles
    wallHeatFluxKappa.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallHeatFluxKappa_H
#define functionObjects_wallHeatFluxKappa_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallHeatFluxKappa Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFluxKappa
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Optional list of wall patches to process
        labelHashSet patchSet_;

        //- Name of radiative heat flux name
        word qrName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the heat-flux
        void calcHeatFlux
        (
            const volScalarField& kappa,
            const volScalarField& T,
            volScalarField& wallHeatFlux
        );


public:

    //- Runtime type information
    TypeName("wallHeatFluxKappa");


    // Constructors

        //- Construct from Time and dictionary
        wallHeatFluxKappa
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        wallHeatFluxKappa(const wallHeatFluxKappa&) = delete;

        //- No copy assignment
        void operator=(const wallHeatFluxKappa&) = delete;


    //- Destructor
    virtual ~wallHeatFluxKappa() = default;


    // Member Functions

        //- Read the wallHeatFluxKappa data
        virtual bool read(const dictionary& dict);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
