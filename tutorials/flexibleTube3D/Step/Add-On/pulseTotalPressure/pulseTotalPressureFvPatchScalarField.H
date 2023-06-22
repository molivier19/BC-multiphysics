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
    pulseTotalPressureFvPatchScalarField

Description
    Inlet pressure boundary condition to impose a total pressure pulsation.

Usage
    \table
        Property            | Description                       | Req'd? | Default
        U                   | Name of the velocity field        |        |
        phi                 | Name of the phi field             |        |
        rho                 | Name of the rho field             |        |
        psi                 | Name of the psi field             |        |
        gamma               | Gamma value                       |        |
        time                | Pulsation time                    |        |
        pPulse              | Peak pressure                     |        |
        p0                  | Reference pressure                |        |
    \endtable

    Example of the boundary conditions specification for a fluid region:
    \verbatim
    obstacleSolid
    {
        type            pulseTotalPressure;
        U               U;
        phi             phi;
        rho             none;
        psi             none;
        gamma           1;
        time            0.003;
        pPulse          uniform 1.3332;
        p0              uniform 0;
        value           uniform 1.3332;
    }
    \endverbatim

SourceFiles
    pulseTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pulseTotalPressureFvPatchScalarField_H
#define pulseTotalPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pulseTotalPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class pulseTotalPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;

        scalar time_;
        scalarField pPulse_;

        //- Total pressure
        scalarField p0_;


public:

    //- Runtime type information
    TypeName("pulseTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        pulseTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pulseTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pulseTotalPressureFvPatchScalarField
        //  onto a new patch
        pulseTotalPressureFvPatchScalarField
        (
            const pulseTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pulseTotalPressureFvPatchScalarField
        (
            const pulseTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pulseTotalPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pulseTotalPressureFvPatchScalarField
        (
            const pulseTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new pulseTotalPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the heat capacity ratio
            scalar gamma() const
            {
                return gamma_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& gamma()
            {
                return gamma_;
            }

            //- Return the total pressure
            const scalarField& p0() const
            {
                return p0_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& p0()
            {
                return p0_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  using the given patch velocity field
            virtual void updateCoeffs(const vectorField& Up);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
