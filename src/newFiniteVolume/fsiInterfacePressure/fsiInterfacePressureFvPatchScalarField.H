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
    fsiInterfacePressure

Group
    fsiBoundaryConditions

Description
    This boundary condition is used to stabilize the FSI simulation 
    with an artificial compressibility method. It introduces terms in
    the matrix which correspond to artificial compressibility
    located in the cells next to the FSI interface.

Usage
    \table
        Property            | Description                       | Req'd? | Default
        DeltaP              | Test-load pressure                | yes    |
    \endtable

    Example of the boundary conditions specification for a fluid region:
    \verbatim
    obstacleFluid
    {
        type            fsiInterfacePressure;
    DeltaP            1e-6;     
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    fsiInterfacePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fsiInterfacePressureFvPatchScalarField_H
#define fsiInterfacePressureFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class fsiInterfacePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fsiInterfacePressureFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
private:

    // Private data

        label iterator_;

        scalarField psi_;

        scalarField pPrevIter_;

        scalarField phipPrevIter_;

        bool updatePsi_;

        Switch noStabOnFirstIteration_;

        label testLoadIndex_;

        scalar relaxFactor_;

        scalar deltaP_;


public:

    //- Runtime type information
    TypeName("fsiInterfacePressure");


    // Constructors

        //- Construct from patch and internal field
        fsiInterfacePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fsiInterfacePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fsiInterfacePressureFvPatchScalarField
        //  onto a new patch
        fsiInterfacePressureFvPatchScalarField
        (
            const fsiInterfacePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fsiInterfacePressureFvPatchScalarField
        (
            const fsiInterfacePressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fsiInterfacePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fsiInterfacePressureFvPatchScalarField
        (
            const fsiInterfacePressureFvPatchScalarField&,
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
                new fsiInterfacePressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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
            virtual void updateCoeffs();

            //- Evaluate the patch field
            // virtual void evaluate
            // (
            //     const Pstream::commsTypes commsType =
            //         Pstream::commsTypes::blocking
            // );

            ////- Return the matrix diagonal coefficients corresponding to the
            ////  evaluation of the value of this patchField with given weights
            //virtual tmp<scalarField> valueInternalCoeffs
            //(
            //    const tmp<scalarField>&
            //) const;

            ////- Return the matrix source coefficients corresponding to the
            ////  evaluation of the value of this patchField with given weights
            //virtual tmp<scalarField> valueBoundaryCoeffs
            //(
            //    const tmp<scalarField>&
            //) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientBoundaryCoeffs() const;


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
