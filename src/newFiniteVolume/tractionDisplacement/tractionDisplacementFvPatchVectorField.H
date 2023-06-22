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
    tractionDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation. The equation was adapted for the large 
    deformation formulation and the St.Venant-Kirchhoff constitutive law.

    calculated as:

        \f[
           \nabla D = \frac{\left( T + P \frac{n \cdot F^{-1}}{|(n \cdot F^{-T})(J n \cdot F^{-T})|}\right)\cdot F^{-T}+2\mu\lambda*\nabla D-(n \cdot \sigma)}{2\mu\lambda};
        \f]

    where
    \vartable
        T       | Other traction forces     [forces]
        P       | Other pressure forces     [forces]
        n       | Surface normal vectors    [-]
        F       | Deformation Tensor        [-]
        J       | Det of Deformation Tensor [-]
        \lambda | first lamé coefficient    [Pa]
        \mu     | second lamé coefficient   [Pa]
        D       | Deformation vector        [length]
        \sigma  | Cauchy tensor             [Pa]
    \endtablei

Usage
    \table
        Property            | Description                       | Req'd? | Default
    \endtable

    Example of the boundary conditions specification:
    \verbatim
    obstacleSolid
    {
        type            tractionDisplacement;
        traction        uniform ( 0 0 0 );
        pressure        uniform 0;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    tractionDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tractionDisplacementFvPatchVectorField_H
#define tractionDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class tractionDisplacementFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        vectorField traction_;
        scalarField pressure_;


public:

    //- Runtime type information
    TypeName("tractionDisplacement");


    // Constructors

        //- Construct from patch and internal field
        tractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  tractionDisplacementFvPatchVectorField onto a new patch
        tractionDisplacementFvPatchVectorField
        (
            const tractionDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionDisplacementFvPatchVectorField
        (
            const tractionDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionDisplacementFvPatchVectorField
        (
            const tractionDisplacementFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual  scalarField& pressure()
            {
                return pressure_;
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
                const fvPatchVectorField&,
                const labelList&
            );


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
