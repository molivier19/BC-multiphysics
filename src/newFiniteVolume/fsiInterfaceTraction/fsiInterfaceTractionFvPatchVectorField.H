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
    fsiInterfaceTraction

Group
    fsiPointBoundaryConditions

Description
    FSI interface traction boundary condition for a solid solver used within 
    a FSI problem. The formulation uses the St.Venant-Kirchhoff constitutive 
    law. The boundary condition computes the traction from the interacting
    fluid flow on the solid side of the FSI interface. It is also possible to
    add an additional pressure or traction on the interface.

    calculated as:

        \f[
           \nabla D = \frac{\left(T_{fsi} + T + P \frac{n \cdot F^{-1}}{|(n \cdot F^{-T})(J n \cdot F^{-T})|}\right)\cdot F^{-T}+2\mu\lambda*\nabla D-(n \cdot \sigma)}{2\mu\lambda};
        \f]

    where
    \vartable
        T_{fsi} | Aerodynamic forces        [forces]
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
        fluidRegion         | Corresponding fluid region name   | yes    |
        fluidPatch          | Corresponding fluid patch name    | yes    |
        radialExtrapolation | Selection of the method           | yes    |
        pointSelection      | Selection of the method           | yes    |
    \endtable

    Example of the boundary conditions specification for a fluid region:
    \verbatim
    obstacleSolid
    {
        type            fsiInterfaceTraction;
        fluidRegion     fluid; // Name of the region on the other side
        fluidPatch      obstacleFluid; // Name of the patch on the other side
    traction        uniform ( 0 0 0 );
        pressure        uniform 0;
        value           uniform (0 0 0);
        radialExtrapolation // See the radialExtrapolation folder for more info
        {
            method  "radialExtrapolationNNI";
            pointSelectionDict // See the pointSelection folder for more info
            {
                method pointSelectionNearestPoint;
            }
        }
    }
    \endverbatim

SourceFiles
    fsiInterfaceDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fsiInterfaceTractionFvPatchVectorField_H
#define fsiInterfaceTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "radialExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fsiInterfaceTractionFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fsiInterfaceTractionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        dictionary dict_;
        vectorField traction_;
        scalarField pressure_;
        vectorField fsiTraction_;
        word fluidRegionName_;
        word fluidPatchName_;
        autoPtr<radialExtrapolation> interpolatorFluidSolid_;
        label timeIndex_;
        label iterator_;
        bool updateInterpolator_;

    // Private member function

        void setInterpolator();

public:

    //- Runtime type information
    TypeName("fsiInterfaceTraction");


    // Constructors

        //- Construct from patch and internal field
        fsiInterfaceTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fsiInterfaceTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fsiInterfaceTractionFvPatchVectorField onto a new patch
        fsiInterfaceTractionFvPatchVectorField
        (
            const fsiInterfaceTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fsiInterfaceTractionFvPatchVectorField
        (
            const fsiInterfaceTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fsiInterfaceTractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fsiInterfaceTractionFvPatchVectorField
        (
            const fsiInterfaceTractionFvPatchVectorField&,
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
                new fsiInterfaceTractionFvPatchVectorField(*this, iF)
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
