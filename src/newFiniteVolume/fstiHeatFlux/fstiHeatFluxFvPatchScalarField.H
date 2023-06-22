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
    Foam::compressible::fstiHeatFluxFvPatchScalarField

Description
    Boundary conditions for temperature. The temperature is computed with
    information from both sides of the interface. To be used for heat-transfer
    problems. This algorithm also takes into account the large deformation
    formulation and the St-Venant-Kirchhoff constitutive law. Was design
    specifically for FSI simulations with heat transfer. Optional thin thermal
    layer resistance can be specified through "thicknessLayers" and "kappaLayers"
    entries.

        \f[
            \frac{(T_{s}-T_{f}) \frac{\kappa_{f} \Delta X_{s}}{\Delta x_{f}}}
            {\frac{\kappa_{s}|n \cdot F^{-T}|}{\Delta X_{s}} + \frac{\kappa_{f}}{\Delta x_{f}}}
        \f]

    where
    \vartable

        _{s}     | indice for the solid region                 [-]
        _{f}     | indice for the fluid region                 [-]
        T        | Cell temperature                            [K]
        \kappa   | Thermal conductivity                        [Wm^{-1}K^{-1}]
        F        | Deformation Gradient                        [-]
        n        | Surface normal vector                       [-]
        \Delta X | Surface to centre distance ref. config.     [m]
        \Delta x | Surface to centre distance deformed config. [m]
    \endtable

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class "temperatureCoupledBase".

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivities per layer [W/m/K] | no |
        kappaMethod  | inherited from temperatureCoupledBase | inherited |
        kappa        | inherited from temperatureCoupledBase | inherited |
        solidRegion         | Name of the region the other side | yes    |
        solidPatch          | Name of the patch the other side  | yes    |
        radialExtrapolation | Selection of the method           | yes    |
        pointSelection      | Selection of the method           | yes    |
    \endtable

    Example of the boundary conditions specification:
    \verbatim
    obstacleSolid
    {
        type            compressible::fstiHeatFlux;
        Tnbr            T;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        kappaMethod     lookup;
        kappa           kappa;
        solidRegion     fluid;
        solidPatch      obstacleFluid;
        radialExtrapolation
        {
            method "radialExtrapolationNNI";
            pointSelectionDict
            {
                method pointSelectionNearestPoint;
            }
        }
        value           uniform 300;
    }
    \endverbatim

SourceFiles
    fstiHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fstiHeatFluxFvPatchScalarField_H
#define fstiHeatFluxFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "radialExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
     Class fstiHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fstiHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    public temperatureCoupledBase
{
    // Private data

        //- temperature specific data
        const word TnbrName_; //- Name of field on the neighbour region
        scalarList thicknessLayers_; //- Thickness of layers
        scalarList kappaLayers_; //- Conductivity of layers
        scalar contactRes_; //- Total contact resistance
        
        //- multiFieldSolver specific data
        dictionary dict_;
        scalarField TFluid_;
        scalarField DkFluid_;
        word fluidRegionName_;
        word fluidPatchName_;
        autoPtr<radialExtrapolation> interpolatorFluidSolid_;
        bool updateInterpolator_;
        label iterator_;
        label timeIndex_;
        
    // Private Member Functions

        void setInterpolator();

public:

    //- Runtime type information
    TypeName("compressible::fstiHeatFlux");

    // Constructors

        //- Construct from patch and internal field
        fstiHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fstiHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        fstiHeatFluxFvPatchScalarField
        (
            const fstiHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fstiHeatFluxFvPatchScalarField
        (
            const fstiHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fstiHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fstiHeatFluxFvPatchScalarField
        (
            const fstiHeatFluxFvPatchScalarField&,
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
                new fstiHeatFluxFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

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

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
