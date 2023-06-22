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
    fsiInterfaceDisplacement

Group
    fsiPointBoundaryConditions

Description
    This boundary condition is used to move the fluid interface according 
    to the solid deformation that is computed in the FSI process.
    
    The displacement value computed in the solid region is used to specify 
    the mesh position for the fluid region.
    
    calculated as:

        \f[
            pointD_{f} = pointD_{s} 
        \f]

    where
    \vartable
        pointD_{f}  | point displacement for the fluid patch [length]
        pointD_{s}  | point displacement from the solid patch [length]
    \endtable

Usage
    \table
        Property            | Description                       | Req'd? | Default
        solidRegion         | Corresponding solid region name   | yes    |
        solidPatch          | Corresponding solid patch name    | yes    |
        radialExtrapolation | Selection of the method           | yes    |
        pointSelection      | Selection of the method           | yes    |
    \endtable

    Example of the boundary conditions specification for a fluid region:
    \verbatim
    obstacleFluid
    {
        type            fsiInterfaceDisplacement;
        solidRegion     solid; // Name of the region on the other side  
    solidPatch      obstacleSolid; // Name of the patch on the other side 
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

#ifndef fsiInterfaceDisplacementPointPatchVectorField_H
#define fsiInterfaceDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchFields.H"
#include "fvMesh.H"

// #define PATCHTOPATCHINTERPOLATION patchToPatchInterpolation
// #include "patchToPatchInterpolation.H"
// #define PATCHTOPATCHINTERPOLATION newPatchToPatchInterpolation
// #include "newPatchToPatchInterpolation.H"
// #include "setPatchToPatchInterpolation.H"
//#include "patchToPatchDataTransfer.H"
//#define PATCHTOPATCHINTERPOLATION patchToPatchDataTransfer
#include "radialExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class fsiInterfaceDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fsiInterfaceDisplacementPointPatchVectorField
:
    /*public FixedValuePointPatchField
    <
        pointPatchField,
        pointMesh,
        pointPatch,
        DummyMatrix,
        vector
    >*/
    public fixedValuePointPatchField<vector>
{
    // Private data

        dictionary dict_;
        word solidRegionName_;
        word solidPatchName_;
        autoPtr<radialExtrapolation> interpolatorSolidFluid_;
        vectorField fsiD_;
        vectorField fsiDPrevIter_;
        /*vectorField fsiDPrevIterPrevIter_;
        vectorField fsiDOldTime_;
        vectorField fsiDOldTimeOldTime_;
        vectorField fsiDResidual_;
        scalarField omega_;*/
        bool updateInterpolator_;
        //label timeIndex_;
        //label outerIteration_;
        scalar scaleResidual_;

    // Private member functions

        // Constant coefficient relaxation.
        //void relaxDisplacementConstCoeff();

        // Aitken relaxation. Note: the fluid field must 
        // be computed before the solid field.
        //void relaxDisplacementAitken();

public:

    //- Runtime type information
    TypeName("fsiInterfaceDisplacement");


    // Constructors

        //- Construct from patch and internal field
        fsiInterfaceDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fsiInterfaceDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        fsiInterfaceDisplacementPointPatchVectorField
        (
            const fsiInterfaceDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            //const PointPatchFieldMapper&
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new fsiInterfaceDisplacementPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fsiInterfaceDisplacementPointPatchVectorField
        (
            const fsiInterfaceDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new fsiInterfaceDisplacementPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );

        // Evaluation functions

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
