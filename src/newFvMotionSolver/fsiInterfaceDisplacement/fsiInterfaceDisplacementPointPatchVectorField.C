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

#include "fsiInterfaceDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "pointFields.H"
#include "multiFieldIterator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fsiInterfaceDisplacementPointPatchVectorField::
fsiInterfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    dict_(),
    solidRegionName_("0"),
    solidPatchName_("0"),
    interpolatorSolidFluid_(nullptr),
    fsiD_(p.size(), vector::zero),
    fsiDPrevIter_(p.size(), vector::zero),
    updateInterpolator_(true),
    scaleResidual_(1.0)
{}


fsiInterfaceDisplacementPointPatchVectorField::
fsiInterfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    dict_(dict),
    solidRegionName_(dict.lookup("solidRegion")),
    solidPatchName_(dict.lookup("solidPatch")),
    interpolatorSolidFluid_(nullptr),
    fsiD_(p.size(), vector::zero),
    fsiDPrevIter_(p.size(), vector::zero),
    updateInterpolator_(true),
    scaleResidual_(1.0)
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


fsiInterfaceDisplacementPointPatchVectorField::
fsiInterfaceDisplacementPointPatchVectorField
(
    const fsiInterfaceDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    solidRegionName_(ptf.solidRegionName_),
    solidPatchName_(ptf.solidPatchName_),
    interpolatorSolidFluid_(nullptr),
    fsiD_(ptf.fsiD_),
    fsiDPrevIter_(p.size(), vector::zero),
    updateInterpolator_(true),
    scaleResidual_(1.0)
{}


fsiInterfaceDisplacementPointPatchVectorField::
fsiInterfaceDisplacementPointPatchVectorField
(
    const fsiInterfaceDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    dict_(ptf.dict_),
    solidRegionName_(ptf.solidRegionName_),
    solidPatchName_(ptf.solidPatchName_),
    interpolatorSolidFluid_(nullptr),
    fsiD_(ptf.fsiD_),
    fsiDPrevIter_(ptf.size(), vector::zero),
    updateInterpolator_(true),
    scaleResidual_(1.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fsiInterfaceDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    fsiD_.autoMap(m);
    fsiDPrevIter_.autoMap(m);
}


void fsiInterfaceDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const fsiInterfaceDisplacementPointPatchVectorField& aODptf =
        refCast<const fsiInterfaceDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    fsiD_.rmap(aODptf.fsiD_, addr);
    fsiDPrevIter_.rmap(aODptf.fsiDPrevIter_, addr);
}


void fsiInterfaceDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const multiFieldIterator& mfi
    (
        db().time().lookupObject<multiFieldIterator>
        (
            "iterator"
        )
    );

    const polyMesh& meshSolid
    (
        db().time().lookupObject<polyMesh>(solidRegionName_)
    );
    label solidPatchID(meshSolid.boundaryMesh().findPatchID(solidPatchName_));

    const pointVectorField& pointD = 
        meshSolid.lookupObject<pointVectorField>("pointD");

    Info<< "Setting displacement on patch \""
        << patch().name() << "\" for outer-iteration "
        << mfi.index()
        << endl;

    if(updateInterpolator_)
    {
        // Move the mesh in the current configuration.
        {
            Info<< endl << "Moving the mesh " 
                << meshSolid.dbDir()
                << " in the current configuration" << endl;

            // Using pointD.prevIter() for cases where the solid 
            // field is solved first. This allows to have the exact
            // same displacement as the fluid interface.
            vectorField newPoints = meshSolid.points()
                + pointD.prevIter().primitiveField();
            const_cast<polyMesh&>(meshSolid).movePoints(newPoints);
        }
        Info<< "Setting interpolator" << endl;
        interpolatorSolidFluid_.clear();

        interpolatorSolidFluid_ = radialExtrapolation::New
        (
            meshSolid.boundaryMesh()[solidPatchID].localPoints(),
            patch().boundaryMesh().mesh()()
                .boundaryMesh()[patch().index()].localPoints(),
            dict_
        );

        updateInterpolator_ = false;

        fsiD_ = interpolatorSolidFluid_->pointInterpolate
        (
            pointD.boundaryField()[solidPatchID].patchInternalField()()
        );

        // Move the mesh back in the initial configuration.
        {
            Info<< "Moving the mesh " 
                << meshSolid.dbDir()
                << " back in the initial configuration" << endl;

            // Using pointD.prevIter() for cases where the solid 
            // field is solved first. This allows to have the exact
            // same displacement as the fluid interface.
            vectorField newPoints = meshSolid.points()
                - pointD.prevIter().primitiveField();
            const_cast<polyMesh&>(meshSolid).movePoints(newPoints);
        }
    }
    else
    {
        fsiD_ = interpolatorSolidFluid_->pointInterpolate
        (
            pointD.boundaryField()[solidPatchID].patchInternalField()()
        );
    }

    Info<< "Motion magnitude: mean = " 
        << gAverage(mag(fsiD_))
        << " max = " 
        << gMax(mag(fsiD_)) << endl;

    if(mfi.index() == 0)
    {
        scaleResidual_ = gMax(mag(fsiD_ - fsiDPrevIter_));
    }
    else
    {
        scaleResidual_ = max(gMax(mag(fsiD_ - fsiDPrevIter_)),scaleResidual_);
    }

    Info<< "Interface displacement residual: " 
        << gMax(mag(fsiD_ - fsiDPrevIter_))/max(scaleResidual_, VSMALL)
        << endl;

    fsiDPrevIter_ = fsiD_;

    vectorField::operator=(fsiD_);

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void fsiInterfaceDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("solidRegion")
        << solidRegionName_ << token::END_STATEMENT << nl;
    os.writeKeyword("solidPatch")
        << solidPatchName_ << token::END_STATEMENT << nl;
    interpolatorSolidFluid_->write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    fsiInterfaceDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
