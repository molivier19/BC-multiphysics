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

#include "fstiTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "multiFieldIterator.H"

#include "fstiHeatFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fstiTemperatureFvPatchScalarField::setInterpolator()
{
    const fvMesh& meshSolid(db().time().lookupObject<fvMesh>(solidRegionName_));
    label solidPatchID(meshSolid.boundaryMesh().findPatchID(solidPatchName_));

    Info<< "Setting interpolator" << endl;

    interpolatorSolidFluid_ = radialExtrapolation::New
    (
        meshSolid.boundary()[solidPatchID].Cf(),
        patch().Cf(),
        dict_
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fstiTemperatureFvPatchScalarField::
fstiTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    temperatureCoupledBase
    (
        patch(),
        "undefined",
        "undefined",
        "undefined-K",
        "undefined-alpha"
    ),
    TnbrName_("undefined-Tnbr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0),

    dict_(nullptr),
    TSolid_(p.size(), 0.0),
    DkSolid_(p.size(), 0.0),
    solidRegionName_("0"),
    solidPatchName_("0"),
    interpolatorSolidFluid_(nullptr),
    updateInterpolator_(true)
{
    fvPatchScalarField::operator=(patchInternalField());
}


fstiTemperatureFvPatchScalarField::
fstiTemperatureFvPatchScalarField
(
    const fstiTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(),ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),

    dict_(ptf.dict_),
    TSolid_(ptf.TSolid_, mapper),
    DkSolid_(ptf.DkSolid_, mapper),
    solidRegionName_(ptf.solidRegionName_),
    solidPatchName_(ptf.solidPatchName_),
    interpolatorSolidFluid_(nullptr),
    updateInterpolator_(true)
{}


fstiTemperatureFvPatchScalarField::
fstiTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    temperatureCoupledBase(patch(),dict),
    TnbrName_(dict.get<word>("Tnbr")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0),

    dict_(dict),
    TSolid_(p.size(), 0.0),
    DkSolid_(p.size(), 0.0),
    solidRegionName_(dict.lookup("solidRegion")),
    solidPatchName_(dict.lookup("solidPatch")),
    interpolatorSolidFluid_(nullptr),
    updateInterpolator_(true)
{
    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }
    
    fvPatchScalarField::operator=(patchInternalField());
}


fstiTemperatureFvPatchScalarField::
fstiTemperatureFvPatchScalarField
(
    const fstiTemperatureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    temperatureCoupledBase(patch(),ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),

    dict_(ptf.dict_),
    TSolid_(ptf.TSolid_),
    DkSolid_(ptf.DkSolid_),
    solidRegionName_(ptf.solidRegionName_),
    solidPatchName_(ptf.solidPatchName_),
    interpolatorSolidFluid_(nullptr),
    updateInterpolator_(true)
{}


fstiTemperatureFvPatchScalarField::
fstiTemperatureFvPatchScalarField
(
    const fstiTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    temperatureCoupledBase(patch(),ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),

    dict_(ptf.dict_),
    TSolid_(ptf.TSolid_),
    DkSolid_(ptf.DkSolid_),
    solidRegionName_(ptf.solidRegionName_),
    solidPatchName_(ptf.solidPatchName_),
    interpolatorSolidFluid_(nullptr),
    updateInterpolator_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fstiTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    TSolid_.autoMap(m);
    DkSolid_.autoMap(m);
}


void fstiTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const fstiTemperatureFvPatchScalarField& tiptf =
        refCast<const fstiTemperatureFvPatchScalarField>(ptf);

    TSolid_.rmap(tiptf.TSolid_, addr);
    DkSolid_.rmap(tiptf.DkSolid_, addr);
}


void fstiTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& meshFluid(patch().boundaryMesh().mesh());
    const fvMesh& meshSolid(db().time().lookupObject<fvMesh>(solidRegionName_));
    const label& solidPatchID(meshSolid.boundaryMesh().findPatchID(solidPatchName_));
    const fvPatch& solidPatch(meshSolid.boundary()[solidPatchID]);

    const multiFieldIterator& mfi
    (
        db().time().lookupObject<multiFieldIterator>
        (
            "iterator"
        )
    );

    const fstiHeatFluxFvPatchScalarField& solidField =
                                refCast<const fstiHeatFluxFvPatchScalarField>
                                (
                                    solidPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
                                );
                                
    // Note: interpolating directly (T or GradT) causes a memory
    // allocation problem. We need an actual scalar Field, not a pointer.
    // Hence the field tmpTraction.

    vectorField n = solidPatch.nf();
    const fvPatchField<tensor>& gradD = solidPatch.lookupPatchField<volTensorField, tensor>("gradD");
    
    scalarField tmpTSolid(solidField.patchInternalField()); // T_Solid
    scalarField tmpDXSolid(solidPatch.deltaCoeffs());
    scalarField tmpDxSolid(tmpDXSolid*mag(n & inv(I+gradD.T())));
    scalarField tmpDkSolid(solidField.kappa(solidField)*tmpDxSolid); // Delta_KSolid
                     

    Info<< "Setting temperature on patch \""
        << patch().name() << "\" for outer-iteration "
        << mfi.index() << endl;

    if(updateInterpolator_)
    {
        // Move the mesh in the current configuration.
        Info<< endl << "Moving the mesh " << meshSolid.dbDir() << " in the current configuration" << endl;
        {
            const pointVectorField& pointD = meshSolid.lookupObject<pointVectorField>("pointD");
            vectorField newPoints = meshSolid.points() + pointD.prevIter().primitiveField();
            const_cast<fvMesh&>(meshSolid).movePoints(newPoints);
        }

        setInterpolator();

        TSolid_ = interpolatorSolidFluid_->pointInterpolate(tmpTSolid);
        DkSolid_ = interpolatorSolidFluid_->pointInterpolate(tmpDkSolid);

        // Move the mesh back in the initial configuration.
        Info<< "Moving the mesh " << meshSolid.dbDir() << " back in the initial configuration" << endl;
        {
            const pointVectorField& pointD = meshSolid.lookupObject<pointVectorField>("pointD");
            vectorField newPoints = meshSolid.points() - pointD.prevIter().primitiveField();
            const_cast<fvMesh&>(meshSolid).movePoints(newPoints);
        }

        updateInterpolator_ = false;
    }
    else
    {
        TSolid_ = interpolatorSolidFluid_->pointInterpolate(tmpTSolid);
        DkSolid_ = interpolatorSolidFluid_->pointInterpolate(tmpDkSolid);
    }

    scalarField TFluid = patchInternalField();
    scalarField DkFluid = kappa(*this)*patch().deltaCoeffs();

    scalarField::operator=(
                                (
                                    TSolid_ // T_Solid
                                    * DkSolid_ // Delta_KSolid 
                                    + TFluid // T_Fluid
                                    * DkFluid // Delta_KFluid
                                )/(
                                    DkSolid_ // Delta_KSolid
                                    + DkFluid // Delta_KFluid
                                )
                            );

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< meshFluid.name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << meshSolid.name() << ':'
            << solidPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << " Kappa Value "
            << " min:" << gMin(kappa(*this))
            << " max:" << gMax(kappa(*this))
            << " avg:" << gAverage(kappa(*this))
            << endl;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void fstiTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    thicknessLayers_.writeEntry("thicknessLayers", os);
    kappaLayers_.writeEntry("kappaLayers", os);

    os.writeKeyword("solidRegion")
        << solidRegionName_ << token::END_STATEMENT << nl;
    os.writeKeyword("solidPatch")
        << solidPatchName_ << token::END_STATEMENT << nl;
    
    os.writeKeyword(radialExtrapolation::typeName);
    dict_.subDict(radialExtrapolation::typeName).write(os);

    temperatureCoupledBase::write(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fstiTemperatureFvPatchScalarField
);

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
