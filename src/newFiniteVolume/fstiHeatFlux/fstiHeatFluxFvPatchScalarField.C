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

#include "fstiHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "multiFieldIterator.H"

#include "fstiTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fstiHeatFluxFvPatchScalarField::setInterpolator()
{
    const fvMesh& meshFluid(db().time().lookupObject<fvMesh>(fluidRegionName_));
    label fluidPatchID(meshFluid.boundaryMesh().findPatchID(fluidPatchName_));

    Info<< "Setting interpolator" << endl;

    interpolatorFluidSolid_ = radialExtrapolation::New
    (
        meshFluid.boundary()[fluidPatchID].Cf(),
        patch().Cf(),
        dict_
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fstiHeatFluxFvPatchScalarField::
fstiHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
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
    TFluid_(p.size(), 0.0),
    DkFluid_(p.size(), 0.0),
    fluidRegionName_("0"),
    fluidPatchName_("0"),
    interpolatorFluidSolid_(nullptr),
    updateInterpolator_(true),
    iterator_(-1),
    timeIndex_(0)
{
    fvPatchScalarField::operator=(patchInternalField());
    gradient() = 0.0;
}


fstiHeatFluxFvPatchScalarField::
fstiHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    TnbrName_(dict.get<word>("Tnbr")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0),

    dict_(dict),
    TFluid_(p.size(), 0.0),
    DkFluid_(p.size(), 0.0),
    fluidRegionName_(dict.lookup("fluidRegion")),
    fluidPatchName_(dict.lookup("fluidPatch")),
    interpolatorFluidSolid_(nullptr),
    updateInterpolator_(true),
    iterator_(-1),
    timeIndex_(0)
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
    gradient() = 0.0;
}


fstiHeatFluxFvPatchScalarField::
fstiHeatFluxFvPatchScalarField
(
    const fstiHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),

    dict_(ptf.dict_),
    TFluid_(ptf.TFluid_, mapper),
    DkFluid_(ptf.DkFluid_, mapper),
    fluidRegionName_(ptf.fluidRegionName_),
    fluidPatchName_(ptf.fluidPatchName_),
    interpolatorFluidSolid_(nullptr),
    updateInterpolator_(true),
    iterator_(ptf.iterator_),
    timeIndex_(ptf.timeIndex_)
{}


fstiHeatFluxFvPatchScalarField::
fstiHeatFluxFvPatchScalarField
(
    const fstiHeatFluxFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    temperatureCoupledBase(patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),

    dict_(ptf.dict_),
    TFluid_(ptf.TFluid_),
    DkFluid_(ptf.DkFluid_),
    fluidRegionName_(ptf.fluidRegionName_),
    fluidPatchName_(ptf.fluidPatchName_),
    interpolatorFluidSolid_(nullptr),
    updateInterpolator_(true),
    iterator_(ptf.iterator_),
    timeIndex_(ptf.timeIndex_)
{}


fstiHeatFluxFvPatchScalarField::
fstiHeatFluxFvPatchScalarField
(
    const fstiHeatFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    temperatureCoupledBase(patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),

    dict_(ptf.dict_),
    TFluid_(ptf.TFluid_),
    DkFluid_(ptf.DkFluid_),
    fluidRegionName_(ptf.fluidRegionName_),
    fluidPatchName_(ptf.fluidPatchName_),
    interpolatorFluidSolid_(nullptr),
    updateInterpolator_(true),
    iterator_(ptf.iterator_),
    timeIndex_(ptf.timeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fstiHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    TFluid_.autoMap(m);
    DkFluid_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fstiHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const fstiHeatFluxFvPatchScalarField& dmptf =
        refCast<const fstiHeatFluxFvPatchScalarField>(ptf);

    TFluid_.rmap(dmptf.TFluid_, addr);
    DkFluid_.rmap(dmptf.DkFluid_, addr);
}


// Update the coefficients associated with the patch field

void fstiHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (patch().boundaryMesh().mesh().time().timeIndex() != timeIndex_)
    {
        timeIndex_ = db().time().timeIndex();
        iterator_ = -1;
    }

    const fvMesh& meshSolid(patch().boundaryMesh().mesh());
    const fvMesh& meshFluid(db().time().lookupObject<fvMesh>(fluidRegionName_));
    const label& fluidPatchID(meshFluid.boundaryMesh().findPatchID(fluidPatchName_));
    const fvPatch& fluidPatch(meshFluid.boundary()[fluidPatchID]);


    if (db().time().foundObject<multiFieldIterator>("iterator"))
    {
        const multiFieldIterator& mfi
        (
            db().time().lookupObject<multiFieldIterator>
            (
                "iterator"
            )
        );

        if(mfi.index() != iterator_)
        {
            iterator_ = mfi.index();

            const fstiTemperatureFvPatchScalarField& fluidField =
                                refCast<const fstiTemperatureFvPatchScalarField>
                                (
                                    fluidPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
                                );

            // Note: interpolating directly (T or GradT) causes a memory
            // allocation problem. We need an actual scalar Field, not a pointer.
            // Hence the field tmpTraction.

            scalarField tmpTFluid(fluidField.patchInternalField());
            scalarField tmpDKFluid(fluidField.kappa(fluidField)*fluidPatch.deltaCoeffs());

            // ===================== Detection of Load balancing on neighbour mesh ================================
            if (updateInterpolator_ == false && Pstream::nProcs() > 0)
            {
                scalarList procInterpolator(Pstream::nProcs(), 0.0);
                scalarList procTField(Pstream::nProcs(), 0.0);
                scalarList procKField(Pstream::nProcs(), 0.0);

                procInterpolator[Pstream::myProcNo()] = interpolatorFluidSolid_().size();
                procTField[Pstream::myProcNo()] = tmpTFluid.size();
                procKField[Pstream::myProcNo()] = tmpTFluid.size();

                reduce(procInterpolator, sumOp<List<scalar>>());
                reduce(procTField, sumOp<List<scalar>>());
                reduce(procKField, sumOp<List<scalar>>());

                for (int i = 0; i < Pstream::nProcs(); i++)
                {
                    if (procInterpolator[i] != procTField[i] || procInterpolator[i] != procKField[i])
                    {
                        updateInterpolator_ = true;
                        if (debug)
                        {
                            Info << "Update interpolator! " << endl;
                        }
                    }
                }
            }
            // =====================================================================================================
            
            Info << "Setting Heat Flux on patch \""
                 << patch().name() 
                 << "\" for outer-iteration "
                 << mfi.index() << endl;

            if(updateInterpolator_)
            {
                // Move the mesh in the current configuration.
                Info << "Moving the mesh " 
                     << meshSolid.dbDir() 
                     << " in the current configuration" << endl;
                {
                    const pointVectorField& pointD = meshSolid.lookupObject<pointVectorField>("pointD");
                    vectorField newPoints = meshSolid.points() + pointD.primitiveField();
                    const_cast<fvMesh&>(meshSolid).movePoints(newPoints);
                }

                setInterpolator();

                TFluid_ = interpolatorFluidSolid_->pointInterpolate(tmpTFluid);
                DkFluid_ = interpolatorFluidSolid_->pointInterpolate(tmpDKFluid);

                // Move the mesh back in the initial configuration.
                Info << "Moving the mesh " 
                     << meshSolid.dbDir() 
                     << " back in the initial configuration" << endl;
                {
                    const pointVectorField& pointD = meshSolid.lookupObject<pointVectorField>("pointD");
                    vectorField newPoints = meshSolid.points() - pointD.primitiveField();
                    const_cast<fvMesh&>(meshSolid).movePoints(newPoints);
                }

                updateInterpolator_ = false;
            }
            else
            {
                TFluid_ = interpolatorFluidSolid_->pointInterpolate(tmpTFluid);
                DkFluid_ = interpolatorFluidSolid_->pointInterpolate(tmpDKFluid);
            }
        }
    }

    vectorField n = patch().nf();
    const fvPatchField<tensor>& gradD = patch().lookupPatchField<volTensorField, tensor>("gradD");

    scalarField TSolid = patchInternalField();
    scalarField DXSolid = patch().deltaCoeffs();
    scalarField DxSolid = DXSolid*mag(n & inv(I+gradD.T()));
    scalarField DkSolid = kappa(*this)*DxSolid;
    scalarField aASolid = mag(n & inv(I+gradD.T()))*det(I+gradD.T());

    gradient() = (
                    TFluid_ // T_fluid
                    - TSolid // T_Solid
                )*
                (
                    DkFluid_ // Delta_Kfluid 
                    * DXSolid // Delta_X
                )/
                (
                    DkSolid // Delta_KSolid
                    + DkFluid_ // Delta_Kfluid
                );

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*aASolid*snGrad()*mag(n & inv(I+gradD.T())));

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << meshFluid.name() << ':'
            << fluidPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature ::"
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << " Kappa Value ::"
            << " min:" << gMin(kappa(*this))
            << " max:" << gMax(kappa(*this))
            << " avg:" << gAverage(kappa(*this))
            << endl;
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void fstiHeatFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    thicknessLayers_.writeEntry("thicknessLayers", os);
    kappaLayers_.writeEntry("kappaLayers", os);

    os.writeKeyword("fluidRegion")
        << fluidRegionName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluidPatch")
        << fluidPatchName_ << token::END_STATEMENT << nl;

    os.writeKeyword(radialExtrapolation::typeName);
    dict_.subDict(radialExtrapolation::typeName).write(os);

    temperatureCoupledBase::write(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //


    makePatchTypeField
    (
        fvPatchScalarField,
        fstiHeatFluxFvPatchScalarField
    );
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
