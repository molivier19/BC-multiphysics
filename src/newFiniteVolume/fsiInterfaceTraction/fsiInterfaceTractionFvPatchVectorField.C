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

#include "fsiInterfaceTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "multiFieldIterator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fsiInterfaceTractionFvPatchVectorField::setInterpolator()
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

fsiInterfaceTractionFvPatchVectorField::
fsiInterfaceTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    dict_(nullptr),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    fsiTraction_(p.size(), vector::zero),
    fluidRegionName_("0"),
    fluidPatchName_("0"),
    interpolatorFluidSolid_(nullptr),
    timeIndex_(0),
    iterator_(-1),
    updateInterpolator_(true)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


fsiInterfaceTractionFvPatchVectorField::
fsiInterfaceTractionFvPatchVectorField
(
    const fsiInterfaceTractionFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    dict_(tdpvf.dict_),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper),
    fsiTraction_(tdpvf.fsiTraction_, mapper),
    fluidRegionName_(tdpvf.fluidRegionName_),
    fluidPatchName_(tdpvf.fluidPatchName_),
    interpolatorFluidSolid_(nullptr),
    timeIndex_(tdpvf.timeIndex_),
    iterator_(tdpvf.iterator_),
    updateInterpolator_(true)
{}


fsiInterfaceTractionFvPatchVectorField::
fsiInterfaceTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    dict_(dict),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size()),
    fsiTraction_(p.size(), vector::zero),
    fluidRegionName_(dict.lookup("fluidRegion")),
    fluidPatchName_(dict.lookup("fluidPatch")),
    interpolatorFluidSolid_(nullptr),
    timeIndex_(0),
    iterator_(-1),
    updateInterpolator_(true)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


fsiInterfaceTractionFvPatchVectorField::
fsiInterfaceTractionFvPatchVectorField
(
    const fsiInterfaceTractionFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    dict_(tdpvf.dict_),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    fsiTraction_(tdpvf.fsiTraction_),
    fluidRegionName_(tdpvf.fluidRegionName_),
    fluidPatchName_(tdpvf.fluidPatchName_),
    interpolatorFluidSolid_(nullptr),
    timeIndex_(tdpvf.timeIndex_),
    iterator_(tdpvf.iterator_),
    updateInterpolator_(true)
{}


fsiInterfaceTractionFvPatchVectorField::
fsiInterfaceTractionFvPatchVectorField
(
    const fsiInterfaceTractionFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    dict_(tdpvf.dict_),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    fsiTraction_(tdpvf.fsiTraction_),
    fluidRegionName_(tdpvf.fluidRegionName_),
    fluidPatchName_(tdpvf.fluidPatchName_),
    interpolatorFluidSolid_(nullptr),
    timeIndex_(tdpvf.timeIndex_),
    iterator_(tdpvf.iterator_),
    updateInterpolator_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fsiInterfaceTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fsiInterfaceTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const fsiInterfaceTractionFvPatchVectorField& dmptf =
        refCast<const fsiInterfaceTractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field

void fsiInterfaceTractionFvPatchVectorField::updateCoeffs()
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

            const fvMesh& meshFluid
            (
                db().time().lookupObject<fvMesh>(fluidRegionName_)
            );

            label fluidPatchID
            (
                meshFluid.boundaryMesh().findPatchID(fluidPatchName_)
            );

            const volSymmTensorField& sigmaFluid = 
                meshFluid.lookupObject<volSymmTensorField>("sigmaFluid");

            // Note: interpolating directly (sigma & n) causes a memory
            // allocation problem. We need an actual vectorField, not a pointer.
            // Hence the field tmpTraction.
            vectorField tmpTraction
            (
               -sigmaFluid.boundaryField()[fluidPatchID]
              & meshFluid.boundary()[fluidPatchID].nf()
            );

            // ===================== Detection of Load balancing on neighbour mesh ================================
            if (updateInterpolator_ == false && Pstream::nProcs() > 0)
            {
                scalarList procInterpolator(Pstream::nProcs(), 0.0);
                scalarList procTField(Pstream::nProcs(), 0.0);

                procInterpolator[Pstream::myProcNo()] = interpolatorFluidSolid_().size();
                procTField[Pstream::myProcNo()] = tmpTraction.size();

                reduce(procInterpolator, sumOp<List<scalar>>());
                reduce(procTField, sumOp<List<scalar>>());

                for (int i = 0; i < Pstream::nProcs(); i++)
                {
                    if (procInterpolator[i] != procTField[i])
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
            

            Info << "Setting traction on patch \""
                 << patch().name() 
                 << "\" for outer-iteration "
                 << mfi.index() << endl;

            if(updateInterpolator_)
            {
                // Move the mesh in the current configuration.
                Info << "Moving the mesh " 
                     << patch().boundaryMesh().mesh().dbDir() 
                     << " in the current configuration" << endl;
                {
                    const pointVectorField& pointD = 
                        patch().boundaryMesh().mesh().lookupObject<pointVectorField>("pointD");

                    vectorField newPoints = 
                        patch().boundaryMesh().mesh().points() 
                      + pointD.primitiveField();

                    const_cast<fvMesh&>(patch().boundaryMesh().mesh()).movePoints
                    (
                        newPoints
                    );
                }

                setInterpolator();

                fsiTraction_ = interpolatorFluidSolid_->pointInterpolate
                (
                    tmpTraction
                );

                // Move the mesh back in the initial configuration.
                Info << "Moving the mesh " 
                     << patch().boundaryMesh().mesh().dbDir() 
                     << " back in the initial configuration" << endl;
                {
                    const pointVectorField& pointD = 
                        patch().boundaryMesh().mesh().lookupObject<pointVectorField>("pointD");

                    vectorField newPoints = 
                        patch().boundaryMesh().mesh().points() 
                      - pointD.primitiveField();

                    const_cast<fvMesh&>(patch().boundaryMesh().mesh()).movePoints
                    (
                        newPoints
                    );
                }

                updateInterpolator_ = false;
            }
            else
            {
                fsiTraction_ = interpolatorFluidSolid_->pointInterpolate
                (
                    tmpTraction
                );
                
            }
        }
    }

    const dictionary& mechanicalProperties =
        db().lookupObject<IOdictionary>("mechanicalProperties");

    dimensionedScalar rho
    (
        "rho",
        dimDensity,
        mechanicalProperties
    );

    dimensionedScalar E
    (
        "E",
        dimPressure,
        mechanicalProperties
    );

    dimensionedScalar nu
    (
        "nu",
        dimless,
        mechanicalProperties
    );

    dimensionedScalar mu = E/(2.0*(1.0 + nu));
    dimensionedScalar lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
    dimensionedScalar threeK = E/(1.0 - 2.0*nu);

    Switch planeStress(mechanicalProperties.lookup("planeStress"));

    if (planeStress)
    {
        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
        threeK = E/(1.0 - nu);
    }
    scalar twoMuLambda = (2*mu + lambda).value();

    vectorField n = patch().nf();

    const fvPatchField<symmTensor>& sigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");
        
    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>("gradD");

    //scalarField SoS0 = mag(n & invF)*det(F);
    scalarField SoS0 = mag(n & inv(I+gradD.T()))*det(I+gradD.T());

    gradient() =
    (
        (((fsiTraction_ + traction_ + pressure_*((n & inv(I+gradD.T()))
      / mag(n & inv(I+gradD.T()))))*SoS0) & (inv(I+gradD)))
      + twoMuLambda*fvPatchField<vector>::snGrad() - (n & sigma) 
    )/twoMuLambda;

    if (debug)
    {
        Info << "Total traction force:" << gSum(fsiTraction_*patch().magSf()) << endl;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


// Write
void fsiInterfaceTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    os.writeKeyword("fluidRegion")
        << fluidRegionName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluidPatch")
        << fluidPatchName_ << token::END_STATEMENT << nl;

    /*******************************************************/
    // TODO
    // The following line does not work...
    // autoPtr not allocated ?!
    //interpolatorFluidSolid_->write(os);

    // See file fsiInterfaceDisplacementPointPatchVectorField.C...

    // These 2 lines are a not very elegant solution to the allocation problem.
    // Requires a dictionary (dict_) to be a member of this class.
    os.writeKeyword(radialExtrapolation::typeName);
    dict_.subDict(radialExtrapolation::typeName).write(os);
    /*******************************************************/

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, fsiInterfaceTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
