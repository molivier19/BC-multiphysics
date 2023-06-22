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

#include "pulseTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pulseTotalPressureFvPatchScalarField::pulseTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("undefined"),
    phiName_("undefined"),
    rhoName_("undefined"),
    psiName_("undefined"),
    gamma_(0.0),
    time_(0.0),
    pPulse_(p.size(), 0.0),
    p0_(p.size(), 0.0)
{}


pulseTotalPressureFvPatchScalarField::pulseTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookup("U")),
    phiName_(dict.lookup("phi")),
    rhoName_(dict.lookup("rho")),
    psiName_(dict.lookup("psi")),
    gamma_(readScalar(dict.lookup("gamma"))),
    time_(readScalar(dict.lookup("time"))),
    pPulse_("pPulse", dict, p.size()),
    p0_("p0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        if(db().time().value() > time_)
        {
            fvPatchField<scalar>::operator=(p0_);
        }
        else
        {
            fvPatchField<scalar>::operator=(pPulse_);
        }
    }
}


pulseTotalPressureFvPatchScalarField::pulseTotalPressureFvPatchScalarField
(
    const pulseTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    time_(ptf.time_),
    pPulse_(ptf.pPulse_, mapper),
    p0_(ptf.p0_, mapper)
{}


pulseTotalPressureFvPatchScalarField::pulseTotalPressureFvPatchScalarField
(
    const pulseTotalPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    time_(tppsf.time_),
    pPulse_(tppsf.pPulse_),
    p0_(tppsf.p0_)
{}


pulseTotalPressureFvPatchScalarField::pulseTotalPressureFvPatchScalarField
(
    const pulseTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    time_(tppsf.time_),
    pPulse_(tppsf.pPulse_),
    p0_(tppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pulseTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    pPulse_.autoMap(m);
    p0_.autoMap(m);
}


void pulseTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const pulseTotalPressureFvPatchScalarField& tiptf =
        refCast<const pulseTotalPressureFvPatchScalarField>(ptf);

    pPulse_.rmap(tiptf.pPulse_, addr);
    p0_.rmap(tiptf.p0_, addr);
}


void pulseTotalPressureFvPatchScalarField::updateCoeffs(const vectorField& Up)
{
    if (updated())
    {
        return;
    }

    scalarField pTot(size());
    if(db().time().value() > time_)
    {
        pTot = p0_;
    }
    else
    {
        pTot = pPulse_;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    if (psiName_ == "none" && rhoName_ == "none")
    {
        operator==(pTot - 0.5*(1.0 - pos(phip))*magSqr(Up));
    }
    else if (rhoName_ == "none")
    {
        const fvPatchField<scalar>& psip =
            patch().lookupPatchField<volScalarField, scalar>(psiName_);

        if (gamma_ > 1.0)
        {
            scalar gM1ByG = (gamma_ - 1.0)/gamma_;

            operator==
            (
                pTot
               /pow
                (
                    (1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up)),
                    1.0/gM1ByG
                )
            );
        }
        else
        {
            operator==(pTot/(1.0 + 0.5*psip*(1.0 - pos(phip))*magSqr(Up)));
        }
    }
    else if (psiName_ == "none")
    {
        const fvPatchField<scalar>& rho =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(pTot - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));
    }
    else
    {
        FatalErrorIn
        (
            "pulseTotalPressureFvPatchScalarField::updateCoeffs()"
        )   << " rho or psi set inconsitently, rho = " << rhoName_
            << ", psi = " << psiName_ << '.' << nl
            << "    Set either rho or psi or neither depending on the "
               "definition of total pressure." << nl
            << "    Set the unused variables to 'none'."
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void pulseTotalPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs(patch().lookupPatchField<volVectorField, vector>(UName_));
}


void pulseTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    os.writeKeyword("time") << time_ << token::END_STATEMENT << endl;
    pPulse_.writeEntry("pPulse", os);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, pulseTotalPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
