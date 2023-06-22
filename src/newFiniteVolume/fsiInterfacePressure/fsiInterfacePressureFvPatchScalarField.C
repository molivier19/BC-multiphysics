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

#include "fsiInterfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multiFieldIterator.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fsiInterfacePressureFvPatchScalarField::
fsiInterfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    iterator_(-1),
    psi_(p.size(), 1.0),
    pPrevIter_(p.size(), 0.0),
    phipPrevIter_(p.size(), 0.0),
    updatePsi_(true),
    noStabOnFirstIteration_(false),
    testLoadIndex_(1),
    relaxFactor_(1.0),
    deltaP_(1.0)
{}


Foam::fsiInterfacePressureFvPatchScalarField::
fsiInterfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    iterator_(-1),
    psi_(p.size(), 0.0),
    pPrevIter_(p.size(), 0.0),
    phipPrevIter_(p.size(), 0.0),
    updatePsi_(true),
    noStabOnFirstIteration_(false),
    testLoadIndex_(1),
    relaxFactor_(1.0),
    deltaP_(readScalar(dict.lookup("deltaP")))
{
    if (dict.found("psi"))
    {
        psi_ = scalarField(dict.lookup("psi"));
        updatePsi_ = false;
        //
        //TODO: keep?
        Info<< "Min(psi) / Max(psi): " << gMin(psi_)<<" / "<< gMax(psi_) << endl;
        Info<< "gAverage(psi): " << gAverage(psi_) << endl;
    }
}


Foam::fsiInterfacePressureFvPatchScalarField::
fsiInterfacePressureFvPatchScalarField
(
    const fsiInterfacePressureFvPatchScalarField& fiptpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(fiptpsf, p, iF, mapper),
    iterator_(fiptpsf.iterator_),
    psi_(fiptpsf.psi_, mapper),
    pPrevIter_(fiptpsf.pPrevIter_, mapper),
    phipPrevIter_(fiptpsf.phipPrevIter_, mapper),
    updatePsi_(fiptpsf.updatePsi_),
    noStabOnFirstIteration_(fiptpsf.noStabOnFirstIteration_),
    testLoadIndex_(fiptpsf.testLoadIndex_),
    relaxFactor_(fiptpsf.relaxFactor_),
    deltaP_(fiptpsf.deltaP_)
{}


Foam::fsiInterfacePressureFvPatchScalarField::
fsiInterfacePressureFvPatchScalarField
(
    const fsiInterfacePressureFvPatchScalarField& fiptpsf
)
:
    zeroGradientFvPatchScalarField(fiptpsf),
    iterator_(fiptpsf.iterator_),
    psi_(fiptpsf.psi_),
    pPrevIter_(fiptpsf.pPrevIter_),
    phipPrevIter_(fiptpsf.phipPrevIter_),
    updatePsi_(fiptpsf.updatePsi_),
    noStabOnFirstIteration_(fiptpsf.noStabOnFirstIteration_),
    testLoadIndex_(fiptpsf.testLoadIndex_),
    relaxFactor_(fiptpsf.relaxFactor_),
    deltaP_(fiptpsf.deltaP_)
{}


Foam::fsiInterfacePressureFvPatchScalarField::
fsiInterfacePressureFvPatchScalarField
(
    const fsiInterfacePressureFvPatchScalarField& fiptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(fiptpsf, iF),
    iterator_(fiptpsf.iterator_),
    psi_(fiptpsf.psi_),
    pPrevIter_(fiptpsf.pPrevIter_),
    phipPrevIter_(fiptpsf.phipPrevIter_),
    updatePsi_(fiptpsf.updatePsi_),
    noStabOnFirstIteration_(fiptpsf.noStabOnFirstIteration_),
    testLoadIndex_(fiptpsf.testLoadIndex_),
    relaxFactor_(fiptpsf.relaxFactor_),
    deltaP_(fiptpsf.deltaP_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fsiInterfacePressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    psi_.autoMap(m);
    pPrevIter_.autoMap(m);
    phipPrevIter_.autoMap(m);

    zeroGradientFvPatchScalarField::autoMap(m);
}


void Foam::fsiInterfacePressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    zeroGradientFvPatchScalarField::rmap(ptf, addr);

    const fsiInterfacePressureFvPatchScalarField& fiptpsf =
        refCast<const fsiInterfacePressureFvPatchScalarField>(ptf);

    psi_.rmap(fiptpsf.psi_, addr);
    pPrevIter_.rmap(fiptpsf.pPrevIter_, addr);
    phipPrevIter_.rmap(fiptpsf.phipPrevIter_, addr);
}


void Foam::fsiInterfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
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

            const dictionary& fsiSolutionDict =
                patch().boundaryMesh().mesh().solutionDict().subDict("fsi");

            label stabUpdateInterval = fsiSolutionDict.lookupOrDefault<label>
            (
                "stabUpdateInterval",
                INT_MAX
            );

            if
            (
                (db().time().timeIndex() - 1) % stabUpdateInterval == 0
             || db().time().timeIndex() == 1
            )
            {
                updatePsi_ = true;
            }
            else
            {
                updatePsi_ = false;
            }

            const volVectorField& U
            (
                patch().boundaryMesh().mesh().lookupObject<volVectorField>("U")
            );

            scalarField phip
            (
                patch().patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
            );

            noStabOnFirstIteration_ =
                Switch(fsiSolutionDict.lookup("noStabOnFirstIteration"));

            relaxFactor_ = readScalar
            (
                fsiSolutionDict.subDict("fsiBoundaries").
                    subDict(patch().name()).lookup("relaxFactor")
            );

            // NOTE: This assumes the solid is solved first.
            testLoadIndex_ = 0;
            label updatePsiIndex = 1;

            if (updatePsi_)
            {
                if (mfi.index() == updatePsiIndex)
                {
                    Info<<"Updating psi on patch " << patch().name() <<endl;

                    // TODO: psi could be a scalar!!!
                    psi_ = gAverage
                    (
                        mag
                        (
                            (phip - phipPrevIter_)/
                            // (deltaP_*patch().magSf()*db().time().deltaT().value())
                            ((pPrevIter_ - patchInternalField())
                                *patch().magSf()*db().time().deltaT().value())
                        )
                    );

                    //TODO: keep?
                    Info<< "Min(psi) / Max(psi): " << gMin(psi_)<<" / "<< gMax(psi_) << endl;
                    Info<< "gAverage(psi): " << gAverage(psi_) << endl;
                }
            }

            phipPrevIter_ = phip;
            pPrevIter_ = patchInternalField();
        }
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


// void Foam::fsiInterfacePressureFvPatchScalarField::evaluate
// (
//     const Pstream::commsTypes commsType
// )
// {
//     if (db().time().foundObject<multiFieldIterator>("iterator"))
//     {
//         const multiFieldIterator& mfi
//         (
//             db().time().lookupObject<multiFieldIterator>
//             (
//                 "iterator"
//             )
//         );

//         if (mfi.index() == testLoadIndex_)
//         {
//             if (!this->updated())
//             {
//                 this->updateCoeffs();
//             }
// Info<<"EVAL"<<endl;
//             scalarField::operator=
//             (
//                 pPrevPrevIter_ + deltaP_
//             );

//             fvPatchScalarField::evaluate();
//         }
//         else
//         {
//             Foam::zeroGradientFvPatchScalarField::evaluate(commsType);
//         }
//     }
//     else
//     {
//         Foam::zeroGradientFvPatchScalarField::evaluate(commsType);
//     }
// }




Foam::tmp<Foam::scalarField>
Foam::fsiInterfacePressureFvPatchScalarField::gradientInternalCoeffs() const
{
    if (db().time().foundObject<multiFieldIterator>("iterator"))
    {
        const multiFieldIterator& mfi
        (
            db().time().lookupObject<multiFieldIterator>
            (
                "iterator"
            )
        );

        if (mfi.index() == 0 && noStabOnFirstIteration_)
        {
            return tmp<scalarField>
            (
                new scalarField(this->size(), 0)
            );
        }
        else if (mfi.index() == testLoadIndex_ && updatePsi_)
        {
            return -(patch().deltaCoeffs());
        }
        else
        {
            return -psi_/relaxFactor_;
        }
    }
    else
    {
        return zeroGradientFvPatchScalarField::gradientInternalCoeffs();
    }
}


Foam::tmp<Foam::scalarField>
Foam::fsiInterfacePressureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    if (db().time().foundObject<multiFieldIterator>("iterator"))
    {
        const multiFieldIterator& mfi
        (
            db().time().lookupObject<multiFieldIterator>
            (
                "iterator"
            )
        );

        if (mfi.index() == 0 && noStabOnFirstIteration_)
        {
            return tmp<scalarField>
            (
                new scalarField(this->size(), 0)
            );
        }
        else if (mfi.index() == testLoadIndex_ && updatePsi_)
        {
            // TODO: test both deltaP and -deltaP for nonlinear effects
            return this->patch().deltaCoeffs()*(pPrevIter_ + deltaP_);
        }
        else
        {
            return psi_*pPrevIter_/relaxFactor_;
        }
    }
    else
    {
        return zeroGradientFvPatchScalarField::gradientBoundaryCoeffs();
    }
}


void Foam::fsiInterfacePressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    os.writeKeyword("deltaP") << deltaP_ << token::END_STATEMENT << nl;

    os.writeEntry("psi", psi_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fsiInterfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
