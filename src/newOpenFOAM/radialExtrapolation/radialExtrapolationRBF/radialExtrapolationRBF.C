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

#include "radialExtrapolationRBF.H"
#include "simpleMatrix.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(radialExtrapolationRBF, 0);
addToRunTimeSelectionTable
(
    radialExtrapolation,
    radialExtrapolationRBF,
    patches
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
radialExtrapolationRBF::radialExtrapolationRBF
(
    const pointField& fromPoints,
    const pointField& toPoints,
    const dictionary& dict
)
:
    radialExtrapolation(fromPoints, toPoints, dict),
    RBFdict_(subDict("RBFparameters")),
    radius_(readScalar(RBFdict_.lookup("radius"))),
    xTerm_(RBFdict_.lookupOrDefault("xTerm", true)),
    yTerm_(RBFdict_.lookupOrDefault("yTerm", true)),
    zTerm_(RBFdict_.lookupOrDefault("zTerm", true))
{
    //TODO: nClosestPointsUsed schould not be defined here...
//    nClosestPointsUsed_ = readLabel(RBFdict_.lookup("nClosestPointsUsed"));
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

radialExtrapolationRBF::~radialExtrapolationRBF(){}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// TODO: selection mechanism to choose the RBF. 
// TODO: Using RBF from OF-1.6-ext ?
scalar radialExtrapolationRBF::RBF(scalar r)
{
    //return exp(-r/radius_);
    //return exp(-r*r/(radius_*radius_));
    return r*r/(radius_*radius_)*log10(r/radius_ + VSMALL);
    //return sqrt(r*r + 0.001);
    //return pow(pos(1 - r/radius_)*(1 - r/radius_), 4.0)*(4.0*r/radius_ + 1.0);
}


void radialExtrapolationRBF::generatePointWeights()
{
    Info<< "Generating weights using RBF..." << endl;


    List< pointField > fromPointsAll(Pstream::nProcs());
    fromPointsAll[Pstream::myProcNo()] = fromPoints_;
    Pstream::gatherList(fromPointsAll);
    Pstream::scatterList(fromPointsAll);

    label nCoeffs(1);
    if(xTerm_){nCoeffs++;}
    if(yTerm_){nCoeffs++;}
    if(zTerm_){nCoeffs++;}


    if(fromPointsSubsetIsUnique_)
    {
        // les tables sont identiques pour tous les noeuds internes si
        // "pointSelectionAll" (ou toute autre méthode qui fournit un ensemble unique).
        // L'indice 0 est utilisé arbitrairement...
        // TODO: optimiser la mémoire de "pointSelection" dans de tels cas:
        // seulement besoin de stocker l'indice 0.
        labelList fromPointsProcIDsPI = fromPointsProcIDs_[0];
        labelList fromPointsIDsPI = fromPointsIDs_[0];

        simpleMatrix<scalar> A(fromPointsIDsPI.size() + nCoeffs, 0.0, 0.0);
        A.source()[0] = 1.0;

        for(label i = 0; i < fromPointsIDsPI.size(); i++)
        {
            for(label j = 0; j < fromPointsIDsPI.size(); j++)
            {
                A[i+nCoeffs][j+nCoeffs] = RBF
                (
                    mag
                    (
                        fromPointsAll
                            [fromPointsProcIDsPI[i]]
                            [fromPointsIDsPI[i]]
                      - fromPointsAll
                            [fromPointsProcIDsPI[j]]
                            [fromPointsIDsPI[j]]
                    )
                );
            }

            A[0][i+nCoeffs] = 1.0;
            A[i+nCoeffs][0] = 1.0;

            label xyzIndex(1);

            if(xTerm_)
            {
                A[xyzIndex][i+nCoeffs] = fromPointsAll
                    [fromPointsProcIDsPI[i]]
                    [fromPointsIDsPI[i]].x();
                A[i+nCoeffs][xyzIndex] = fromPointsAll
                    [fromPointsProcIDsPI[i]]
                    [fromPointsIDsPI[i]].x();
                xyzIndex++;
            }

            if(yTerm_)
            {
                A[xyzIndex][i+nCoeffs] = fromPointsAll
                    [fromPointsProcIDsPI[i]]
                    [fromPointsIDsPI[i]].y();
                A[i+nCoeffs][xyzIndex] = fromPointsAll
                    [fromPointsProcIDsPI[i]]
                    [fromPointsIDsPI[i]].y();
                xyzIndex++;
            }

            if(zTerm_)
            {
                A[xyzIndex][i+nCoeffs] = fromPointsAll
                    [fromPointsProcIDsPI[i]]
                    [fromPointsIDsPI[i]].z();
                A[i+nCoeffs][xyzIndex] = fromPointsAll
                    [fromPointsProcIDsPI[i]]
                    [fromPointsIDsPI[i]].z();
                //xyzIndex++;
            }
        }

    //......     Compute inverse... 
        Info<< "Inversing matrix... ";

        labelList pivotIndices(A.n());
        LUDecompose(A, pivotIndices);
        simpleMatrix<scalar> invA(A.n(), 0.0, 0.0);

        for(label k = 0; k < A.n(); k++)
        {

            scalarField sourceSol(A.n(),0.0);
            sourceSol[k] = 1.0;
            LUBacksubstitute(A, pivotIndices, sourceSol);
            for(label m = 0; m < A.n(); m++)
            {
                invA[m][k] = sourceSol[m];
            }
        }
        Info<< "done." << endl;

        for(label pI = 0; pI < toPoints_.size(); pI++)
        {
            label xyzIndex(1);

            if(xTerm_)
            {
                A.source()[xyzIndex] = toPoints_[pI].x();
                xyzIndex++;
            }

            if(yTerm_)
            {
                A.source()[xyzIndex] = toPoints_[pI].y();
                xyzIndex++;
            }

            if(zTerm_)
            {
                A.source()[xyzIndex] = toPoints_[pI].z();
                //xyzIndex++;
            }

            for(label i = 0; i < fromPointsIDsPI.size(); i++)
            {
                A.source()[i+nCoeffs] = RBF
                (
                    mag
                    (
                        toPoints_[pI]
                      - fromPointsAll
                            [fromPointsProcIDsPI[i]]
                            [fromPointsIDsPI[i]]
                    )
                );
            }

    //...... Multiply source with inverse... scalarField w = ...
            scalarField w(A.n(),0.0);
            for(label i = 0; i < A.n(); i++)
            {
                for(label j = 0; j < A.n(); j++)
                {
                    w[i] += A.source()[j]*invA[j][i];
                }
            }

            for(label i = 0; i < fromPointsWeights_[pI].size(); i++)
            {
                fromPointsWeights_[pI][i] = w[i+nCoeffs];
            }
        }

    }
    else
    {
        for(label pI = 0; pI < toPoints_.size(); pI++)
        {
            labelList fromPointsProcIDsPI = fromPointsProcIDs_[pI];
            labelList fromPointsIDsPI = fromPointsIDs_[pI];

            simpleMatrix<scalar> A(fromPointsIDsPI.size() + nCoeffs, 0.0, 0.0);
            for(label i = 0; i < fromPointsIDsPI.size(); i++)
            {
                for(label j = 0; j < fromPointsIDsPI.size(); j++)
                {
                    A[i+nCoeffs][j+nCoeffs] = RBF
                    (
                        mag
                        (
                            fromPointsAll
                                [fromPointsProcIDsPI[i]]
                                [fromPointsIDsPI[i]]
                          - fromPointsAll
                                [fromPointsProcIDsPI[j]]
                                [fromPointsIDsPI[j]]
                        )
                    );
                }

                A[0][i+nCoeffs] = 1.0;
                A[i+nCoeffs][0] = 1.0;
                A.source()[0] = 1.0;

                label xyzIndex(1);

                if(xTerm_)
                {
                    A[xyzIndex][i+nCoeffs] = fromPointsAll
                        [fromPointsProcIDsPI[i]]
                        [fromPointsIDsPI[i]].x();
                    A[i+nCoeffs][xyzIndex] = fromPointsAll
                        [fromPointsProcIDsPI[i]]
                        [fromPointsIDsPI[i]].x();
                    A.source()[xyzIndex] = toPoints_[pI].x();
                    xyzIndex++;
                }

                if(yTerm_)
                {
                    A[xyzIndex][i+nCoeffs] = fromPointsAll
                        [fromPointsProcIDsPI[i]]
                        [fromPointsIDsPI[i]].y();
                    A[i+nCoeffs][xyzIndex] = fromPointsAll
                        [fromPointsProcIDsPI[i]]
                        [fromPointsIDsPI[i]].y();
                    A.source()[xyzIndex] = toPoints_[pI].y();
                    xyzIndex++;
                }

                if(zTerm_)
                {
                    A[xyzIndex][i+nCoeffs] = fromPointsAll
                        [fromPointsProcIDsPI[i]]
                        [fromPointsIDsPI[i]].z();
                    A[i+nCoeffs][xyzIndex] = fromPointsAll
                        [fromPointsProcIDsPI[i]]
                        [fromPointsIDsPI[i]].z();
                    A.source()[xyzIndex] = toPoints_[pI].z();
                    //xyzIndex++;
                }

                A.source()[i+nCoeffs] = RBF
                (
                    mag
                    (
                        toPoints_[pI]
                      - fromPointsAll
                            [fromPointsProcIDsPI[i]]
                            [fromPointsIDsPI[i]]
                    )
                );
            }

            scalarField w = A.LUsolve();

            for(label i = 0; i < fromPointsWeights_[pI].size(); i++)
            {
                fromPointsWeights_[pI][i] = w[i+nCoeffs];
            }
        }
    }
    Info<< "done." << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
