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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Interpolate point field
template<class Type>
Field<Type>
radialExtrapolation::pointInterpolate
(
    const Field<Type>& pf
)
{
    if (pf.size() != size())
    {
        FatalErrorIn
        (
            "radialExtrapolation::pointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << fromPoints_.size() << " field size: " << pf.size()
            << abort(FatalError);
    }

    if(IDsNeedUpdate_)
    {
        generatePointIDs();
    }

    if(weightsNeedUpdate_)
    {
        fromPointsWeights_.clear();
        fromPointsWeights_.setSize(toPoints().size());

        for(int n = 0; n < toPoints_.size(); n++)
        {
            fromPointsWeights_[n].clear();
            fromPointsWeights_[n].setSize(fromPointsIDs_[n].size());

            for(label nn = 0; nn < fromPointsWeights_[n].size(); nn++)
            {
                fromPointsWeights_[n][nn] = 0;
            }
        }

        generatePointWeights();
        weightsNeedUpdate_ = false;
    }

    Field<Type> result
    (
        toPoints_.size(),
        pTraits<Type>::zero
    );

    // Assign values
    List< Field<Type> > pfAll(Pstream::nProcs());
    pfAll[Pstream::myProcNo()] = pf;
    Pstream::gatherList(pfAll);
    Pstream::scatterList(pfAll);

    /* forAll(result,toI) */
    for (label toI = 0; toI < result.size(); toI++)
    {
        result[toI] = pTraits<Type>::zero;

        /* forAll(fromPointsWeights_[toI], pI) */
        for (label pI = 0; pI < fromPointsWeights_[toI].size(); pI++)
        {
            result[toI] += fromPointsWeights_[toI][pI]
                *pfAll[fromPointsProcIDs_[toI][pI]][fromPointsIDs_[toI][pI]];
        }
    }

    return result;
}


template<class Type>
Field<Type>
radialExtrapolation::pointInterpolate
(
    const tmp<Field<Type> >& tpf
)
{
    tmp<Field<Type> > tint = pointInterpolate<Type>(tpf());
    tpf.clear();
    return tint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
