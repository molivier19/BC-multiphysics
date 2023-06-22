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
    Foam::newDisplacementLaplacianFvMotionSolver

Description
    Mesh motion solvers for an "fvMesh". Based on solving the cell Centre
    Laplacian for the motion displacement.

SourceFiles
    newDisplacementLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/


#ifndef newDisplacementLaplacianFvMotionSolver_H
#define newDisplacementLaplacianFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class newDisplacementLaplacianFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class newDisplacementLaplacianFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Optionally read point-position field. Used only for position
        //  boundary conditions.
        mutable autoPtr<pointVectorField> pointLocation_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Frozen points (that are not on patches). -1 or points that are
        //  fixed to be at points0_ location
        label frozenPointsZone_;
        
        //- Maximum number of outer iterations
        int nCorrectors_;
        
        //- Convergence tolerance
        scalar convergenceTolerance_;

        //- Grad-based volume-to-point interpolation
//        Switch useGradBasedVPI_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        newDisplacementLaplacianFvMotionSolver
        (
            const newDisplacementLaplacianFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const newDisplacementLaplacianFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("newDisplacementLaplacian");


    // Constructors

        //- Construct from polyMesh and data stream
        newDisplacementLaplacianFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        newDisplacementLaplacianFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );



    // Destructor

        ~newDisplacementLaplacianFvMotionSolver();


    // Member Functions


        //- Return reference to the cell motion displacement field
        volVectorField& cellDisplacement()
        {
            return cellDisplacement_;
        }

        //- Return const reference to the cell motion displacement field
        const volVectorField& cellDisplacement() const
        {
            return cellDisplacement_;
        }

        //- Return reference to the diffusivity field
        motionDiffusivity& diffusivity();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
