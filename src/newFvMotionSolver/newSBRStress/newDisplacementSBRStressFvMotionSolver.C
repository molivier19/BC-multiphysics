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

#include "newDisplacementSBRStressFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "mapPolyMesh.H"
//#include "fvMatrices.H"
//#include "gradBasedVolPointInterpolation.H"
//#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(newDisplacementSBRStressFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        newDisplacementSBRStressFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        newDisplacementSBRStressFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newDisplacementSBRStressFvMotionSolver::
newDisplacementSBRStressFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement().dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement().boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    nCorrectors_(readInt(coeffDict().lookup("nCorrectors"))),
    convergenceTolerance_
    (
        readScalar(coeffDict().lookup("convergenceTolerance"))
    )
//    useGradBasedVPI_(lookup("useGradBasedVPI"))
{}


Foam::newDisplacementSBRStressFvMotionSolver::
newDisplacementSBRStressFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            displacementMotionSolver::pointDisplacement().dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>
        (
            displacementMotionSolver::pointDisplacement().boundaryField()
        )
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    nCorrectors_(readInt(coeffDict().lookup("nCorrectors"))),
    convergenceTolerance_
    (
        readScalar(coeffDict().lookup("convergenceTolerance"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::newDisplacementSBRStressFvMotionSolver::
~newDisplacementSBRStressFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::newDisplacementSBRStressFvMotionSolver::curPoints() const
{
//    if(useGradBasedVPI_)
//    {
//        gradBasedVolPointInterpolation interpolator(fvMesh_);
//        pointDisplacement_ =
//            gradBasedVolPointInterpolation::New(fvMesh_).interpolate
//        (
//            cellDisplacement_,
//            fvc::grad(cellDisplacement_)()
//        );
//        Warning<< 
//            "gradBasedVolPointInterpolation not supported in parallel runs."
//            << endl;
//    }
//    else
//    {
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );
//    }

    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement().primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::newDisplacementSBRStressFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    const surfaceScalarField Df
    (
        dimensionedScalar("viscosity", dimViscosity, 1.0)
       *diffusivityPtr_->operator()()
    );

    volTensorField gradCd("gradCd", fvc::grad(cellDisplacement_));

    int iCorr = 0;
    scalar initialResidual = 0;
    scalar firstInitialResidual = 0;

    do
    {
        //cellDisplacement_.storePrevIter();
        //gradCd.storePrevIter();

        /*Foam::fvVectorMatrix dispEqn
        (
            fvm::laplacian
            (
                2*Df,
                cellDisplacement_,
                "laplacian(diffusivity,cellDisplacement)"
            )

          - fvc::laplacian
            (
                1*Df,
                cellDisplacement_,
                "laplacian(diffusivity,cellDisplacement)"
            )

          + fvc::div
            (
                Df
               *(
                   (
                       cellDisplacement_.mesh().Sf()
                     & fvc::interpolate(gradCd.T())
                   )

                   // Solid-body rotation "lambda" term
                 - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
               )
            )
        );*/

        Foam::fvVectorMatrix dispEqn
        (
            fvm::laplacian
            (
                2*Df,
                cellDisplacement_,
                "laplacian(diffusivity,cellDisplacement)"
            )

          + fvc::div
            (
                Df
               *(
                   fvc::dotInterpolate
                   (
                       cellDisplacement_.mesh().Sf(),
                       gradCd.T() - gradCd
                   )
                   // (
                   //     cellDisplacement_.mesh().Sf()
                   //   & fvc::interpolate(gradCd.T() - gradCd)
                   // )

                   // Solid-body rotation "lambda" term
                 - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
               )
            )
        );

        /*scalar nu = 1.0/3.0;
        Foam::fvVectorMatrix dispEqn
        (
            fvm::laplacian
            (
                (2.0 + 2.0*nu/(1 - 2.0*nu))*Df,
                cellDisplacement_,
                "laplacian(diffusivity,cellDisplacement)"
            )
            -fvc::laplacian
            (
                (1.0 + 2.0*nu/(1 - 2.0*nu))*Df,
                cellDisplacement_,
                "laplacianExp(diffusivity,cellDisplacement)"
            )

          + fvc::div
            (
                Df
               *(
                   (
                       cellDisplacement_.mesh().Sf()
                     & fvc::interpolate(gradCd.T())
                   )
                   // Solid-body rotation "lambda" term
                 + (2.0*nu/(1 - 2.0*nu))*cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
                 //- cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
                )
            )
        );*/
        initialResidual = cmptMax
        (
            dispEqn.solveSegregatedOrCoupled
            (
                dispEqn.solverDict()
            ).initialResidual()
        );
        if(iCorr == 0)
          {
            firstInitialResidual = initialResidual;
          }
        iCorr++;

        gradCd = fvc::grad(cellDisplacement_);
        //cellDisplacement_.relax();
        //gradCd.relax();
    }
    while (initialResidual > convergenceTolerance_ && iCorr < nCorrectors_);
    Info << endl << "Mesh Motion Solver Monitoring Summary" << endl
         << "cellDisplacement initial residual : " << firstInitialResidual << endl
         << "cellDisplacement final residual : " << initialResidual << endl
         << "No iteration : " << iCorr << endl << endl;
}


void Foam::newDisplacementSBRStressFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


// ************************************************************************* //
