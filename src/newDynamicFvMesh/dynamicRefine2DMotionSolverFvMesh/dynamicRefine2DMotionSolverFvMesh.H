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
    Foam::dynamicRefine2DMotionSolverFvMesh

Description
    A FvMesh with a bulit in 2D dynamic refinement method combined with 
    a motion solver algorithm. The algorithm works with an Hex 2D splitting 
    method.

    \verbatim
    dynamicFvMesh   dynamicRefine2DMotionSolverFvMesh;

    motionSolver newDisplacementLaplacian;

    newDisplacementLaplacianCoeffs
    {
        nCorrectors 10;
        convergenceTolerance 1e-10;
        diffusivity  quadratic inverseDistance 1(InterfaceFluid);
    }

    dynamicRefine2DMotionSolverFvMeshCoeffs{
        // How often to refine
        refineInterval  1;
        // Field to be refinement on
        field           alpha.water;
        // Refine field inbetween lower..upper
        lowerRefineLevel 0.001;
        upperRefineLevel 0.999;
        // If value < unrefineLevel (default=GREAT) unrefine
        //unrefineLevel   10;
        // Have slower than 2:1 refinement
        nBufferLayers   1;
        // Refine cells only up to maxRefinement levels
        maxRefinement   2;
        // Stop refinement if maxCells reached
        maxCells        200000;
        // Flux field and corresponding velocity field. Fluxes on changed
        // faces get recalculated by interpolating the velocity. Use 'none'
        // on surfaceScalarFields that do not need to be reinterpolated, use
        // NaN to detect use of mapped variable
        correctFluxes
        (
            (phi none)  //NaN)   //none)
            (nHatf none)   //none)
            (rho*phi none)   //none)
            (ghf none)  //NaN)   //none)
        );

    // Tag all wall patch as protected cell. Those cells won't be refined.
        protectWall true;
        protectedWall (fluidWall solidWall);
        distance 0.0045;

        // Write the refinement level as a volScalarField
        dumpLevel       true;
    }
    \endverbatim

SourceFiles
    dynamicRefine2DMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefine2DMotionSolverFvMesh_H
#define dynamicRefine2DMotionSolverFvMesh_H

#include "dynamicFvMesh.H"
#include "hexRef2D.H"
#include "PackedBoolList.H"
#include "Switch.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;                                                                             // Motion Solver

/*---------------------------------------------------------------------------*\
                      Class dynamicRefine2DMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefine2DMotionSolverFvMesh
    : public dynamicFvMesh
{

  protected:
    //- Mesh cutting engine
    hexRef2D meshCutter_;

    //- Dump cellLevel for postprocessing
    Switch dumpLevel_;

    //- Fluxes to map
    List<Pair<word>> correctFluxes_;

    //- Number of refinement/unrefinement steps done so far.
    label nRefinementIterations_;

    //- Protected cells (usually since not hexes)
    PackedBoolList protectedCell_;

  private:
    // Private Member Functions
    autoPtr<motionSolver> motionPtr_;                                                                             // Motion Solver

  protected:
    //- Count set/unset elements in packedlist.
    static label count(const PackedBoolList &, const unsigned int);

    //- Calculate cells that cannot be refined since would trigger
    //  refinement of protectedCell_ (since 2:1 refinement cascade)
    void calculateProtectedCells(PackedBoolList &unrefineableCell) const;

    //- Read the projection parameters from dictionary
    void readDict();

    //- Refine cells. Update mesh and fields.
    autoPtr<mapPolyMesh> refine(const labelList &);

    //- Unrefine cells. Gets passed in centre points of cells to combine.
    autoPtr<mapPolyMesh> unrefine(const labelList &);

    // Selection of cells to un/refine

    //- Calculates approximate value for refinement level so
    //  we don't go above maxCell
    scalar getRefineLevel(
        const label maxCells,
        const label maxRefinement,
        const scalar refineLevel,
        const scalarField &) const;

    //- Get per cell max of connected point
    scalarField maxPointField(const scalarField &) const;

    //- Get point min of connected cell
    scalarField minCellField(const volScalarField &) const;

    scalarField cellToPoint(const scalarField &vFld) const;

    scalarField error(
        const scalarField &fld,
        const scalar minLevel,
        const scalar maxLevel) const;

    //- Select candidate cells for refinement
    virtual void selectRefineCandidates(
        const scalar lowerRefineLevel,
        const scalar upperRefineLevel,
        const scalarField &vFld,
        PackedBoolList &candidateCell) const;

    //- Subset candidate cells for refinement
    virtual labelList selectRefineCells(
        const label maxCells,
        const label maxRefinement,
        const PackedBoolList &candidateCell) const;

    //- Select edges that can be unrefined (2D case).
    virtual labelList selectUnrefineEdges(
        const scalar unrefineLevel,
        const PackedBoolList &markedCell,
        const scalarField &pFld) const;

    //- Extend markedCell with cell-face-cell.
    void extendMarkedCells(PackedBoolList &markedCell) const;

  private:
    //- Disallow default bitwise copy construct
    dynamicRefine2DMotionSolverFvMesh(const dynamicRefine2DMotionSolverFvMesh &);

    //- Disallow default bitwise assignment
    void operator=(const dynamicRefine2DMotionSolverFvMesh &);

  public:
    //- Runtime type information
    TypeName("dynamicRefine2DMotionSolverFvMesh");

    // Constructors

    //- Construct from IOobject
    explicit dynamicRefine2DMotionSolverFvMesh(const IOobject &io);

    //- Construct from components without boundary.
    //  Boundary is added using addFvPatches() member function
    explicit dynamicRefine2DMotionSolverFvMesh (
                                                    const IOobject& io,
                                                    pointField&& points,
                                                    faceList&& faces,
                                                    labelList&& allOwner,
                                                    labelList&& allNeighbour,
                                                    const bool syncPar = true
                                                );

    // Destructor

    virtual ~dynamicRefine2DMotionSolverFvMesh();

    // Member Functions

    //- Direct access to the refinement engine
    const hexRef2D &meshCutter() const
    {
        return meshCutter_;
    }

    //- Cells which should not be refined/unrefined
    const PackedBoolList &protectedCell() const
    {
        return protectedCell_;
    }

    //- Cells which should not be refined/unrefined
    PackedBoolList &protectedCell()
    {
        return protectedCell_;
    }

    //- Update the mesh for both mesh motion and topology change
    virtual bool update();

    // Writing

    //- Write using given format, version and compression
    virtual bool writeObject(
        IOstreamOption streamOpt,
        const bool valid) const;

    //- Return the motionSolver
        const motionSolver& motion() const;                                                                             // Motion Solver

    //- Combine Field
    template<class Type>
        void combineFields(Field<Type>& field)
        {
            List<Field<Type> > allValues(Pstream::nProcs());

            allValues[Pstream::myProcNo()] = field;

            Pstream::gatherList(allValues);
            Pstream::scatterList(allValues);

            field =
            ListListOps::combine<Field<Type> >
            (
                allValues,
                accessOp<Field<Type> >()
            );
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
