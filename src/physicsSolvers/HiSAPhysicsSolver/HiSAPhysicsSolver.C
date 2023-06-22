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

#include "HiSAPhysicsSolver.H"

#include "addToRunTimeSelectionTable.H"

#include "bound.H"
#include "preconditioner.H"
#include "solver.H"
#include "fvcSmooth.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "dynamicRefineFvMesh.H"
#include "dynamicRefine2DFvMesh.H"
#include "dynamicRefine3DMotionSolverFvMesh.H"
#include "dynamicRefine2DMotionSolverFvMesh.H"
#include "processorFvPatchField.H"
#include "wallDist.H"
#include "jacobian.H"
#include "characteristicWallPressureFvPatchScalarField.H"


#include "pointFields.H"
#include "motionSolver.H"
#include "multiFieldIterator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(HiSAPhysicsSolver, 0);

addToRunTimeSelectionTable(physicsSolver, HiSAPhysicsSolver, IOobject);

// * * * * * * * * * * * * Protected member functions  * * * * * * * * * * * //

void Foam::HiSAPhysicsSolver::setPseudoDeltaT()
{
    const fvMesh& mesh = mesh_;
    volScalarField& rPseudoDeltaT = rPseudoDeltaT_();
    {
        tmp< volScalarField > gamma = pThermo_->gamma();
        tmp< volScalarField > psi = pThermo_->psi();
        const volScalarField& ppsi = pThermo_->psi();
        tmp< surfaceScalarField > lambda;
        if (mesh.moving())
        {
            //lambda = fvc::interpolate(sqrt(gamma()/psi())) + mag((fvc::interpolate(U_())&mesh.Sf())-fvc::meshPhi(U_()))/mesh.magSf();
            lambda = fvc::interpolate(sqrt(gamma()/ppsi)) + mag(fvc::interpolate(U_())-fvc::meshPhi(U_())*mesh.Sf()/sqr(mesh.magSf()));
        }
        else
        {
            // lambda = fvc::interpolate(sqrt(gamma()/psi())) + mag(fvc::interpolate(U_())&mesh.Sf()/mesh.magSf());
            lambda = fvc::interpolate(sqrt(gamma()/ppsi)) + fvc::interpolate(mag(U_()));
        }

        if (localTimestepping_)
        {
            rPseudoDeltaT == dimensionedScalar("0", dimless/dimTime, 0.0);
            // surfaceScalarField frdt = mesh.deltaCoeffs()*lambda;
            surfaceScalarField frdt = mesh.nonOrthDeltaCoeffs()*(lambda());
            const labelList& own = mesh.owner();
            const labelList& nei = mesh.neighbour();
            forAll(own,facei)
            {
                rPseudoDeltaT[own[facei]] =
                    max(rPseudoDeltaT[own[facei]],frdt[facei]);
                rPseudoDeltaT[nei[facei]] =
                    max(rPseudoDeltaT[nei[facei]],frdt[facei]);
            }
            forAll(frdt.boundaryField(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                if (mesh.boundary()[patchi].coupled())
                {
                    forAll(fc, bfacei)
                    {
                        rPseudoDeltaT[fc[bfacei]] =
                            max
                            (
                                rPseudoDeltaT[fc[bfacei]],
                                frdt.boundaryField()[patchi][bfacei]
                            );
                    }
                }
                else if 
                (
                    mesh.boundary()[patchi].type() == "wall"
                )
                {
                    const scalarField pNonOrthDeltaCoef(mesh.nonOrthDeltaCoeffs().boundaryField()[patchi]);

                    forAll(fc, bfacei)
                    {
                        // Value based on cell value lambda and twice the patch edge length
                        scalar pLambda = 
                            0.5*pNonOrthDeltaCoef[bfacei]*
                            (
                                sqrt(gamma()[fc[bfacei]]/psi()[fc[bfacei]]) + 
                                mag(U_()[fc[bfacei]])
                            );
                    
                        rPseudoDeltaT[fc[bfacei]] =
                            max
                            (
                                rPseudoDeltaT[fc[bfacei]], // Cell value based on max of internal neighbours
                                pLambda 
                            );
                    }
                }
            }
            gamma.clear();
            psi.clear();
            lambda.clear();

            rPseudoDeltaT /= pseudoCoField_();
            rPseudoDeltaT.correctBoundaryConditions();
        }
        else
        {
            rPseudoDeltaT = max(mesh.deltaCoeffs()*lambda); // Global timestep
            rPseudoDeltaT /= pseudoCoNum_();
        }
    }

    if (localTimestepping_)
    {
        scalar totCells = mesh.globalData().nTotalCells();
        Info<< "Pseudo Courant No: "
            << "Min: "
            << min(pseudoCoField_()).value()
            << " Mean: "
            << sum(pseudoCoField_()).value()/totCells
            << " Max: "
            << max(pseudoCoField_()).value() << endl;
        Info<< "Pseudo deltaT: "
            << "Min: "
            << 1.0/max(rPseudoDeltaT).value()
            << " Mean: "
            << sum(1.0/rPseudoDeltaT).value()/totCells
            << " Max: "
            << 1.0/min(rPseudoDeltaT).value() << endl;
    }
    else
    {
        Info << "Pseudo Courant No: " << pseudoCoNum_->value() << endl;
        Info << "Pseudo deltaT: " << 1.0/min(rPseudoDeltaT).value() << endl;
    }
}


void Foam::HiSAPhysicsSolver::setPseudoCoNum()
{
    // Switched Evolution Relaxation
    if (initRes_.valid())
    {
        if (!solnControl_->firstIter() && prevRes_.valid())
        {
            residualIO& initRes(initRes_());
            residualIO& prevRes(prevRes_());

            scalar normInit = sqrt(magSqr(initRes.getScalar(0)) +
                                magSqr(initRes.getScalar(1)) +
                                magSqr(initRes.getVector(0)));
            scalar normPrev = sqrt(magSqr(prevRes.getScalar(0)) +
                                magSqr(prevRes.getScalar(1)) +
                                magSqr(prevRes.getVector(0)));
            scalar coNumRatio = normPrev/normInit;
            coNumRatio = max(min(coNumRatio, pseudoCoNumMaxIncr_), pseudoCoNumMinDecr_);
            if (localTimestepping_)
            {
                pseudoCoField_() *= coNumRatio;
                pseudoCoField_() = max(min(pseudoCoField_(), pseudoCoNumMax_), pseudoCoNumMin_);
                pseudoCoField_().correctBoundaryConditions();
            }
            else
            {
                pseudoCoNum_() *= coNumRatio;
                pseudoCoNum_().value() = max(min(pseudoCoNum_().value(), pseudoCoNumMax_), pseudoCoNumMin_);
            }
        }
        prevRes_.clear();
        prevRes_.set(new residualIO(initRes_()));
    }
}

void Foam::HiSAPhysicsSolver::createPreconditioners
(
    PtrList<preconditioner<2,1>>& preconditioners,
    PtrList<jacobian>& jacobians,
    const dictionary& parentDict
)
{
    if (parentDict.found("preconditioner"))
    {
        word preconditionerName(parentDict.lookup("preconditioner"));
        label nj = jacobians.size();
        label np = preconditioners.size();

        if (parentDict.found(preconditionerName))
        {
            const dictionary& dict = parentDict.subDict(preconditionerName);

            // Create preconditioner-specific Jacobian if specified
            if (dict.found("inviscidJacobian"))
            {
                jacobians.append
                (
                    new jacobian
                    (
                        parentDict.subDict(preconditionerName),
                        mesh_,
                        scalarVars_[0],
                        vectorVars_[0],
                        scalarVars_[1],
                        pThermo_(),
                        U_(),
                        ddtCoeff_(),
                        inviscid_,
                        turbulence_
                    )
                );
                nj = jacobians.size();
            }

            // Recursively create further preconditioners if applicable
            createPreconditioners(preconditioners, jacobians, dict);
        }

        // Create this preconditioner; pass newly created Jacobian if created,
        // else previous one, and pointer to child preconditioner if one
        // was created otherwise NULL.
        preconditioners.append
        (
            preconditioner<2,1>::New
            (
                parentDict,
                jacobians[nj-1].matrix(),
                np > 0 ? preconditioners(np-1) : NULL
            )
        );
    }
}

template<class FieldType, class Type>
void Foam::HiSAPhysicsSolver::parallelSyncFields(const wordList& fields)
{
    forAll(fields, i)
    {
        FieldType& f = const_cast<FieldType&>(meshPtr_->lookupObject<FieldType>(fields[i]));
        forAll(f.boundaryField(), patchI)
        {
            if (isA<processorFvPatchField<Type>>(f.boundaryField()[patchI]))
            {
                f.boundaryFieldRef()[patchI].initEvaluate();
            }
        }
        forAll(f.boundaryField(), patchI)
        {
            if (isA<processorFvPatchField<Type>>(f.boundaryField()[patchI]))
            {
                f.boundaryFieldRef()[patchI].evaluate();
            }
        }
    }
}

void Foam::HiSAPhysicsSolver::redistributePar()
{
    fvMesh& mesh = mesh_;

    // Check load balance
    scalarList procLoad(Pstream::nProcs(), 0.0);
    procLoad[Pstream::myProcNo()] = mesh.nCells();

    reduce(procLoad, sumOp<List<scalar>>());

    scalar overallLoad = sum(procLoad);
    scalar averageLoad = overallLoad/double(Pstream::nProcs());

    bool balanced = true;
    for (int i = 0; i < Pstream::nProcs(); i++)
    {
        if (Foam::mag(procLoad[i] - averageLoad)/averageLoad > maxLoadImbalance_)
        {
            balanced = false;
        }
    }
    Info << "Max load imbalance: " << max(Foam::mag(procLoad-averageLoad)/averageLoad)*100.0 << " %" << endl;

    // Redistribute
    if (!balanced)
    {
        if (rebalance_)
        {
            Info << "Redistributing parallel decomposition" << endl;
            scalar timeBeforeDist = mesh.time().elapsedCpuTime();

            labelList decomposition(mesh.nCells(), 0);

            {
                IOdictionary decomposeParDict
                (
                    IOobject
                    (
                        "decomposeParDict",
                        mesh.time().system(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
                // For convenience, try parallel scotch
                if (word(decomposeParDict.lookup("method")) == "scotch")
                {
                    decomposeParDict.set("method", "ptscotch");
                }
                autoPtr<decompositionMethod> decomposer
                (
                    decompositionMethod::New
                    (
                        decomposeParDict
                    )
                );
                decomposition = decomposer().decompose(mesh, mesh.cellCentres());
            }

            // Demand-driven data is cleared by fvMeshDistribute, but turbulence model has
            // stored a reference to the wall distance field. We preserve the object
            // (prevent a dangling reference) by checking it out of the object registry
            // and checking back in afterwards. The wall distance field itself stays
            // checked in and is redistributed together with all the other fields.
            wallDist* pwd = 0;
            if (mesh.foundObject<wallDist>(wallDist::typeName))
            {
                pwd = &const_cast<wallDist&>(mesh.lookupObject<wallDist>(wallDist::typeName));
                pwd->release();
            }

            {
                Info << "Distributing mesh..." << endl;

                // Create mesh re-distribution engine
                const scalar tolDim = 1e-6*mesh.bounds().mag();
                fvMeshDistribute distributor(mesh, tolDim);

                // Re-distribute the mesh
                autoPtr<mapDistributePolyMesh> map = distributor.distribute(decomposition);

                if (isA<dynamicRefineFvMesh>(mesh) )
                {
                    Info << "hexRef8 refinement with dynamicRefineFvMesh" << endl;
                    // Update the refinement
                    dynamicRefineFvMesh& refineMesh = dynamic_cast<dynamicRefineFvMesh&>(mesh);
                    
                    const_cast<hexRef8&>(refineMesh.meshCutter()).distribute(map);
                    // Update protected-cell list

                    bitSet& pc = refineMesh.protectedCell();

                    List<bool> upc(pc.size());
                    forAll(pc, cellI)
                    {
                        upc[cellI] = pc[cellI];
                    }
                    map->distributeCellData(upc);

                    pc = bitSet(upc);

                }
                else if (isA<dynamicRefine3DMotionSolverFvMesh>(mesh))
                {
                    Info << "hexRef8 refinement with dynamicRefine3DMotionSolverFvMesh" << endl;
                    // Update the refinement
                    dynamicRefine3DMotionSolverFvMesh& refineMesh = dynamic_cast<dynamicRefine3DMotionSolverFvMesh&>(mesh);

                    const_cast<hexRef8&>(refineMesh.meshCutter()).distribute(map);
                    // Update protected-cell list

                    bitSet& pc = refineMesh.protectedCell();

                    List<bool> upc(pc.size());
                    forAll(pc, cellI)
                    {
                        upc[cellI] = pc[cellI];
                    }
                    map->distributeCellData(upc);

                    pc = bitSet(upc);

                }
                else if (isA<dynamicRefine2DFvMesh>(mesh) )
                {
                    Info << "hexRef2D refinement with dynamicRefine2DFvMesh" << endl;
                    // Update the refinement
                    dynamicRefine2DFvMesh& refineMesh = dynamic_cast<dynamicRefine2DFvMesh&>(mesh);
                    
                    const_cast<hexRef2D&>(refineMesh.meshCutter()).distribute(map);
                    // Update protected-cell list

                    PackedBoolList& pc = refineMesh.protectedCell();

                    List<bool> upc(pc.size());
                    forAll(pc, cellI)
                    {
                        upc[cellI] = pc[cellI];
                    }
                    map->distributeCellData(upc);

                    pc = PackedBoolList(upc);

                }
                else if (isA<dynamicRefine2DMotionSolverFvMesh>(mesh))
                {
                    Info << "hexRef2D refinement with dynamicRefine2DMotionSolverFvMesh" << endl;
                    // Update the refinement
                    dynamicRefine2DMotionSolverFvMesh& refineMesh = dynamic_cast<dynamicRefine2DMotionSolverFvMesh&>(mesh);
                    
                    const_cast<hexRef2D&>(refineMesh.meshCutter()).distribute(map);
                    // Update protected-cell list

                    PackedBoolList& pc = refineMesh.protectedCell();

                    List<bool> upc(pc.size());
                    forAll(pc, cellI)
                    {
                        upc[cellI] = pc[cellI];
                    }
                    map->distributeCellData(upc);

                    pc = PackedBoolList(upc);

                }
            }

            if (pwd)
            {
                pwd->store();
            }

            Info << "Parallel distribution..." << endl;
            // The distributor zeros the processor patches, so re-sync

            const wordList psfs(mesh.names(pointScalarField::typeName)); // Adding pointField Compatibility for redistribution
            const wordList pvfs(mesh.names(pointVectorField::typeName)); // Adding pointField Compatibility for redistribution

            const wordList vsfs(mesh.names(volScalarField::typeName));
            const wordList vvfs(mesh.names(volVectorField::typeName));
            const wordList vtfs(mesh.names(volTensorField::typeName));
            const wordList vstfs(mesh.names(volSymmTensorField::typeName));
            const wordList vsphtfs(mesh.names(volSphericalTensorField::typeName));

            parallelSyncFields<pointScalarField, scalar>(psfs); // Adding pointField Compatibility for redistribution
            parallelSyncFields<pointVectorField, vector>(pvfs); // Adding pointField Compatibility for redistribution

            parallelSyncFields<volScalarField, scalar>(vsfs);
            parallelSyncFields<volVectorField, vector>(vvfs);
            parallelSyncFields<volTensorField, tensor>(vtfs);
            parallelSyncFields<volSymmTensorField, symmTensor>(vstfs);
            parallelSyncFields<volSphericalTensorField, sphericalTensor>(vsphtfs);

            scalarList procLoadNew (Pstream::nProcs(), 0.0);
            procLoadNew[Pstream::myProcNo()] = mesh.nCells();

            reduce(procLoadNew, sumOp<List<scalar> >());

            scalar overallLoadNew = sum(procLoadNew);
            scalar averageLoadNew = overallLoadNew/double(Pstream::nProcs());

            Info << "Max load imbalance: " << max(Foam::mag(procLoadNew-averageLoadNew)/averageLoadNew)*100.0 << " %" << endl;

            Info<< "Execution time for mesh redistribution = "
                << mesh.time().elapsedCpuTime() - timeBeforeDist << " s\n" << endl;
        }
        else
        {
            Info << "Parallel mesh redistribution not selected" << endl;
        }
    }
}

void Foam::HiSAPhysicsSolver::findDebugCell()
{
    // Code adapted from findRefCell.C
    const dictionary& dict = meshPtr_->time().controlDict();
    if (dict.found("debugCell"))
    {
        cellDebugging_ = true;

        if (!Pstream::parRun() || Pstream::myProcNo() == readLabel(dict.lookup("debugProcessor")))
        {
            debugCell_ = readLabel(dict.lookup("debugCell"));
            if (debugCell_ < 0 || debugCell_ >= meshPtr_->nCells())
            {
                debugCell_ = -1;
                Pout << "Debug cell " << debugCell_ << " is out of range." << endl;
            }
        }
        else
        {
            debugCell_ = -1;
        }

        label hasRef = (debugCell_ >= 0 ? 1 : 0);
        if (returnReduce<label>(hasRef, sumOp<label>()) != 1)
        {
            WarningInFunction
                << "debugCell not found." << endl;
            cellDebugging_ = false;
        }
    }
    else if (dict.found("debugPoint"))
    {
        point debugPoint(dict.lookup("debugPoint"));

        // Try fast approximate search avoiding octree construction
        debugCell_ = meshPtr_->findCell(debugPoint, polyMesh::FACE_PLANES);

        label hasRef = (debugCell_ >= 0 ? 1 : 0);
        label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

        // If reference cell no found use octree search
        // with cell tet-decompositoin
        if (sumHasRef != 1)
        {
            debugCell_ = meshPtr_->findCell(debugPoint);

            hasRef = (debugCell_ >= 0 ? 1 : 0);
            sumHasRef = returnReduce<label>(hasRef, sumOp<label>());
        }

        if (sumHasRef != 1)
        {
            WarningInFunction
                << "Unable to set debug cell at point " << debugPoint << ":"
                << nl << "Found on " << sumHasRef << " domains (should be one)"
                << endl;
            cellDebugging_ = false;
        }
        else
        {
            cellDebugging_ = true;
        }
    }
    else
    {
        cellDebugging_ = false;
    }
    if (cellDebugging_)
    {
        if (debugCell_ >= 0)
        {
            Pout << "Found debug cell " << debugCell_ << nl << endl;
        }
    }
}

void Foam::HiSAPhysicsSolver::cellDebug()
{
    // ================================================ PtrInit ===============================================================

    dynamicFvMesh& mesh = mesh_;
    // const Time& runTime = mesh.time();
    // pseudotimeControl& solnControl = solnControl_();

    psiThermo& thermo = pThermo_();

    // volScalarField& p = thermo.p();
    // volScalarField& e = thermo.he();
    // volScalarField& T = thermo.T();
    // const volScalarField& psi = thermo.psi();

    surfaceVectorField& phiUp = phiUp_();
    surfaceScalarField& phiEp = phiEp_();

    volVectorField& U = U_();

    volScalarField& rho = scalarVars_[0];
    volVectorField& rhoU = vectorVars_[0];
    volScalarField& rhoE = scalarVars_[1];

    surfaceScalarField& phi = phi_();

    volScalarField& rPseudoDeltaT = rPseudoDeltaT_();
    // scalarField& ddtCoeff = ddtCoeff_();

    volScalarField& rhoR = scalarResiduals_[0];
    volVectorField& rhoUR = vectorResiduals_[0];
    volScalarField& rhoER = scalarResiduals_[1];

    // fluxScheme& flux = flux_();

    // ================================================ PtrInit (end) ===============================================================

    if (cellDebugging_)
    {
        volScalarField ee = (rhoE-0.5*magSqr(rhoU)/rho)/rho;
        surfaceScalarField ef = fvc::interpolate(ee, phi, "reconstruct(T)");
        surfaceScalarField phie = phi*ef;
        surfaceVectorField Uf = Up_()/mesh.magSf();
        surfaceScalarField phiEk = phi*0.5*magSqr(Uf);
        tmp<surfaceVectorField> phiUSnGradMuU;
        tmp<surfaceVectorField> phiUTauMC;
        if (!inviscid_)
        {
            volScalarField muEff("muEff", turbulence_->muEff());
            phiUSnGradMuU = -fvc::interpolate(muEff)*fvc::snGrad(U)*mesh.magSf();
            phiUTauMC = -mesh.Sf() & fvc::interpolate(tauMC_());
        }

        if (debugCell_ > -1)
        {
            Pout << "Debug information at cell " << debugCell_ << endl;
            scalar rhoMagUSqr = magSqr(rhoU[debugCell_])/rho[debugCell_];
            Pout<< "rho:" << tab << rho[debugCell_] << nl
                << "rhoU:" << tab << rhoU[debugCell_] << nl
                << "rhoE:" << tab << rhoE[debugCell_] << nl
                << "U:" << tab << rhoU[debugCell_]/rho[debugCell_] << nl
                << "p:" << tab << thermo.p()[debugCell_] << nl
                << "0.5*rho|U|^2:" << tab << 0.5*rhoMagUSqr << nl
                << "rhoe:" << tab << rhoE[debugCell_]-0.5*rhoMagUSqr << nl
                << "T:" << tab << thermo.T()[debugCell_] << nl
                << "Volume:" << tab << mesh.V()[debugCell_] << nl
                << "deltaT:" << tab << 1.0/(rPseudoDeltaT[debugCell_]+SMALL) << endl;
            if (localTimestepping_)
            {
                Pout << "pseudoCo: " << tab << pseudoCoField_()[debugCell_] << endl;
            }
            if (mesh.foundObject<volScalarField>("k"))
            {
                Pout << "k:" << tab << mesh.lookupObject<volScalarField>("k")[debugCell_] << endl;
            }
            if (mesh.foundObject<volScalarField>("omega"))
            {
                Pout << "omega:" << tab << mesh.lookupObject<volScalarField>("omega")[debugCell_] << endl;
            }
            Pout << "Residuals:" << nl;
            Pout << rhoR[debugCell_] << " " << rhoUR[debugCell_] << " " << rhoER[debugCell_] << endl;
            Pout << "Neighbouring faces:" << endl;
            forAll(phi, faceI)
            {
                label own = 0;
                label otherCellI = -1;
                if (mesh.owner()[faceI] == debugCell_)
                {
                    own = 1;
                    otherCellI = mesh.neighbour()[faceI];
                }
                else if (mesh.neighbour()[faceI] == debugCell_)
                {
                    own = -1;
                    otherCellI = mesh.owner()[faceI];
                }
                if (own != 0)
                {
                    Pout<< "face:" << tab << faceI << nl
                        << "\tSf:" << tab << mesh.Sf()[faceI]*own << nl
                        << "\tphi:" << tab << phi[faceI]*own << nl
                        << "\tphiUp:" << tab << phiUp[faceI]*own << nl
                        << "\tphiEp:" << tab << phiEp[faceI]*own << nl
                        << "\tphie:" << tab << phie[faceI]*own << nl
                        << "\tphiEk:" << tab << phiEk[faceI]*own << nl
                        << "\tUf:" << tab << Uf[faceI] << nl
                        << "\tFace ratio:" << tab << (own < 0 ? 1 : 0) + own*mesh.weights()[faceI] << endl;
                    if (!inviscid_)
                    {
                        Pout << "\tphiUSnGradMuU:" << tab << phiUSnGradMuU()[faceI]*own << endl;
                        Pout << "\tphiUtauMC:" << tab << phiUTauMC()[faceI]*own << endl;
                    }
                    scalar rhoMagUSqr = magSqr(rhoU[otherCellI])/rho[otherCellI];
                    Pout<< " neighbouring cell:" << tab << otherCellI << nl
                        << "\trho:" << tab << rho[otherCellI] << nl
                        << "\trhoU:" << tab << rhoU[otherCellI] << nl
                        << "\trhoE:" << tab << rhoE[otherCellI] << nl
                        << "\tU:" << tab << rhoU[otherCellI]/rho[otherCellI] << nl
                        << "\tp:" << tab << thermo.p()[otherCellI] << nl
                        << "\t0.5*rho|U|^2:" << tab << 0.5*rhoMagUSqr << nl
                        << "\trhoe:" << tab << rhoE[otherCellI]-0.5*rhoMagUSqr << nl
                        << "\tT:" << tab << thermo.T()[otherCellI] << nl
                        << "\tVolume:" << tab << mesh.V()[otherCellI] << nl
                        << "\tdeltaT:" << tab << 1.0/(rPseudoDeltaT[otherCellI]+SMALL) << endl;
                    if (localTimestepping_)
                    {
                        Pout << "\tpseudoCo: " << tab << pseudoCoField_()[otherCellI] << endl;
                    }
                    if (mesh.foundObject<volScalarField>("k"))
                    {
                        Pout << "\tk:" << tab << mesh.lookupObject<volScalarField>("k")[otherCellI] << endl;
                    }
                    if (mesh.foundObject<volScalarField>("omega"))
                    {
                        Pout << "\tomega:" << tab << mesh.lookupObject<volScalarField>("omega")[otherCellI] << endl;
                    }
                }
            }
            forAll(phi.boundaryField(), patchI)
            {
                forAll(phi.boundaryField()[patchI], bfaceI)
                {
                    if (mesh.boundary()[patchI].faceCells()[bfaceI] == debugCell_)
                    {
                        Pout<< "patch:" << tab << patchI << tab
                            << "bface:" << tab << bfaceI << nl
                            << "\tSf:" << tab << mesh.Sf().boundaryField()[patchI][bfaceI] << nl
                            << "\tphi:" << tab << phi.boundaryField()[patchI][bfaceI] << nl
                            << "\tphiUp:" << tab << phiUp.boundaryField()[patchI][bfaceI] << nl
                            << "\tphiEp:" << tab << phiEp.boundaryField()[patchI][bfaceI] << nl
                            << "\tU_n:" << tab <<
                            (
                                U[debugCell_] & mesh.Sf().boundaryField()[patchI][bfaceI]
                            )/mesh.magSf().boundaryField()[patchI][bfaceI] << endl;
                    }
                }
            }
        }
    }
}

/* TODO - Check if interpolation function is related to physicsSolver or rhoCentralFoam.
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::HiSAPhysicsSolver::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& dir,
    const word& reconFieldName
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        fvc::interpolate
        (
            vf,
            dir,
            "reconstruct("
          + (reconFieldName != word::null ? reconFieldName : vf.name())
          + ')'
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    sf.rename(vf.name() + '_' + dir.name());

    return tsf;
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HiSAPhysicsSolver::HiSAPhysicsSolver
(
    const IOobject& io, 
    const word& regionName
)
:
    physicsSolver(io, regionName),

    meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region_,
                io.time().timeName(),
                io.time(),
                IOobject::MUST_READ
            )
        )
    ),

    mesh_(meshPtr_()),

    solnControl_(nullptr),

    pThermo_(nullptr),
    
    U_(nullptr),
    tauMC_(nullptr),
    phi_(nullptr),
    phiUp_(nullptr),
    phiEp_(nullptr),
    Up_(nullptr),
    eZero_(nullptr),
    rPseudoDeltaT_(nullptr),
    finished_(true),
    ddtCoeff_(nullptr),
    
    bounded_(false),
    
    scalarVars_(PtrList<volScalarField>()),
    scalarVarsPrevIter_(PtrList<volScalarField>()),
    vectorVars_(PtrList<volVectorField>()),
    vectorVarsPrevIter_(PtrList<volVectorField>()),
    scalarResiduals_(PtrList<volScalarField>()),
    vectorResiduals_(PtrList<volVectorField>()),

    sigmaFluidptr_(nullptr),
    
    turbulence_(nullptr),
    
    flux_(nullptr),
    
    inviscid_(true),

    moveMeshOuterCorrectors_(0),

    steadyState_(true),
    localTimestepping_(true),
    localTimesteppingBounding_(true),
    localTimesteppingLowerBound_(0.0),
    localTimesteppingUpperBound_(0.0),

    initRes_(nullptr),
    prevRes_(nullptr),

    pseudoCoNum_(nullptr), 
    pseudoCoField_(nullptr),
    
    pseudoCoNumMin_(0.0),
    pseudoCoNumMax_(1.0),
    pseudoCoNumMaxIncr_(0.1),
    pseudoCoNumMinDecr_(0.1),

    cellDebugging_(false),
    debugCell_(0),

    rebalance_(false),
    maxLoadImbalance_(0.1)
{
    // ================================== Initialise Simulation parameter ==================================
    
    dynamicFvMesh& mesh = mesh_;
    const Time& runTime = mesh.time();

    // Detect steady-state analysis
    const dictionary& ddtControls = mesh.schemesDict().subDict("ddtSchemes");
    wordList ddtToc (ddtControls.toc());
    steadyState_ = false;
    forAll(ddtToc,s)
    {
        word ddtScheme(ddtToc[s]);
        word ddtSchemeLastWord;
        const tokenList& tokens = ddtControls.lookup(ddtToc[s]);
        if (tokens.last().isWord() && tokens.last().wordToken() == "steadyState")
        {
            if (ddtToc[s] == "default" || ddtToc[s] == "rhoU")
            {
                steadyState_ = true;
            }
        }
    }
    if (steadyState_)
    {
        Info << "Steady-state analysis detected" << nl << endl;
    }
    else
    {
        Info << "Transient analysis detected" << nl << endl;
    }

    residualIO defaultPseudoTol(2, 1, residualOrdering, 1e-4);
    residualIO defaultPseudoTolRel(2, 1, residualOrdering, 0.0);

    solnControl_.set
    (
        new pseudotimeControl
        (
            mesh,
            steadyState_,
            2,
            1,
            defaultPseudoTol,
            defaultPseudoTolRel
        )
    );
    if (steadyState_)
    {
        solnControl_->setCorr(runTime.startTimeIndex());
    }

    localTimestepping_ =
        solnControl_->dict().lookupOrDefault<Switch>
        (
            "localTimestepping",
            true
        );
    if (localTimestepping_)
    {
        Info << "Local timestepping selected" << nl << endl;
        localTimesteppingBounding_ =
            solnControl_->dict().lookupOrDefault<Switch>
            (
                "localTimesteppingBounding",
                true
            );

        localTimesteppingLowerBound_ =
            solnControl_->dict().lookupOrDefault<scalar>
            (
                "localTimesteppingLowerBound",
                0.95
            );
        localTimesteppingLowerBound_ =
            (localTimesteppingLowerBound_ > 0 ? localTimesteppingLowerBound_ : 0.0);
        localTimesteppingLowerBound_ =
            (localTimesteppingLowerBound_ < 0.99 ? localTimesteppingLowerBound_ : 0.99);

        // localTimesteppingUpperBound_ =
        //     solnControl_->dict().lookupOrDefault<scalar>
        //     (
        //         "localTimesteppingUpperBound",
        //         1.5
        //     );
        // localTimesteppingUpperBound_ =
        //     (localTimesteppingUpperBound_ < 1.01 ? localTimesteppingUpperBound_ : 1.01);
    }
    else
    {
        Info << "Global timestepping selected" << nl << endl;
        localTimesteppingBounding_ = false;
    }

    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================

    // Read or initialise pseudo Co number
    // NOTE: It is not necessarily read before every outer iteration
    // (see resetPseudo)
    if (!localTimestepping_)
    {
        pseudoCoNum_.set
        (
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "pseudoCoNum",
                    runTime.timeName(),
                    "uniform",
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                dimensionedScalar
                (
                    "pseudoCoNumInit",
                    dimless,
                    solnControl_().dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
                )
             )
        );

        Info << "Initial pseudo Courant No: " << pseudoCoNum_->value() << nl << endl;
    }
    else
    {
        pseudoCoField_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "pseudoCoField",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "pseudoCoNumInit",
                    dimless,
                    solnControl_().dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
                )
            )
        );
        pseudoCoField_();

        scalar totCells = mesh.globalData().nTotalCells();
        Info<< "Initial pseudo Courant No: "
            << "Min: "
            << min(pseudoCoField_()).value()
            << " Mean: "
            << sum(pseudoCoField_()).value()/totCells
            << " Max: "
            << max(pseudoCoField_()).value() << nl << endl;
    }

    findDebugCell();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HiSAPhysicsSolver::~HiSAPhysicsSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HiSAPhysicsSolver::createFields()
{
    dynamicFvMesh& mesh = mesh_;
    const Time& runTime = mesh.time();

    Info<< "Reading thermophysical properties\n" << endl;
    pThermo_ =
    (
        psiThermo::New(mesh)
    );
    Info << endl;
    psiThermo& thermo = pThermo_();

    // volScalarField& p = thermo.p();
    volScalarField& e = thermo.he();
    if (e.name() != "e")
    {
        FatalErrorInFunction
            << "Only energy type internalEnergy supported."
            << nl << exit(FatalError);
    }
    const volScalarField& mu = thermo.mu();

    if (gMax(mu.internalField()) > 0.0)
    {
        inviscid_ = false;
        Info << "Viscous analysis detected" << nl << endl;
    }
    else
    {
        inviscid_ = true;
        Info << "Inviscid analysis detected" << nl << endl;
    }

    Info<< "Reading field U" << nl << endl;
    U_.set
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    volVectorField& U = U_();

    scalarVars_.setSize(2);
    vectorVars_.setSize(1);

    /*// =================================== rhoBoundaryType ===================================
    const volScalarField::Boundary& pbf = p.boundaryField();
    wordList rhoBoundaryTypes = pbf.types();

    forAll(rhoBoundaryTypes, patchi)
    {
        if (rhoBoundaryTypes[patchi] == "waveTransmissive")
        {
            rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
        }
        else if (rhoBoundaryTypes[patchi] == "characteristicWallPressure")
        {
            rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
        }
        // For characteristic-based calculations fixedRho is used which computes rho
        // directly from p & T
        else if (pbf[patchi].fixesValue() || (rhoBoundaryTypes[patchi] == "characteristicFarfieldPressure") )
        {
            rhoBoundaryTypes[patchi] = fixedRhoFvPatchScalarField::typeName; // NOTE: Based on current values (rho = p * psi)
        }
    }
    // =================================== rhoBoundaryType(end) ===================================*/

    scalarVars_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            thermo.rho(),
            calculatedFvPatchScalarField::typeName
            //rhoBoundaryTypes // NOTE: fixRho does a psi*p update (Check if working with nozzle BCs)
        )
    );
    volScalarField& rho = scalarVars_[0];

    vectorVars_.set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "rhoU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            rho*U
        )
    );
    volVectorField& rhoU = vectorVars_[0];

    volScalarField TZero
    (
        IOobject("TZero", runTime.timeName(), mesh), 
        mesh, 
        dimensionedScalar("0", thermo.T().dimensions(), SMALL)
    );
    volScalarField pZero
    (
        IOobject("pZero", runTime.timeName(), mesh), 
        mesh, 
        dimensionedScalar("0", thermo.p().dimensions(), SMALL)
    );
    eZero_.set(new volScalarField("eZero", thermo.he(pZero, TZero)));

    scalarVars_.set
    (
        1,
        new volScalarField
        (
            IOobject
            (
                "rhoE",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            rho*(e - eZero_() + 0.5*magSqr(U))
        )
    );
    volScalarField& rhoE = scalarVars_[1];

    phi_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh.Sf() & fvc::interpolate(rhoU)
        )
    );
    surfaceScalarField& phi = phi_();


    Info<< "Creating turbulence model\n" << endl;
    turbulence_ =
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );
    Info << endl;

    Info<< "Creating Sigma Fluid Field\n" << endl;
    sigmaFluidptr_.reset
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigmaFluid",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - thermo.rho()*(
          thermo.p()*symmTensor(1,0,0,1,0,1)
          + turbulence_->devRhoReff())

        )
    );

    Info<< "Creating flux scheme\n" << endl;
    flux_ =
    (
        fluxScheme::New(mesh.schemesDict(), thermo, rho, U, rhoU, rhoE)
    );
    Info << endl;

    rPseudoDeltaT_.set
    (
        new volScalarField
        (
            IOobject
            (
                "rPseudoDeltaT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("rPseudoDeltaT", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    ddtCoeff_.set(new scalarField);

    scalarResiduals_.setSize(2);
    vectorResiduals_.setSize(1);

    // Continuity residual
    scalarResiduals_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "rhoR",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", rho.dimensions()/dimTime, 0.0)
        )
    );


    // Mom Residual
    vectorResiduals_.set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "rhoUR",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", rhoU.dimensions()/dimTime, vector::zero)
        )
    );

    // Energy residual
    scalarResiduals_.set
    (
        1,
        new volScalarField
        (
            IOobject
            (
                "rhoER",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", rhoE.dimensions()/dimTime, 0.0)
        )
    );

    if (!inviscid_)
    {
        // tauMC is stored for the benefit of the maxwellSlip BC - maybe better
        // to modify the BC to recalculate it rather
        tauMC_.set
        (
            new volTensorField
            (
                IOobject
                (
                    "tauMC",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                turbulence_->muEff()*dev2(Foam::T(fvc::grad(U)))
            )
        );
    }

    //rho.correctBoundaryConditions();
}


bool Foam::HiSAPhysicsSolver::isConverged()
{
    // The implementation of the pseudo timeStep manage the convergence criteria. 
    // In the original solver the convergence is check before solving the domain. 
    return solnControl_().loop();
}


void Foam::HiSAPhysicsSolver::solveFields()
{
    // ================================================ PtrInit ===============================================================
    dynamicFvMesh& mesh = mesh_;
    const Time& runTime = mesh.time();
    pseudotimeControl& solnControl = solnControl_();

    psiThermo& thermo = pThermo_();

    volScalarField& p = thermo.p();
    volScalarField& e = thermo.he();
    volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

    volSymmTensorField& sigmaFluid = sigmaFluidptr_();

    volVectorField& U = U_();

    volScalarField& rho = scalarVars_[0];
    volVectorField& rhoU = vectorVars_[0];
    volScalarField& rhoE = scalarVars_[1];

    surfaceScalarField& phi = phi_();

    // volScalarField& rPseudoDeltaT = rPseudoDeltaT_();
    scalarField& ddtCoeff = ddtCoeff_();

    volScalarField& rhoR = scalarResiduals_[0];
    volVectorField& rhoUR = vectorResiduals_[0];
    volScalarField& rhoER = scalarResiduals_[1];

    fluxScheme& flux = flux_();
    // ================================================ PtrInit (end) ===============================================================
    
    // ========================================== solnControl Adjustement ===========================================================
    
    // Get Iterator value from multiFieldSolver
    const multiFieldIterator& mfi
    (
        mesh.time().lookupObject<multiFieldIterator>
        (
            "iterator"
        )
    );

    // Assign multiField iterator value to pseudoTimeControl
    solnControl.setCorr(mfi.index());

    // Initialise pseudoTimeControl if first iteration
    if (solnControl.firstIter())
    {
        solnControl.init();
    }

    // ========================================== solnControl Adjustement (end) ========================================================

    if 
    (
        solnControl.firstIter() 
    || steadyState_ 
    || (moveMeshOuterCorrectors_ && !(solnControl.corr() % moveMeshOuterCorrectors_))
    )
    {
        if (isA<dynamicFvMesh>(mesh))
        {
            // Do any mesh changes
            
            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

            // Work around issues with mapping of oriented fields
            surfaceVectorField Un(phi*mesh.Sf()/sqr(mesh.magSf()));
            #if OPENFOAM >= 1712
            phi.setOriented(false);
            #endif

            dynamicCast<dynamicFvMesh&>(mesh).update();

            if (mesh.changing())
            {
                Info<< "Execution time for mesh.update() = "
                    << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                    << " s" << endl;
            }
            if (mesh.topoChanging())
            {
                redistributePar();
            }

            #if OPENFOAM >= 1712
            phi.setOriented(true);
            #endif
            phi = (Un & mesh.Sf());
        }
    }

    // Store value before solve for bounding. Don't use storePrevIter here
    // as it introduces too many problems for mesh refinement/redistribution
    volScalarField rhoPrevIter("rhoPrevIter", rho);
    bounded_ = false;

    phiUp_.set
    (
        new surfaceVectorField
        (
            IOobject
            (
                "phiUp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*rhoU.dimensions()*U.dimensions()
        )
    );
    phiEp_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiEp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*rhoE.dimensions()*U.dimensions()
        )
    );
    Up_.set
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Up",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*U.dimensions()
        )
    );

    surfaceVectorField& phiUp = phiUp_();
    surfaceScalarField& phiEp = phiEp_();
    surfaceVectorField& Up = Up_();

    // Interpolation of primitive fields on faces
    flux.calcFlux(phi, phiUp, phiEp, Up);

    // =================================== residualsUpdate ===================================
    
    // Construct cont eqn: [A] \Delta W = [R]
    rhoR =
    (
        -fvc::div(phi)
    );

    // Construct mom eqn: [A] \Delta W = [R]
    rhoUR =
    (
        -fvc::div(phiUp)
    );

    // Construct energ eqn: [A] \Delta W = [R]
    rhoER =
    (
        -fvc::div(phiEp)
    );

    // Viscous flow
    if (!inviscid_)
    {
        volScalarField muEff("muEff", turbulence_->muEff());
        volScalarField alphaEff("alphaEff", turbulence_->alphaEff());

        volTensorField& tauMC = tauMC_();
        tauMC = muEff*dev2(Foam::T(fvc::grad(U)));
        rhoUR += fvc::laplacian(muEff,U);
        rhoUR += fvc::div(tauMC);

        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
            fvc::interpolate(muEff)*fvc::snGrad(U)
            + (mesh.Sf()/mesh.magSf() & fvc::interpolate(tauMC))
            )
            & Up
        );
        rhoER += fvc::div(sigmaDotU);

        volScalarField eCalc("eCalc", rhoE/rho-0.5*magSqr(U)); // prevent e inheriting BC from T in thermo
        rhoER += fvc::laplacian(alphaEff, eCalc, "laplacian(alphaEff,e)");
    }

    cellDebug();

    // Pseudo and real-time contribution

    tmp<fvScalarMatrix> ddtRho = fvm::ddt(rho);
    tmp<fvVectorMatrix> ddtRhoU = fvm::ddt(rhoU);
    tmp<fvScalarMatrix> ddtRhoE = fvm::ddt(rhoE);
    
#if OPENFOAM >= 1712
    // Work around bug causing timeIndex not to be stored for ddt0 field (CrankNicolson)
    // Seems to be caused by commented out line in GeometricField::storeOldTimes(),
    // commit 44ce956798
    if (mesh.foundObject<volScalarField>("ddt0(rho)"))
    {
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>("ddt0(rho)")
        ).timeIndex() = mesh.time().timeIndex();
    }
    if (mesh.foundObject<volVectorField>("ddt0(rhoU)"))
    {
        const_cast<volVectorField&>
        (
            mesh.lookupObject<volVectorField>("ddt0(rhoU)")
        ).timeIndex() = mesh.time().timeIndex();
    }
    if (mesh.foundObject<volScalarField>("ddt0(rhoE)"))
    {
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>("ddt0(rhoE)")
        ).timeIndex() = mesh.time().timeIndex();
    }
#endif

    if (!steadyState_)
    {
        // NOTE: Solving the system A \Delta W^\tau = B and therefore the
        // source term needs to be adjusted as below.
        // For steady state analysis this term reduces it to zero.
        //
        // NOTE: The residual is in strong form and therefore the division by V
        rhoR.primitiveFieldRef() -= (ddtRho().diag()*rho-ddtRho().source())/mesh.V();
        rhoUR.primitiveFieldRef() -= (ddtRhoU().diag()*rhoU-ddtRhoU().source())/mesh.V();
        rhoER.primitiveFieldRef() -= (ddtRhoE().diag()*rhoE-ddtRhoE().source())/mesh.V();
    }

    // NOTE: These should be equal if ddtScheme is consistent
    ddtCoeff = max(max(ddtRho().diag(),ddtRhoU().diag()),ddtRhoE().diag())/mesh.V();

    // =================================== residualsUpdate(end) ===================================

    // Pseudo Courant number relaxation
    setPseudoCoNum();

    // Update pseudo time step
    setPseudoDeltaT();

    // Set based on latest rPseudoDeltaT
    ddtCoeff = max(max(fvm::ddt(rho)->diag(),fvm::ddt(rhoU)->diag()),fvm::ddt(rhoE)->diag())/mesh.V();

    cellDebug();

    // Store previous iteration values. Don't use storePrevIter here as it
    // introduces too many problems for mesh refinement/redistribution
    scalarVarsPrevIter_.clear();
    vectorVarsPrevIter_.clear();
    forAll(scalarVars_, i)
    {
        scalarVarsPrevIter_.append(new volScalarField(scalarVars_[i].name() + "PrevIter", scalarVars_[i]));
    }
    forAll(vectorVars_, i)
    {
        vectorVarsPrevIter_.append(new volVectorField(vectorVars_[i].name() + "PrevIter", vectorVars_[i]));
    }

    // Solver
    {
        residualIO defaultSolverTol(2, 1, residualOrdering, 1e-12);

        // Storage of jacobians and preconditioners
        PtrList<jacobian> jacobians;
        PtrList<preconditioner<2,1>> preconditioners;

        const dictionary& dict = mesh.solutionDict().subDict("flowSolver");
        const word solverType(dict.lookup("solver"));

        // Create main Jacobian
        jacobians.append
        (
            new jacobian
            (
                dict.subOrEmptyDict(solverType),
                mesh,
                rho,
                rhoU,
                rhoE,
                thermo,
                U,
                ddtCoeff_(),
                inviscid_,
                turbulence_
            )
        );

        // Recursively create needed preconditioners and their jacobians (if
        // applicable)
        createPreconditioners
        (
            preconditioners,
            jacobians,
            dict.subOrEmptyDict(solverType)
        );

        // Create solver
        autoPtr<solver<2,1>> sol =
            solver<2,1>::New
            (
                dict,
                jacobians[0].matrix(),
                preconditioners.size() ? preconditioners(0) : NULL,
                defaultSolverTol
            );

        // Call solver
        initRes_.clear();
        sol->solve
        (
            scalarVars_,
            vectorVars_,
            scalarResiduals_,
            vectorResiduals_,
            initRes_
        );
        solnControl_->setResidual(initRes_());
    }

    if (localTimesteppingBounding_)
    {
        volScalarField rhoMin =
            localTimesteppingLowerBound_*scalarVarsPrevIter_[0];
        // volScalarField rhoMax = 2*scalarVarsPrevIter_[0];
        volScalarField eMin =
            localTimesteppingLowerBound_*
            (scalarVarsPrevIter_[1]/scalarVarsPrevIter_[0]
            - 0.5*magSqr(vectorVarsPrevIter_[0]/scalarVarsPrevIter_[0]));
        // volScalarField eMax =
        //     100.0*(scalarVarsPrevIter_[1]/scalarVarsPrevIter_[0]
        //     - 0.5*magSqr(vectorVarsPrevIter_[0]/scalarVarsPrevIter_[0]));
        volScalarField eTemp(rhoE/rho - 0.5*magSqr(rhoU/rho));

        volScalarField factor
        (
            IOobject
            (
                "factor",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("factor", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(mesh.owner(), facei)
        {
            label own = mesh.owner()[facei];
            label nei = mesh.neighbour()[facei];

            if
            (
                (rho[own] < rhoMin[own]) //|| (rho[own] > rhoMax[own])
            || (eTemp[own] < eMin[own]) //|| (eTemp[own] > eMax[own])
            || (eTemp[own] < SMALL)
            )
            {
                factor[own] = min(0.5, factor[own]);
                factor[nei] = min(0.75, factor[nei]);
            }
            if
            (
                (rho[nei] < rhoMin[nei]) //|| (rho[nei] > rhoMax[nei])
            || (eTemp[nei] < eMin[nei]) //|| (eTemp[nei] > eMax[nei])
            || (eTemp[nei] < SMALL)
            )
            {
                factor[nei] = min(0.5, factor[nei]);
                factor[own] = min(0.75, factor[own]);
            }
        }

        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                scalarField rhoNei =
                    rho.boundaryField()[patchi].patchNeighbourField();
                scalarField rhoMinNei =
                    rhoMin.boundaryField()[patchi].patchNeighbourField();
                // scalarField rhoMaxNei =
                //     rhoMax.boundaryField()[patchi].patchNeighbourField();
                scalarField eTempNei =
                    eTemp.boundaryField()[patchi].patchNeighbourField();
                scalarField eMinNei =
                    eMin.boundaryField()[patchi].patchNeighbourField();
                // scalarField eMaxNei =
                //     eMax.boundaryField()[patchi].patchNeighbourField();
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                forAll(fc, bfacei)
                {
                    if
                    (
                        (rho[fc[bfacei]] < rhoMin[fc[bfacei]])
                    //  || (rho[fc[bfacei]] > rhoMax[fc[bfacei]])
                    || (eTemp[fc[bfacei]] < eMin[fc[bfacei]])
                    //  || (eTemp[fc[bfacei]] > eMax[fc[bfacei]])
                    || (eTemp[fc[bfacei]] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] = min(0.5, factor[fc[bfacei]]);
                    }
                    if
                    (
                        (rhoNei[bfacei] < rhoMinNei[bfacei])
                    //  || (rhoNei[bfacei] > rhoMaxNei[bfacei])
                    || (eTempNei[bfacei] < eMinNei[bfacei])
                    //  || (eTempNei[bfacei] > eMaxNei[bfacei])
                    || (eTempNei[bfacei] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] = min(0.75, factor[fc[bfacei]]);
                    }
                }
            }
        }

        pseudoCoField_() *= factor;
    }

    phiUp_.clear();
    phiEp_.clear();
    Up_.clear();

    // =================================== UpdateFields ===================================

    // Optional failsafe to avoid -ve rho. Avoid dividing by very small quantity
    const dictionary& dict = mesh.solutionDict().subDict("pseudoTime");
    dimensionedScalar rhoMin
    (
        "rhoMin", 
        dimDensity, 
        dict.lookupOrDefault<scalar>("rhoMin", -GREAT)
    );
    bound(rho, rhoMin);

    // Compute internal fields

    U.ref() = rhoU.internalField()/rho.internalField();

    e.ref() = rhoE.internalField()/rho.internalField()
            - 0.5*magSqr(U.internalField()) + eZero_->internalField();

    // Temperature bound
    dimensionedScalar TMin
    (
        "TMin", 
        dimTemperature, 
        dict.lookupOrDefault<scalar>("TMin", ROOTSMALL)
    );
    dimensionedScalar TMax
    (
        "TMax", 
        dimTemperature, 
        dict.lookupOrDefault<scalar>("TMax", GREAT)
    );
    
    // Bound energy
    volScalarField TBound
    (
        IOobject("TBound", runTime.timeName(), mesh), 
        mesh,
        TMin
    );
    volScalarField eBound = thermo.he(thermo.p(), TBound);

    if (max(neg(e-eBound)).value() > 0.5)
    {
        e = max(e, eBound);
        bounded_ = true;

        Info<< "Bounding " << e.name()
            << " to TMin: " << TMin.value()
            << endl;
    }

    // Only do max bound if it was specified
    if (TMax.value() < GREAT)
    {
        TBound == TMax;
        eBound == thermo.he(thermo.p(), TBound);
        if (max(pos(e-eBound)).value() > 0.5)
        {
            e = min(e, eBound);
            Info<< "Bounding " << e.name()
                << " to TMax: " << TMax.value()
                << endl;
        }
    }

    // Calc T and psi from e
    thermo.correct();

    p.ref() = rho.internalField()/psi.internalField();

    // Recalc rhoU and rhoE in case rho or e were bounded
    rhoU.ref() = rho.internalField()*U.internalField();
    rhoE.ref() = 
        rho.internalField()*
        (
            e.internalField() - eZero_() + 0.5*magSqr(U.internalField())
        );

    // Correct boundary fields
    p.correctBoundaryConditions();
    U.correctBoundaryConditions();
    T.correctBoundaryConditions();

    T.boundaryFieldRef() == max(T.boundaryField(), TMin.value());
    if (TMax.value() < GREAT)
    {
        T.boundaryFieldRef() == max(T.boundaryField(), TMax.value());
    }

    e.boundaryFieldRef() == thermo.he(p, T)->boundaryField();

    thermo.correct();                       // NOTE: Correct psi boundaryField

    /*//rho.boundaryField() == psi.boundaryField()*p.boundaryField(); // Slight difference compared to when using rhoFix BC to update
    rho.correctBoundaryConditions();        // NOTE: rhoFix does a psi*p update (See eg charactPress)*/

    tmp<volScalarField> trho = thermo.rho();
    forAll(rho.boundaryField(), patchi)
    {
        // Substitute wall pressure with zero grad to avoid unphysical values
        if
        (
            isA<characteristicWallPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            rho.boundaryFieldRef()[patchi] =
                trho().boundaryField()[patchi].patchInternalField();
        }
        else
        {
            rho.boundaryFieldRef()[patchi] = trho().boundaryField()[patchi];
        }
    }
    rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();
    rhoE.boundaryFieldRef() =
        rho.boundaryField()*
        (
            e.boundaryField() 
        - eZero_->boundaryField() 
        + 0.5*magSqr(U.boundaryField())
        );


    // =================================== UpdateFields(end) ===================================

    // Correct turbulence
    turbulence_->correct();

    // TODO: check
    sigmaFluid = -rho*
        (
            p*symmTensor(1,0,0,1,0,1)
        + turbulence_->devRhoReff()
        );
}


void Foam::HiSAPhysicsSolver::startTimeStep()
{
    // Info<< "Start Time Step\n" << endl;
    // Clear out residuals from previous time step
    if (!steadyState_)
    {
        initRes_.clear();
        prevRes_.clear();
    }

    pseudotimeControl& solnControl = solnControl_();
    moveMeshOuterCorrectors_ =
        solnControl.dict().lookupOrDefault<label>
        (
            "moveMeshOuterCorrectors",
            0
        );

    // Reset pseudo Co number before every time step
    if (solnControl.dict().lookupOrDefault<bool>("resetPseudo", false) )
    {
        if (!localTimestepping_)
        {
            pseudoCoNum_->value() =
                solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
        }
        else
        {
            pseudoCoField_() == solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
        }
    }
    pseudoCoNumMin_ =
        solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMin", 0.1);
    pseudoCoNumMax_ =
        solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMax", 25);
    pseudoCoNumMaxIncr_ = solnControl.dict().lookupOrDefault<scalar>
    (
        "pseudoCoNumMaxIncreaseFactor",
        1.25
    );
    pseudoCoNumMinDecr_ = solnControl.dict().lookupOrDefault<scalar>
    (
        "pseudoCoNumMinDecreaseFactor",
        0.1
    );

    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    rebalance_ =
        solverControls.lookupOrDefault<Switch>("rebalance", false);
    maxLoadImbalance_ =
        solverControls.lookupOrDefault<scalar>("maxLoadImbalance", 0.1);
}

void Foam::HiSAPhysicsSolver::endTimeStep()
{
    // Info<< "End Time Step\n" << endl;
}

bool Foam::HiSAPhysicsSolver::computeDeltaT()
{
    //const Time& runTime = mesh_.time();

    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    bool computeDeltaT =
        solverControls.lookupOrDefault<bool>("computeDeltaT", false);

    return computeDeltaT;
}

scalar Foam::HiSAPhysicsSolver::returnDeltaT()
{
    const Time& runTime = mesh_.time();
    dynamicFvMesh& mesh = mesh_;
    psiThermo& thermo = pThermo_();
    volVectorField& U = U_();

    const dictionary& solverControls =
        mesh_.solutionDict().subDict("solverControls");

    scalar maxCo =
    solverControls.getOrDefault<scalar>("maxCo", 1);

    Info << "Mesh region: " << name() << endl;

    // ======================== Compute Courant Number ========================

    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    if (mesh.nInternalFaces())
    {
        tmp< surfaceScalarField > lambda;
        if (mesh.moving())
        {
            lambda = fvc::interpolate(sqrt(thermo.gamma()/thermo.psi())) + mag((fvc::interpolate(U)&mesh.Sf())-fvc::meshPhi(U))/mesh.magSf();
        }
        else
        {
            lambda = fvc::interpolate(sqrt(thermo.gamma()/thermo.psi())) + mag(fvc::interpolate(U)&mesh.Sf()/mesh.magSf());
        }

        surfaceScalarField amaxSfbyDelta
        (
            lambda*mesh.magSf()*mesh.deltaCoeffs()
        );

        CoNum = max(amaxSfbyDelta/mesh.magSf()).value()*runTime.deltaTValue();

        meanCoNum =
            (sum(amaxSfbyDelta)/sum(mesh.magSf())).value()
        *runTime.deltaTValue();
    }

    Info<< "Mean and max Courant Numbers = "
    << meanCoNum << " " << CoNum << endl;

    // ======================== Compute Courant Number(end) ========================

    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    return deltaTFact*runTime.deltaTValue();
}

// ************************************************************************* //
