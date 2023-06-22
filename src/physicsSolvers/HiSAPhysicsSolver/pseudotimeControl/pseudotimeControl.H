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
    pseudotimeControl

Description
    Adaptation of the pseudo-time control library for the "multiFieldSolver" 
    application. Adapted and build with respect to the library HiSA version 
    "hisa 1.4.8" for OpenFOAM-v2006. The source of the HISA solver can be 
    found here:

    https://hisa.gitlab.io
    https://gitlab.com/hisa/hisa

SourceFiles
    pseudotimeControl.C
    
\*---------------------------------------------------------------------------*/

#ifndef pseudotimeControl_H
#define pseudotimeControl_H


#include "solutionControl.H"
#include "residualIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pseudotimeControl Declaration
\*---------------------------------------------------------------------------*/

class pseudotimeControl
:
public solutionControl
{
protected:

    // Protected data

        // Solution controls

            //- Steady state flag
            bool steadyState_;

            //- Maximum number of outer (pseudotime) correctors
            int nCorrOuter_; // Will be adjusted by multiFieldSolver 

            //- Minimum number of outer (pseudotime) correctors
            int nCorrOuterMin_; // Will be adjusted by multiFieldSolver

            //- FieldDictionary from the multiFieldSolver
            IOdictionary fieldsDictionary_;

            //- Flag to indicate whether to only solve turbulence on final iter
            bool turbOnFinalIterOnly_;

            //- Converged flag
            bool converged_;

            //- First iteration or first iteration after restart (no residual computed yet)
            bool firstIteration_;

        // Residuals

            residualIO residualTols_;
            residualIO residualTolsRel_;
            residualIO residuals_;
            residualIO initResiduals_;

            // For persistent storage of initial residuals
            IOdictionary state_;

    // Protected Member Functions

        //- Read controls from fvSolution dictionary
#if OPENFOAM >= 1906
        virtual bool read();
#else
        virtual void read();
#endif

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Disallow default bitwise copy construct
        pseudotimeControl(const pseudotimeControl&);

        //- Disallow default bitwise assignment
        void operator=(const pseudotimeControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("pseudotimeControl");


    // Constructors

        //- Construct from mesh
        pseudotimeControl
        (
            fvMesh& mesh,
            const bool steadyState,
            const label nScalars,
            const label nVectors,
            const residualIO& defaultTol,
            const residualIO& defaultTolRel
        );


    //- Destructor
    virtual ~pseudotimeControl();


    // Member Functions

        // Access

            //- Maximum number of pseudotime correctors
            inline label nCorrOuter() const; // Will be adjusted by multiFieldSolver


        // Solution control

            //- pseudotime loop
            virtual bool loop();

            //- Helper function to identify when to store the intial residuals
            inline bool storeInitialResiduals() const;

            //- Helper function to identify first outer (pseudotime) iteration
            inline bool firstIter() const;

            //- Helper function to identify final outer (pseudotime) iteration
            inline bool finalIter() const;

            //- Helper function to identify whether to solve for turbulence
            inline bool turbCorr() const;

            //- Helper function to store residual (initial solver residual)
            inline void setResidual(const residualIO& res);

            inline void setCorr(const int& corr)
            {
                corr_ = corr;
            }

            //- initialise variable for a new time-step
            inline void init();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pseudotimeControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
