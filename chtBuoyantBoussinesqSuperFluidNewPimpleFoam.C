/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    chtBuoyantBoussinesqPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseHeliumTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent flow"
        " of incompressible fluids.\n"
        "Uses the Boussinesq approximation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createFluidMesh.H"
    #include "createSolidMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createSolidFields.H"
    #include "createTimeControls.H"
    #include "getCellDims.H"
    #include "CourantNo.H"
    #include "SuekaneNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
	    #include "createTimeControls.H"
		#include "SuekaneNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
			laminarTransport.correct();

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

		// updates q
		// teraz q liczone na bazie stare rhos i sHe i nowe T, Un i Us
		q = fvc::interpolate(rhos*sHe*T*(Un.component(1) - Us.component(1)));
		const surfaceScalarField::Boundary& patchHeatFlux = q.boundaryField();
        Info<< "\nWall heat fluxes " << endl;
        forAll(patchHeatFlux, patchi)
        {
           if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
            {
                Info<< mesh.boundary()[patchi].name()
                    << ": Total "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )
                    << " [W] over "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [m2] ("
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )/
                       gSum 
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [W/m2])"
                    << endl;
            }
	}


		// Tutaj pewna niescislosc bo nowo liczony gradT a stare alphaEff, cp i rho
        surfaceScalarField heatFluxFromAlphaEff =
			fvc::interpolate(alphaEff*cp*(rho/rhok))*(fvc::snGrad(T));
			//fvc::interpolate(alphaEff*cp*rho)*(gradT + gradp);

        const surfaceScalarField::Boundary& patchHeatFlux2 =
                 heatFluxFromAlphaEff.boundaryField();

        Info<< "\nWall heat fluxes from alphaEff" << endl;
        forAll(patchHeatFlux2, patchi)
        {
           if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
            {
                Info<< mesh.boundary()[patchi].name()
                    << ": Total "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux2[patchi]
                       )
                    << " [W] over "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [m2] ("
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux2[patchi]
                       )/
                       gSum 
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [W/m2])"
                    << endl;
            }
        }
        runTime.write();

        runTime.printExecutionTime(Info);
    }
	
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
