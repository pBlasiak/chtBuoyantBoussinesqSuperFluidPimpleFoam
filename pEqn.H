{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());
    surfaceScalarField phiM(-rAUf*(fvc::interpolate(M) & mesh.Sf()));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi))
      + phig
      + phiM
    );
	//Info<< "max(phig) = " << max(phig) << endl;
	//Info<< "max(phiM) = " << max(phiM) << endl;

    MRF.makeRelative(phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phiM + phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

	//Un liczone jako nowe U odjac zmienne z poprzedniego dt???
	Un = U - pow
		     (
				pow(rhos,3)*sHe/AGM/pow(rho,3)/rhon/max
				(
					//pow(mag(fvc::grad(T)),2), 
					magGradT2,
					dimensionedScalar("small", dimTemperature*dimTemperature/dimLength/dimLength, SMALL)
				//), 1./3)*fvc::grad(T);
				), 1./3)*gradT;
	Un.correctBoundaryConditions();
	Us = U + pow
		     (
				pow(rhon,3)*sHe/AGM/pow(rho,3)/rhon/max
				(
					//pow(mag(fvc::grad(T)),2), 
					magGradT2,
					dimensionedScalar("small", dimTemperature*dimTemperature/dimLength/dimLength, SMALL)
				//), 1./3)*fvc::grad(T);
				), 1./3)*gradT;
	Us.correctBoundaryConditions();
	UnmUs = Un - Us;
	
    Info<< "min/max(U) = " << min(mag(U)).value() << ", "
        << max(mag(U)).value() <<endl;

    p = p_rgh + rhok*gh;
	P = p*(rho/rhok);

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
}
