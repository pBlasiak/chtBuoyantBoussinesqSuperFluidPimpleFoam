// Solve the momentum equation

MRF.correctBoundaryVelocity(U);

//czy ta zmienna jest widoczna w TEqn czy znika zaraz po 
//opuszczeniu UEqn?
//volVectorField gradT(fvc::grad(T));

// Update properties for M
rho = laminarTransport.rhoHe();
rhon = laminarTransport.rhon();
rhos = laminarTransport.rhos();
AGM = laminarTransport.AGMHe();
sHe = laminarTransport.sHe();
beta = laminarTransport.betaHe();
cp = laminarTransport.cpHe();
rhok = 1.0 - beta*(T - TRef);
if(gradP)
{
	Info<< "Licze z gradP!!!!" << endl;
	gradT = fvc::grad(T)-fvc::grad(p*(rho/rhok))/rho/sHe;	
}
else
{
	gradT = fvc::grad(T);	
}
magGradT = max(mag(gradT), dimensionedScalar("small", dimTemperature/dimLength, SMALL)); 
//Info<< "magGradT = " << magGradT << endl;
magGradT2 = pow(magGradT,2);


switch (kHeModel)
{
	case HeliumThermalConductivityModel::Hepak:
		 {
			kHe = pow(laminarTransport.onebyf()/magGradT2, 1./3);  // HEPAK
		 	break;
		 }

	case HeliumThermalConductivityModel::Kitamura:
		 {
			kHe = pow(1./(AGM*rhon/pow(rhos*sHe*T,3)/sHe)/magGradT2, 1./3); // KITAMURA
		 	break;
		 }

	case HeliumThermalConductivityModel::Suekane:
		 {
			const dimensionedScalar sLambda("sLambda", dimensionSet(0,2,-2,-1,0,0,0), 1560);
		 	const dimensionedScalar Ac("Ac", dimensionSet(-1,1,1,0,0,0,0), 1297);
		 	const dimensionedScalar TLambda("TLambda", dimTemperature, 2.168);
		 	const dimensionedScalar c("c", dimTemperature, 0.01);
		 	kHe = pow(1./(Ac/rho/rho/pow(sLambda,4)/pow(TLambda,3)*
		 				pow(pow((T+c)/(TLambda+c),5.7)*(1-pow((T+c)/(TLambda+c),5.7)),-3))
		 			/magGradT2, 1./3); // SUEKANE
		 	//Info<< "f^-1 = " << 1./(Ac/rho/rho/pow(sLambda,4)/pow(TLambda,3)*
		 	//			pow(pow((T+c)/(TLambda+c),5.7)*(1-pow((T+c)/(TLambda+c),5.7)),-3));
		 	break;
		 }

	default:
		 FatalErrorInFunction
    	     << "Unknown helium thermal conductivity model. Valid models are: "
    	     << HeliumThermalConductivityModelNames_ << nl << exit(FatalError);
		 break;
}




//volScalarField onebyf
//(
//    IOobject
//    (
//        "onebyf",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//	mesh,
//    dimensionedScalar("onebyf", dimensionSet(3,1,-9,-1,0,0,0), 1.169728e13)
//);
//	kHe = pow(onebyf/magGradT2, 1./3); // VAN SCIVER
	
	// To jest potrzebne do alphatJayatille....
	Pr = cp*laminarTransport.etaHe()/kHe;

	//Info<< "SFTC = " << laminarTransport.onebyf() << endl;
	//Info<< "1/f from Kitamura = " << 1./(AGM*rhon/pow(rhos*sHe*T,3)/sHe) << endl;
	volVectorField M
	(
		fvc::div
		(
			rhon*rhos/rho*pow(sHe/AGM/rhon/magGradT2, 2./3)*gradT*gradT
		)/(rho/rhok)
	  + laminarTransport.etaHe()/(rho/rhok)*pow
	    (
			pow(rhos,3)*sHe/AGM/pow(rho,3)/rhon/magGradT2, 1./3
		)*
		(
			fvc::laplacian(gradT) + 1./3*fvc::grad(fvc::div(gradT))
		)
	);
	//Info<< "max(M) = " << max(M) << endl;

	// TODO: zastanowic sie czy nie powinno byc
	// rho pod pochodna jak jest w buoyantPimpleFoam
	// ale w sumie stosujemy aproksymacje Boussinesqa wiec jest ok:
	// In flows accompanied by heat transfer, the fluid properties are normally functions of temperature. The variations may be small and yet be the cause of the fluid motion. If the density variation is not large, one may treat the density as constant in the unsteady and convection terms, and treat it as variable only in the gravitational term. This is called the Boussinesq approximation.
	// We have to know when these approximations are valid. Ferziger and Peric (2001, p.15) states that:

    //This approximation introduces errors of the order of 1% if the temperature differences are below e.g. 2℃ for water and 15℃ for air. The error may be more substantial when temperature differences are larger; the solution may even be qualitatively wrong.
	// https://caefn.com/openfoam/solvers-buoyantboussinesqpimplefoam
    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
	  //- fvm::laplacian(laminarTransport.etaHe()/(rho/rhok), U)
      //- fvc::div(laminarTransport.etaHe()/(rho/rhok)*dev(fvc::grad(U)().T()))
     ==
        fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
				  - (fvc::interpolate(M) & mesh.Sf())
            )
        );

        fvOptions.correct(U);
    }
