{
	// TODO: trzeba sprawdzic wallFunction alphatJayatille....
	// czy dobrze jest ten zmienny Pr
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

	alphaEff = kHe/(rho/rhok)/cp + alphat;
	//Info<< "alphaEff w TEqn.H = " << alphaEff << endl;
	//Info<< "alphat w TEqn.H = " << alphat << endl;
	//Info<< "rho w TEqn.H = " << rho << endl;
	//Info<< "rhok w TEqn.H = " << rhok << endl;
	//Info<< "cp w TEqn.H = " << cp << endl;
	//

	//volScalarField JT = (U & fvc::grad(p*(rho/rhok)))/(rho/rhok)/cp/T;
	//Info<< "JT = " << JT << endl;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
	//  + AGM*rhon*rhos*pow(sHe/AGM/rhon*magGradT,4./3)/(rho/rhok)/cp
	  + (
			mutualFrictionGM
		  ? AGM*rhon*rhos*pow(sHe/AGM/rhon*magGradT,4./3)/(rho/rhok)/cp
	      : AGM0
	      //: AGM*rhon*rhos*pow(sHe/AGM/rhon*magGradT,4./3)/(rho/rhok)/cp*0
		)
      + fvOptions(T)
	  // Additional term 
	  //+ fvm::SuSp(- (U & fvc::grad(p*(rho/rhok)))/(rho/rhok)/cp/T, T)
	  - (
			JT
		  ? fvm::SuSp((U & fvc::grad(p*(rho/rhok)))/(rho/rhok)/cp/T, T)
		  : fvm::Sp(JT0, T)
		)
	  //- ((U & fvc::grad(p*(rho/rhok)))/(rho/rhok)/cp)
	  // The rest of the JT effect
	  + (
			JT
		  ? fvm::SuSp( (U & fvc::grad(p*(rho/rhok)))*beta/(rho/rhok)/cp, T)
		  : fvm::Sp(JT0, T)
		)
	  //- fvm::SuSp(- (U & fvc::grad(p*(rho/rhok)))*beta/(rho/rhok)/cp, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);
    //T.correctBoundaryConditions();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;


	    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
		

	fvScalarMatrix TSolidEqn
            (
               rho_solid*cp_solid*fvm::ddt(Tsolid) - fvm::laplacian(KTsolid, Tsolid)
            );
    TSolidEqn.solve();
	

    Info<< "min/max(Tsolid) = " << min(Tsolid).value() << ", "
        << max(Tsolid).value() <<endl;

	// beta musi byc updated
    //rhok = 1.0 - beta*(T - TRef);
}
