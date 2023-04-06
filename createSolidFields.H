	Info<< "Reading field Tsolid\n" << endl;
    volScalarField Tsolid
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    Info<< "Reading solid diffusivity KT\n" << endl;

    volScalarField KTsolid
    (
        IOobject
        (
		"KT",
		runTime.timeName(),
		solidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
        ),
        solidMesh
    );
    
    IOdictionary thermophysicalProperties
(
	IOobject
	(
		"thermophysicalProperties",
		runTime.constant(),
		solidMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);



    // On creation of both diffusivities, couple the sides
    alphaEff.correctBoundaryConditions();
    KTsolid.correctBoundaryConditions();

  
