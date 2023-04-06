Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
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
);

#include "createPhi.H"

#include "readTransportProperties.H"


volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	laminarTransport.cpHe()
);

Info<< "Creating parameters for solid\n" << endl;

IOdictionary transportProperties(
    IOobject(
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

dimensionedScalar rho_solid(
    transportProperties.lookup("rho_solid"));

dimensionedScalar cp_solid(
    transportProperties.lookup("Cp_solid"));



autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Thermal expansion coefficient [1/K]
volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh
    ),
    laminarTransport.betaHe()
);

// Kinematic density for buoyancy force
// rhok = rho/rho0
// rho0 is a reference density at TRef
// rho = rho0*(1-beta*(T-TRef))
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta*(T - TRef)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	laminarTransport.rhoHe()
);

volScalarField rhon
(
    IOobject
    (
        "rhon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	rho
);
	rhon = laminarTransport.rhon();

volScalarField rhos
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	rho
);
	rhos = rho - rhon;

volScalarField sHe
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	laminarTransport.sHe()
);

volScalarField AGM
(
    IOobject
    (
        "AGM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	laminarTransport.AGMHe()
);

tmp<volScalarField> AGM0
(
	new volScalarField 
	(
		IOobject
    	(
    	    "AGM0",
    	    runTime.timeName(),
    	    mesh,
    	    IOobject::NO_READ,
    	    IOobject::NO_WRITE
    	),
		mesh,
		dimensionedScalar("AGM0", dimensionSet(0,0,-1,1,0,0,0), 0)
	)
);

const dimensionedScalar JT0("JT0", dimensionSet(0,0,-1,0,0,0,0), 0);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (p_rgh + rhok*gh)*(rho/rhok)
);

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"
#include "createFvOptions.H"

// normal velocity component
volVectorField Un
(
    IOobject
    (
        "Un",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// super velocity component
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField UnmUs
(
    IOobject
    (
        "UnmUs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Un - Us
);

volScalarField alphaEff
(
    IOobject
    (
        "alphaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("small", dimensionSet(0,2,-1,0,0,0,0), SMALL)
);

volScalarField kHe
(
    IOobject
    (
        "kHe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
    dimensionedScalar("small", dimensionSet(1,1,-3,-1,0,0,0), SMALL)
);

surfaceScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
fvc::interpolate(rhos*sHe*T*(Un.component(1) - Us.component(1)))
);

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("small", dimensionSet(0,-1,0,1,0,0,0), vector(0,0,0))
);
gradT = fvc::grad(T);

volScalarField magGradT
(
    IOobject
    (
        "magGradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
    dimensionedScalar("small", dimensionSet(0,-1,0,1,0,0,0), SMALL)
);
magGradT = mag(gradT);

volScalarField magGradT2
(
    IOobject
    (
        "magGradT2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
    dimensionedScalar("small", dimensionSet(0,-2,0,2,0,0,0), SMALL)
);

