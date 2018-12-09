Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
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


	Info<< "Solid volume η\n" << endl;
    volScalarField eta
    (
        IOobject
        (
            "eta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
		,
		dimensionedScalar("eta", dimless, scalar(0.0))
    );	
	
	Info<< "Body force bodyforce\n" << endl;
    volVectorField bodyforce
    (
        IOobject
        (
            "bodyforce",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
		,
		dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0),vector::zero)
    );	
	
	
    volScalarField x
    (
        IOobject
        (
            "x",
            runTime.constant(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
		dimensionedScalar("zero",dimensionSet(0,1,0,0,0,0,0),0.0)
    );
	
    volScalarField y
    (
        IOobject
        (
            "y",
            runTime.constant(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
		dimensionedScalar("zero",dimensionSet(0,1,0,0,0,0,0),0.0)
    );
	
	volScalarField z
    (
        IOobject
        (
            "z",
            runTime.constant(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
		dimensionedScalar("zero",dimensionSet(0,1,0,0,0,0,0),0.0)
    );
	
	x = mesh.C().component(vector::X);
	y = mesh.C().component(vector::Y);
	z = mesh.C().component(vector::Z);



#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
