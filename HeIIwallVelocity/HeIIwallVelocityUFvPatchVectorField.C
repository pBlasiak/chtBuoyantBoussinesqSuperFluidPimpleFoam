/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "HeIIwallVelocityUFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HeIIwallVelocityUFvPatchVectorField::
HeIIwallVelocityUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::HeIIwallVelocityUFvPatchVectorField::
HeIIwallVelocityUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::HeIIwallVelocityUFvPatchVectorField::
HeIIwallVelocityUFvPatchVectorField
(
    const HeIIwallVelocityUFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::HeIIwallVelocityUFvPatchVectorField::
HeIIwallVelocityUFvPatchVectorField
(
    const HeIIwallVelocityUFvPatchVectorField& sfwpvf
)
:
    fixedValueFvPatchVectorField(sfwpvf)
{}


Foam::HeIIwallVelocityUFvPatchVectorField::
HeIIwallVelocityUFvPatchVectorField
(
    const HeIIwallVelocityUFvPatchVectorField& sfwpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(sfwpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HeIIwallVelocityUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Accessing face normal vectors
    const vectorField nf(patch().nf());

	// Accessing G vector at the patch
	const fvPatchField<vector>& Gwall = patch().lookupPatchField<volVectorField, vector>("gradT");

    // Calculated the component of the G vector parallel to the wall
    vectorField GwallTang = Gwall - (nf & Gwall)*nf;
	//Info<< patch().name() << endl;
	//Info<< GwallTang << endl;

	// Accessing other variables at the patch
	const fvPatchField<scalar>& rhoWall = patch().lookupPatchField<volScalarField, scalar>("rho");
	const fvPatchField<scalar>& rhosWall = patch().lookupPatchField<volScalarField, scalar>("rhos");
	const fvPatchField<scalar>& rhonWall = patch().lookupPatchField<volScalarField, scalar>("rhon");
	const fvPatchField<scalar>& sWall = patch().lookupPatchField<volScalarField, scalar>("sHe");
	const fvPatchField<scalar>& AGMWall = patch().lookupPatchField<volScalarField, scalar>("AGM");
	const fvPatchField<scalar>& magGradT2Wall = patch().lookupPatchField<volScalarField, scalar>("magGradT2");

    vectorField::operator==(pow
			(
				pow(rhosWall,3)*sWall/AGMWall/pow(rhoWall,3)/rhonWall/magGradT2Wall, 
				1./3
			)*GwallTang);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::HeIIwallVelocityUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        HeIIwallVelocityUFvPatchVectorField
    );
}

// ************************************************************************* //
