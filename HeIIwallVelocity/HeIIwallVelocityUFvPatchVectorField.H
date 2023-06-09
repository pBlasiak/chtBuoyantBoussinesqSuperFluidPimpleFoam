/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::HeIIwallVelocityUFvPatchVectorField

Group
    grpWallBoundaryConditions

Description
    This boundary condition provides a velocity condition for cases with
    simplified super fluid helium models. 

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            HeIIwallVelocityU;
        value           uniform (0 0 0);    // Initial value
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchVectorField

SourceFiles
    HeIIwallVelocityUFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef HeIIwallVelocityUFvPatchVectorField_H
#define HeIIwallVelocityUFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class HeIIwallVelocityUFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class HeIIwallVelocityUFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("HeIIwallVelocityU");


    // Constructors

        //- Construct from patch and internal field
        HeIIwallVelocityUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        HeIIwallVelocityUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given HeIIwallVelocityUFvPatchVectorField
        //  onto a new patch
        HeIIwallVelocityUFvPatchVectorField
        (
            const HeIIwallVelocityUFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        HeIIwallVelocityUFvPatchVectorField
        (
            const HeIIwallVelocityUFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new HeIIwallVelocityUFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        HeIIwallVelocityUFvPatchVectorField
        (
            const HeIIwallVelocityUFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new HeIIwallVelocityUFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
	
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
