/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "JinZeroFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
/*
Foam::scalar Foam::JinZeroFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//
Foam::JinZeroFvPatchScalarField::
JinZeroFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    //beta_(p.size(), Zero),
    diffCoeffName_("DiffCoeffName")
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1.0;
}


Foam::JinZeroFvPatchScalarField::
JinZeroFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict              //
)
:
    mixedFvPatchScalarField(p, iF),
    //beta_("beta", dict, p.size()),
    diffCoeffName_(dict.lookupOrDefault<word>("DiffCoeffName", "wordDefault"))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    refValue() = *this;
    */
}


Foam::JinZeroFvPatchScalarField::
JinZeroFvPatchScalarField
(
    const JinZeroFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    //beta_(mapper(ptf.beta_)),
    diffCoeffName_(ptf.diffCoeffName_)
{}


Foam::JinZeroFvPatchScalarField::
JinZeroFvPatchScalarField
(
    const JinZeroFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    //beta_(ptf.beta_),
    diffCoeffName_(ptf.diffCoeffName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//void Foam::JinZeroFvPatchScalarField::autoMap
//(
//    const fvPatchFieldMapper& m
//)
//{
//    mixedFvPatchScalarField::autoMap(m);
//    //m(beta_, beta_);
//}
//
//
//void Foam::JinZeroFvPatchScalarField::rmap
//(
//    const fvPatchScalarField& ptf,
//    const labelList& addr
//)
//{
//    mixedFvPatchScalarField::rmap(ptf, addr);
//
//    const JinZeroFvPatchScalarField& tiptf =
//        refCast<const JinZeroFvPatchScalarField>(ptf);
//
//    //beta_.rmap(tiptf.beta_, addr);
//}


void Foam::JinZeroFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& diffCoeffp =
        patch().lookupPatchField<volScalarField, scalar>(diffCoeffName_);

    refGrad() = 0;
    refValue() = 0;
    valueFraction() = 0.5/(0.5 + diffCoeffp*patch().deltaCoeffs());

    mixedFvPatchScalarField::updateCoeffs();
}



void Foam::JinZeroFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);


    //writeEntry(os, "beta", beta_);
    writeEntry(os, "diffCoeffName", diffCoeffName_);
    writeEntry(os, "refValue", refValue());
    writeEntry(os, "refGradient", refGrad());
    writeEntry(os, "valueFraction", valueFraction());
//    writeEntry(os, "value", value);

}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        JinZeroFvPatchScalarField
    );
}

// ************************************************************************* //
