/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "shearStressMarangoniFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H" //Important!!! Otherwise: 'patch' has not been declared

#include "volFields.H"
#include "transformFvPatchFields.H"
#include "symmTransformField.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//template<class Type>
//Foam::shearStressMarangoniFvPatchVectorField<Type>::shearStressMarangoniFvPatchVectorField
  shearStressMarangoniFvPatchVectorField::shearStressMarangoniFvPatchVectorField
(
    const fvPatch& p,
      const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchVectorField(p, iF),
    gradTName_(),
    MAcoef_()
{}


//template<class Type>
//Foam::shearStressMarangoniFvPatchVectorField<Type>::shearStressMarangoniFvPatchVectorField
shearStressMarangoniFvPatchVectorField::shearStressMarangoniFvPatchVectorField
(
      const shearStressMarangoniFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
      transformFvPatchVectorField(ptf, p, iF, mapper),
    gradTName_(ptf.gradTName_),
    MAcoef_(ptf.MAcoef_)
{}

shearStressMarangoniFvPatchVectorField::shearStressMarangoniFvPatchVectorField
(
    const fvPatch& p,
    //const DimensionedField<Type, volMesh>& iF,
      const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
      transformFvPatchVectorField(p, iF),
    gradTName_( dict.lookup("gradTName") ),
    MAcoef_(readScalar( dict.lookup("MAcoef") ) )
{

    evaluate();

}


shearStressMarangoniFvPatchVectorField::shearStressMarangoniFvPatchVectorField
(
      const shearStressMarangoniFvPatchVectorField& ptf
)
:
    transformFvPatchVectorField(ptf),
    gradTName_(ptf.gradTName_),
    MAcoef_(ptf.MAcoef_)
{}


  shearStressMarangoniFvPatchVectorField::shearStressMarangoniFvPatchVectorField
(
      const shearStressMarangoniFvPatchVectorField& ptf,
      const DimensionedField<vector, volMesh>& iF
)
:
      transformFvPatchVectorField(ptf, iF),
    gradTName_(ptf.gradTName_),
    MAcoef_(ptf.MAcoef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


  void shearStressMarangoniFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
      transformFvPatchVectorField::autoMap(m);
}


  void shearStressMarangoniFvPatchVectorField::rmap
(
      const fvPatchVectorField& ptf,
    const labelList& addr
)
{
      transformFvPatchVectorField::rmap(ptf, addr);

      const shearStressMarangoniFvPatchVectorField& dmptf =
          refCast<const shearStressMarangoniFvPatchVectorField>(ptf);

    gradTName_ = dmptf.gradTName_;
    MAcoef_ = dmptf.MAcoef_;
}


tmp<vectorField>
shearStressMarangoniFvPatchVectorField::snGrad() const
{
    const vectorField nHat = this->patch().nf();
    const vectorField pif = this->patchInternalField();
    vectorField result;

    if( db().foundObject<vectorField>(gradTName_))
    {
      const fvPatchField<vector>& gradT = patch().lookupPatchField<volVectorField, vector>(gradTName_);
      // Info << "gradT found: " << gradTName_;
      // Info << " gradT max " << max(gradT);
      // Info << " gradT min " << min(gradT) << endl;
      vectorField gradT_internal = gradT.patchInternalField();
      vectorField gradTplane = transform(I - sqr(nHat), gradT_internal);
      vectorField pifplane = transform(I - sqr(nHat), pif);
      result = pifplane + ( MAcoef_*gradTplane )/this->patch().deltaCoeffs();
    }
    else 
    { 
      result = transform(I - sqr(nHat), pif);
    }

    return (result-pif)*this->patch().deltaCoeffs();
}

void shearStressMarangoniFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    // Info << "entering  marangoniFvPatchVectorField::evaluate()" << endl;
    if (!this->updated())
    {
      // Info << "marangoniFvPatchVectorField::evaluate(): calling updatecoeffs" << endl;
      this->updateCoeffs();
    }

    const vectorField nHat = this->patch().nf();
    const vectorField pif = this->patchInternalField();

    vectorField result;
    scalarField gradTn;

    if( db().foundObject<vectorField>(gradTName_))
    {
      const fvPatchField<vector>& gradT = patch().lookupPatchField<volVectorField, vector>(gradTName_);
      // Info << "gradT found: " << gradTName_;
      // Info << " gradT max " << max(gradT);
      // Info << " gradT min " << min(gradT) << endl;
      vectorField gradT_internal = gradT.patchInternalField();
      vectorField gradTplane = transform(I - sqr(nHat), gradT_internal);
      vectorField pifplane = transform(I - sqr(nHat), pif);
      // Info << "gradTplane: " << gradTName_;
      // Info << " gradTplane max " << max(gradTplane);
      // Info << " gradTplane min " << min(gradTplane) << endl;


      result = pifplane + ( MAcoef_*gradTplane)/(this->patch().deltaCoeffs());
    }
    else
    {
      result = transform(I - sqr(nHat), pif);
    }
    vectorField::operator=
    (
      result
    );

    transformFvPatchVectorField::evaluate();
}

tmp<vectorField>
shearStressMarangoniFvPatchVectorField::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return
      transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
}

  void shearStressMarangoniFvPatchVectorField::write(Ostream& os) const
{
    transformFvPatchVectorField::write(os);
    os.writeKeyword("gradTName") << gradTName_ << token::END_STATEMENT << nl;
    os.writeKeyword("MAcoef") << MAcoef_ << token::END_STATEMENT << nl;
    writeEntry("value", os);

}

// ************************************************************************* //

makePatchTypeField(fvPatchVectorField, shearStressMarangoniFvPatchVectorField);

} // End namespace Foam

// ************************************************************************* //
