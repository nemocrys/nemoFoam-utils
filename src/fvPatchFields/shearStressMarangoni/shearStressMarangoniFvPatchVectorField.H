/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2011-2017 OpenCFD Ltd.
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
    Foam::shearStressMarangoniFvPatchVectorField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a partial slip condition.  The amount of
    slip is controlled by a user-supplied field.

Usage
    \table
        Property      | Description             | Required    | Default value
        gradTName     | Name of Temperature gradient field (has to be calculated in solver) | yes     | 
        MAcoef | Factor for Marangoni force: Marangoni coefficient / dynamic viscosity | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            shearStressMarangoni;
        gradTName       gradT;
        MAcoef          -1.321; 
    }
    \endverbatim

See also
    Foam::transformFvPatchField

SourceFiles
    shearStressMarangoniFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef shearStressMarangoniFvPatchVectorField_H
#define shearStressMarangoniFvPatchVectorField_H

//#include "transformFvPatchField.H" //NOT THIS ONE!!!
#include "fvPatchFields.H"
#include "transformFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class shearStressMarangoniFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

//template<class Type>
class shearStressMarangoniFvPatchVectorField
:
    //public transformFvPatchField<Type>
      public transformFvPatchVectorField
{
    // Private data

        word gradTName_;
        scalar MAcoef_;

public:

    //- Runtime type information
    TypeName("shearStressMarangoni");


    // Constructors

        //- Construct from patch and internal field
        shearStressMarangoniFvPatchVectorField
        (
            const fvPatch&,
              const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        shearStressMarangoniFvPatchVectorField
        (
            const fvPatch&,
              const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given shearStressMarangoniFvPatchVectorField onto a new patch
        shearStressMarangoniFvPatchVectorField
        (
              const shearStressMarangoniFvPatchVectorField&,
            const fvPatch&,
              const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        shearStressMarangoniFvPatchVectorField
        (
              const shearStressMarangoniFvPatchVectorField&
        );

        //- Construct and return a clone
          virtual tmp<fvPatchVectorField> clone() const
        {
              return tmp<fvPatchVectorField>
            (
                  new shearStressMarangoniFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        shearStressMarangoniFvPatchVectorField
        (
              const shearStressMarangoniFvPatchVectorField&,
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
                  new shearStressMarangoniFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<vectorField> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return face-gradient transform diagonal
            virtual tmp<vectorField> snGradTransformDiag() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "shearStressMarangoniFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
