/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::waveTransmissiveCompressibleFluidFvPatchField

Description
    Foam::waveTransmissiveCompressibleFluidFvPatchField

SourceFiles
    waveTransmissiveCompressibleFluidFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef waveTransmissiveCompressibleFluidFvPatchField_H
#define waveTransmissiveCompressibleFluidFvPatchField_H

#include "advectiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class waveTransmissiveCompressibleFluidFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class waveTransmissiveCompressibleFluidFvPatchField
:
    public advectiveFvPatchField<Type>
{

    // Private data

        //- Fluid bulk modulus
        const scalar K_;

        //- Fluid density
		const scalar rho_;


public:

    //- Runtime type information
    TypeName("waveTransmissiveCompressibleFluid");


    // Constructors

        //- Construct from patch and internal field
        waveTransmissiveCompressibleFluidFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveTransmissiveCompressibleFluidFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given waveTransmissiveCompressibleFluidFvPatchField
        //  onto a new patch
        waveTransmissiveCompressibleFluidFvPatchField
        (
            const waveTransmissiveCompressibleFluidFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waveTransmissiveCompressibleFluidFvPatchField
        (
            const waveTransmissiveCompressibleFluidFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new waveTransmissiveCompressibleFluidFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveTransmissiveCompressibleFluidFvPatchField
        (
            const waveTransmissiveCompressibleFluidFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new waveTransmissiveCompressibleFluidFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access


        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "waveTransmissiveCompressibleFluidFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
