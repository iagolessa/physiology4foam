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

\*---------------------------------------------------------------------------*/

#include "waveTransmissiveCompressibleFluidFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
waveTransmissiveCompressibleFluidFvPatchField<Type>::waveTransmissiveCompressibleFluidFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    K_(0.0),
    rho_(0.0)
{}


template<class Type>
waveTransmissiveCompressibleFluidFvPatchField<Type>::waveTransmissiveCompressibleFluidFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    K_(readScalar(dict.lookup("K"))),
    rho_(readScalar(dict.lookup("rho")))
{}


template<class Type>
waveTransmissiveCompressibleFluidFvPatchField<Type>::waveTransmissiveCompressibleFluidFvPatchField
(
    const waveTransmissiveCompressibleFluidFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    K_(ptf.K_),
    rho_(ptf.rho_)
{}


template<class Type>
waveTransmissiveCompressibleFluidFvPatchField<Type>::waveTransmissiveCompressibleFluidFvPatchField
(
    const waveTransmissiveCompressibleFluidFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    K_(ptpsf.K_),
    rho_(ptpsf.rho_)
{}


template<class Type>
waveTransmissiveCompressibleFluidFvPatchField<Type>::waveTransmissiveCompressibleFluidFvPatchField
(
    const waveTransmissiveCompressibleFluidFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    K_(ptpsf.K_),
    rho_(ptpsf.rho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<scalarField> waveTransmissiveCompressibleFluidFvPatchField<Type>::advectionSpeed() const
{

    const surfaceScalarField& phi =
        this->db().objectRegistry::template lookupObject<surfaceScalarField>
        (this->phiName_);

    fvsPatchField<scalar> phip = this->lookupPatchField
    (
        this->phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop = this->lookupPatchField
        (
            this->rhoName_,
            reinterpret_cast<const volScalarField*>(0),
            reinterpret_cast<const scalar*>(0)
        );

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return phip/this->patch().magSf() + sqrt(K_/rho_);
}

template<class Type>
void waveTransmissiveCompressibleFluidFvPatchField<Type>::write(Ostream& os) const
{
    advectiveFvPatchField<Type>::write(os);

    os.writeKeyword("K") << K_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
