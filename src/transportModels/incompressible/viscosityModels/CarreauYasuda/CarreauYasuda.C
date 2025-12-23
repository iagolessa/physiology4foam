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

#include "CarreauYasuda.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(CarreauYasuda, 0);
    addToRunTimeSelectionTable
    (
        viscosityModel,
        CarreauYasuda,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::CarreauYasuda::calcNu() const
{
    return
        nuInf_
      + (nu0_ - nuInf_)
       *pow(scalar(1) + pow(k_*strainRate(), a_), (n_ - 1.0)/a_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::CarreauYasuda::CarreauYasuda
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    CarreauYasudaCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nu0_(CarreauYasudaCoeffs_.lookup("nu0")),
    nuInf_(CarreauYasudaCoeffs_.lookup("nuInf")),
    k_(CarreauYasudaCoeffs_.lookup("k")),
    n_(CarreauYasudaCoeffs_.lookup("n")),
    a_(CarreauYasudaCoeffs_.lookup("a")),
    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::CarreauYasuda::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    CarreauYasudaCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    CarreauYasudaCoeffs_.lookup("nu0") >> nu0_;
    CarreauYasudaCoeffs_.lookup("nuInf") >> nuInf_;
    CarreauYasudaCoeffs_.lookup("k") >> k_;
    CarreauYasudaCoeffs_.lookup("n") >> n_;
    CarreauYasudaCoeffs_.lookup("a") >> a_;

    return true;
}


// ************************************************************************* //
