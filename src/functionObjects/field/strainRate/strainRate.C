/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "strainRate.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(strainRate, 0);
    addToRunTimeSelectionTable(functionObject, strainRate, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::strainRate::calc()
{
    if (foundObject<volVectorField>(fieldName_))
    {
        // Get cell volumes
        volScalarField V
        (
            IOobject
            (
                mesh_.V().name(),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(mesh_.V().dimensions(), Zero),
            calculatedFvPatchField<scalar>::typeName
        );

        V.ref() = mesh_.V();

        const volVectorField& U = lookupObject<volVectorField>(fieldName_);
        const tmp<volTensorField> tgradU(fvc::grad(U));
        const volTensorField& gradU = tgradU();

        // Definition of the strain-rate based on the modified second
        // invariant of the velocity vector
        volScalarField magSymmGradU = mag(symm(gradU));

        const auto& transportProperties =
            mesh_.lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        //const volScalarField cellVols = mesh_.V();

        // Friction velocity
        const volScalarField sqrUTau = nu*magSymmGradU;

        // Length scales
        const volScalarField lPlus = sqrt(sqrUTau)*pow(V, 1.0/3.0)/nu;

        // Time scales
        const volScalarField tPlus = nu/sqrUTau;

        // Print max and min
        scalar minLPlus = gMin(lPlus);
        scalar maxLPlus = gMax(lPlus);

        scalar minTPlus = gMin(tPlus);
        scalar maxTPlus = gMax(tPlus);

        if (Pstream::master())
        {
            Log << " l+ : min = " << minLPlus 
                << ", max = " << maxLPlus << nl;

            Log << " t+ : min = " << minTPlus 
                << ", max = " << maxTPlus << nl;
        }

        return store
        (
            resultName_,
            sqrt(2.0)*magSymmGradU
        );
    }

    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::strainRate::strainRate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName(typeName, fieldName_);
}


// ************************************************************************* //
