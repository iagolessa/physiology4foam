/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "CauchyStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CauchyStress, 0);
    addToRunTimeSelectionTable(functionObject, CauchyStress, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::CauchyStress::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Wall shear component");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    os << endl;
}


void Foam::functionObjects::CauchyStress::calcShearStress
(
    const volSymmTensorField& Reff,
    volSymmTensorField& shearStress
)
{
    shearStress.dimensions().reset(Reff.dimensions());

    for (const label patchi : patchSet_)
    {
        symmTensorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

		ssp = Reffp;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CauchyStress::CauchyStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchSet_()
{
    read(dict);

    writeFileHeader(file());

    volSymmTensorField* CauchyStressPtr
    (
        new volSymmTensorField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor(sqr(dimLength)/sqr(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(CauchyStressPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CauchyStress::~CauchyStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CauchyStress::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::CauchyStress::execute()
{
    volSymmTensorField& CauchyStress =
        mesh_.lookupObjectRef<volSymmTensorField>(type());

    // Compressible
    {
        typedef compressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcShearStress(modelPtr->devRhoReff(), CauchyStress);
            return true;
        }
    }

    // Incompressible
    {
        typedef incompressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcShearStress(modelPtr->devReff(), CauchyStress);
            return true;
        }
    }

    FatalErrorInFunction
        << "Unable to find turbulence model in the "
        << "database" << exit(FatalError);

    return false;
}


bool Foam::functionObjects::CauchyStress::write()
{
    const volSymmTensorField& CauchyStress =
        obr_.lookupObject<volSymmTensorField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << CauchyStress.name() << endl;

    CauchyStress.write();

    const fvPatchList& patches = mesh_.boundary();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const symmTensorField& ssp = CauchyStress.boundaryField()[patchi];

        symmTensor minSsp = gMin(ssp);
        symmTensor maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }

    return true;
}


// ************************************************************************* //
