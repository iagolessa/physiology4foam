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

Application
    wallShearStress

Description
    Calculates and reports wall shear stress components for all patches, for
    the specified times when using RAS turbulence models.

    Default behaviour assumes operating in incompressible mode.
    Use the -compressible option for compressible RAS cases or -liquid for
    weakly compressible liquids (used with the sonicLiquidFoam solver. In this
    case, the Newtonian model is assumed with laminar flow and the dynamic
    viscosity mst be provided).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"

#include "basicPsiThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& wallShearComponent
)
{
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> model
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    // Define vol. field Reff initialized by the deviatoric part
    // of the viscous stress tensor for the laminar model
    // in this case
    const volSymmTensorField Reff(model->devReff());

    // Wall unit normal
    surfaceVectorField nf = -mesh.Sf()/mesh.magSf();

    // Defining object wallTraction
    surfaceVectorField wallTraction = nf & (fvc::interpolate(Reff));

    // Compute the wall shear stress component (tangential)
    // on the boundary
    forAll(wallShearComponent.boundaryField(), patchI)
    {
        wallShearComponent.boundaryField()[patchI] =
        wallTraction.boundaryField()[patchI]
        -
        (
            wallTraction.boundaryField()[patchI] & nf.boundaryField()[patchI]
        )*nf.boundaryField()[patchI];
    }
}


void calcCompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& wallShearComponent
)
{
    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!rhoHeader.headerOk())
    {
        Info<< "    no rho field" << endl;
        return;
    }

    Info<< "Reading field rho\n" << endl;
    volScalarField rho(rhoHeader, mesh);

    #include "compressibleCreatePhi.H"

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    autoPtr<compressible::RASModel> model
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    const volSymmTensorField Reff(model->devRhoReff());

    // Wall unit normal
    surfaceVectorField nf = -mesh.Sf()/mesh.magSf();

    // Defining object wallTraction
    surfaceVectorField wallTraction = nf & (fvc::interpolate(Reff));

    // Compute the wall shear stress component (tangential)
    // on the boundary
    forAll(wallShearComponent.boundaryField(), patchI)
    {
        wallShearComponent.boundaryField()[patchI] =
        wallTraction.boundaryField()[patchI]
        -
        (
            wallTraction.boundaryField()[patchI] & nf.boundaryField()[patchI]
        )*nf.boundaryField()[patchI];
    }
}


void calcLiquidCompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& wallShearComponent
)
{
    // Read the dynamic viscosity
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar mu(transportProperties.lookup("mu"));

    // Compute the effective shear stress tensor for laminar incompressible
    // flow and the Newtonian model
    const volSymmTensorField Reff(mu*twoSymm(fvc::grad(U)));

    // Wall unit normal
    surfaceVectorField nf = -mesh.Sf()/mesh.magSf();

    // Defining object wallTraction
    surfaceVectorField wallTraction = nf & (fvc::interpolate(Reff));

    // Compute the wall shear stress component (tangential)
    // on the boundary
    forAll(wallShearComponent.boundaryField(), patchI)
    {
        wallShearComponent.boundaryField()[patchI] =
        wallTraction.boundaryField()[patchI]
        -
        (
            wallTraction.boundaryField()[patchI] & nf.boundaryField()[patchI]
        )*nf.boundaryField()[patchI];
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::validOptions.insert("compressible","");
    argList::validOptions.insert("liquid","");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    bool compressible = args.optionFound("compressible");
    bool liquid = args.optionFound("liquid");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);

            // Define wal shear component vector field with per-density units
            // as similarly done in the standard utility
            volVectorField wallShearComponent
            (
                IOobject
                (
                    "wallShearComponent",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector
                (
                    "wallShearComponent",
                    sqr(dimLength)/sqr(dimTime),
                    vector::zero
                )
            );


            if (!liquid)
            {
                if (compressible)
                {
                    calcCompressible
                    (
                        mesh,
                        runTime,
                        U,
                        wallShearComponent
                    );
                }
                else
                {
                    calcIncompressible
                    (
                        mesh,
                        runTime,
                        U,
                        wallShearComponent
                    );
                }
            }
            else
            {
                // Redefine the stress units to *traction* units
                wallShearComponent.dimensions().reset(
                    dimPressure
                );

                calcLiquidCompressible
                (
                    mesh,
                    runTime,
                    U,
                    wallShearComponent
                );

            }

            Info<< "Writing wall shear component to field "
                << wallShearComponent.name()
                << nl << endl;

            wallShearComponent.write();
        }
        else
        {
            Info<< "    no U field" << endl;
        }

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
