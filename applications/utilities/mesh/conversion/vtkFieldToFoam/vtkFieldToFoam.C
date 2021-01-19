/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    vtkUnstructuredToFoam

Group
    grpMeshConversionUtilities

Description
    Convert legacy VTK file (ascii) containing an unstructured grid
    to an OpenFOAM mesh without boundary information.

Note
    The .vtk format does not contain any boundary information.
    It is purely a description of the internal mesh.
    Not extensively tested.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "pointFields.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "vtkUnstructuredReader.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert vtkUnstructuredGrid to OpenFOAM mesh and write the cell "
        "in it to a OpenFOAM volScalarField."
    );

    argList::noParallel();
    argList::addArgument("vtk-file", "The input legacy ascii vtk file");
    argList::addArgument("fieldName", "Name of the cell field to convert");

    #include "setRootCase.H"
    #include "createTime.H"

    IFstream mshStream(args[1]);
    const word fieldName(args[2]);

    vtkUnstructuredReader reader(runTime, mshStream);

    // Convert to polyMesh
    polyMesh vtkPolyMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        std::move(reader.points()),
        reader.cells(),
        faceListList(),
        wordList(),
        wordList(),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< "Writing mesh ..." << endl << endl;

    vtkPolyMesh.removeFiles();
    vtkPolyMesh.write();

    // Read polyMesh to fvMesh
    Info<< nl << "Converting cell data to volScalarField ..." << endl;

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Getting field from cells
    const scalarField T =
        reader.cellData().lookupObject<scalarField>(fieldName);

    // Creating field thickness
    volScalarField Thickness
    (
        IOobject
        (
            fieldName,
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0,1,0,0,0,0,0),
        T,
        zeroGradientFvPatchScalarField::typeName
    );


    Info << "Writing field " << fieldName << " ..." << endl;
    Thickness.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
