/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::functionObjects::CauchyStress

Group
    grpForcesFunctionObjects

Description
    Calculates and write the Cauchy stress tensor at wall patches as
    the volSymmTensorField field 'CauchyStress'.

    The symmetric tensor field is retrieved from the turbulence model.  All
    wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    CauchyStress1
    {
        type        CauchyStress;
        libs        ("libfieldFunctionObjects.so");
        ...
        patches     (".*Wall");
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: CauchyStress | yes        |
        resultName | Name of wall shear stress field | no  | \<function name\>
        patches  | list of patches to process | no         | all wall patches
        log      | Log to standard output     | no         | yes
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    CauchyStress.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_CauchyStress_H
#define functionObjects_CauchyStress_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class CauchyStress Declaration
\*---------------------------------------------------------------------------*/

class CauchyStress
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the shear-stress
        void calcShearStress
        (
            const volSymmTensorField& Reff,
            volSymmTensorField& shearStress
        );


private:

    // Private member functions

        //- No copy construct
        CauchyStress(const CauchyStress&) = delete;

        //- No copy assignment
        void operator=(const CauchyStress&) = delete;


public:

    //- Runtime type information
    TypeName("CauchyStress");


    // Constructors

        //- Construct from Time and dictionary
        CauchyStress
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~CauchyStress();


    // Member Functions

        //- Read the CauchyStress data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
