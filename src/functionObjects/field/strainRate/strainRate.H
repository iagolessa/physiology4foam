/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::strainRate

Group
    grpFieldFunctionObjects

Description
    Computes the strain rate, the magnitude of symetric part of the velocity
    gradient.

    \f[
        S = \sqrt 2 \times symm \nabla \vec U
    \f]

    where
    \vartable
        S             | Strain rate         [1/s]
        \vec U        | Velocity vector     [m/s]
    \endvartable

    Operands:
    \table
      Operand        | Type           | Location
      input          | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volVectorField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    strainRate1
    {
        // Mandatory entries (unmodifiable)
        type        strainRate;
        libs        (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: strainRate               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func strainRate
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::strainRate

SourceFiles
    strainRate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_strainRate_H
#define functionObjects_strainRate_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class strainRate Declaration
\*---------------------------------------------------------------------------*/

class strainRate
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the strainRate field and return true if successful
        virtual bool calc();

public:

    //- Runtime type information
    TypeName("strainRate");


    // Constructors

        //- Construct from Time and dictionary
        strainRate
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        strainRate(const strainRate&) = delete;

        //- No copy assignment
        void operator=(const strainRate&) = delete;


    //- Destructor
    virtual ~strainRate() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
