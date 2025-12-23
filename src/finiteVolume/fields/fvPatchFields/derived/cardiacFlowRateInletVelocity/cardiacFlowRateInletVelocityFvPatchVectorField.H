/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::cardiacFlowRateInletVelocityFvPatchVectorField

Description
    Velocity inlet boundary condition either correcting the extrapolated
    velocity or creating a paraboloid velocity field normal to the patch
    adjusted to match the specified flow rate.

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s

    The condition is recommended for cylindrical inlets, as it generates a
    paraboloid inlet velocity field, which may introduce reverse flow if the
    inlet is not circular. Additionally, the informed flow rate can be scaled
    by a term (Di/Davg)^n, where Di is the diameter of the inlet and Davg is a
    reference diameter where the passed flow rate was measured; n is a power
    coefficient. The motivation for this is Murray's law of cardiovascular
    flows. It is recommended to be used for flows whose inlet flow was
    extracted from average populations, and a patient-specific flow rate should
    be used, intead.

Usage
    \table
        Property     | Description             | Required    | Default value
        massFlowRate | mass flow rate [kg/s]   | no          |
        volumetricFlowRate | volumetric flow rate [m3/s]| no |
        rho          | density field name      | no          | rho
        rhoInlet     | inlet density           | no          |
        extrapolateProfile | Extrapolate velocity profile | no | false
        scaledFlowRate | Scale flow rate | no | false
        avgDiameter | If scaled, reference diameter | no |
        powerCoeff | If scaled, power coeff. | no |
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                cardiacFlowRateInletVelocity;
        volumetricFlowRate  0.2;
        extrapolateProfile  yes;
        value               uniform (0 0 0);
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    <patchName>
    {
        type                cardiacFlowRateInletVelocity;
        massFlowRate        0.2;
        extrapolateProfile  yes;
        rho                 rho;
        rhoInlet            1.0;
        value               uniform (0 0 0);
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive into the domain (as an inlet)
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types
    Foam::flowRateOutletVelocityFvPatchVectorField

SourceFiles
    cardiacFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cardiacFlowRateInletVelocityFvPatchVectorField_H
#define cardiacFlowRateInletVelocityFvPatchVectorField_H

#include "foamCompatibility.H"
#include "fixedValueFvPatchFields.H"

#if defined(FOAMEXTEND)
    #include "interpolationTable.H"
#else
    #include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class cardiacFlowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class cardiacFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

#if defined(FOAMEXTEND)
        //- Time-varying inlet integral flow rate
        const interpolationTable<scalar> flowRate_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Sum areas.  Note: not correctly updated on mesh motion
        scalar gSumArea_;
#else
        //- Time-varying inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoInlet_;

        //- Is volumetric?
        bool volumetric_;
#endif

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Central point of patch
        const vector centre_;

        //- Use Murray's law scaled flow rate?
        bool scaledFlowRate_;

        //- If scaled flow rate, the average diameter
        scalar avgDiameter_;

        //- If scaled flow rate, the power coefficient (n in Murray's law)
        scalar powerCoeff_;

    // Private member functions

#if !defined(FOAMEXTEND)
        template<class RhoType>
        //- Update the patch values given the appropriate density type and value
        void updateValues(const RhoType& rho);
#endif

        // Possible TODO: implement a unified calculateMurrayScaling() here
        // to return the (Di/Davg)^n factor for both versions.

public:

   //- Runtime type information
   TypeName("cardiacFlowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        cardiacFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cardiacFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  cardiacFlowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        cardiacFlowRateInletVelocityFvPatchVectorField
        (
            const cardiacFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cardiacFlowRateInletVelocityFvPatchVectorField
        (
            const cardiacFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new cardiacFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cardiacFlowRateInletVelocityFvPatchVectorField
        (
            const cardiacFlowRateInletVelocityFvPatchVectorField&,
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
                new cardiacFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
