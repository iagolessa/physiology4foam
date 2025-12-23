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

\*---------------------------------------------------------------------------*/

#include "paraboloidFlowRateInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::
paraboloidFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(),
#if defined(FOAMEXTEND)
    phiName_("phi"),
    gSumArea_(gSum(p.magSf())),
#else
    rhoInlet_(0.0),
    volumetric_(false),
#endif
    rhoName_("rho"),
    centre_(),
    scaledFlowRate_(false),
    avgDiameter_(),
    powerCoeff_()
{}


Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::
paraboloidFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
#if defined(FOAMEXTEND)
    fixedValueFvPatchField<vector>(p, iF, dict),
    flowRate_(dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    gSumArea_(gSum(p.magSf())),
#else
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -VGREAT)),
    volumetric_(false),
#endif
    rhoName_("rho"),
    centre_
    (
        dict.lookupOrDefault<vector>("centre", vector::zero)
    ),
    scaledFlowRate_
    (
        dict.lookupOrDefault<bool>("scaledFlowRate", false)
    ),
    avgDiameter_
    (
        dict.lookupOrDefault<scalar>("avgDiameter", 0.0)
    ),
    powerCoeff_
    (
        dict.lookupOrDefault<scalar>("powerCoeff", 0.0)
    )
{
#if defined(OPENFOAM_COM)
    if (dict.found("volumetricFlowRate"))
    {
        volumetric_ = true;
        flowRate_ = Function1<scalar>::New("volumetricFlowRate", dict);
    }
    else if (dict.found("massFlowRate"))
    {
        volumetric_ = false;
        flowRate_ = Function1<scalar>::New("massFlowRate", dict);
        rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << nl
            << exit(FatalIOError);
    }

    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
#endif

    if (scaledFlowRate_ && avgDiameter_ <= 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "'avgDiameter' not found or negative value specified "
            << "for scaled flow rate" << nl
            << exit(FatalIOError);
    }

    if (scaledFlowRate_ && powerCoeff_ == 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "'powerCoeff' not found or negative value specified "
            << "for scaled flow rate" << nl
            << exit(FatalIOError);
    }
}


Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::
paraboloidFlowRateInletVelocityFvPatchVectorField
(
    const paraboloidFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
#if defined(FOAMEXTEND)
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    gSumArea_(gSum(p.magSf())),
#else
    flowRate_(ptf.flowRate_.clone()),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
#endif
    rhoName_(ptf.rhoName_),
    centre_(ptf.centre_),
    scaledFlowRate_(ptf.scaledFlowRate_),
    avgDiameter_(ptf.avgDiameter_),
    powerCoeff_(ptf.powerCoeff_)
{}


Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::
paraboloidFlowRateInletVelocityFvPatchVectorField
(
    const paraboloidFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
#if defined(FOAMEXTEND)
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    gSumArea_(gSum(ptf.patch().magSf())),
#else
    flowRate_(ptf.flowRate_.clone()),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
#endif
    rhoName_(ptf.rhoName_),
    centre_(ptf.centre_),
    scaledFlowRate_(ptf.scaledFlowRate_),
    avgDiameter_(ptf.avgDiameter_),
    powerCoeff_(ptf.powerCoeff_)
{}


Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::
paraboloidFlowRateInletVelocityFvPatchVectorField
(
    const paraboloidFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
#if defined(FOAMEXTEND)
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    gSumArea_(gSum(ptf.patch().magSf())),
#else
    flowRate_(ptf.flowRate_.clone()),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
#endif
    rhoName_(ptf.rhoName_),
    centre_(ptf.centre_),
    scaledFlowRate_(ptf.scaledFlowRate_),
    avgDiameter_(ptf.avgDiameter_),
    powerCoeff_(ptf.powerCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#if defined(OPENFOAM_COM)
template<class RhoType>
void Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    const scalar t = db().time().timeOutputValue();

    // Volumetric flow-rate
    scalar tFlowRate = flowRate_->value(t);

    const vectorField n(patch().nf());

	const scalar gSumArea(gSum(patch().magSf()));

	// Store square radius for efficiency
	// Patch is considered circular
	scalar radiusSqr = gSumArea/constant::mathematical::pi;

	// Radial coordinate
	const scalarField r(mag(patch().Cf() - centre_));

    if (scaledFlowRate_)
    {
        // Updated volumetric flow-rate
        tFlowRate *= pow(2.0*sqrt(radiusSqr)/avgDiameter_, powerCoeff_);
    }

	// Compute mean velocity
	const scalar avgU = -tFlowRate/gSumArea;

	// Volumetric flow-rate
	operator==(n*2.0*(avgU/rho)*(1.0 - (r*r)/(radiusSqr)));
}
#endif


void Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

#if defined(FOAMEXTEND)
	// Store square radius
	// Patch is considered circular
	scalar radiusSqr = gSumArea_/mathematicalConstant::pi;

	// Radial coordinate
	const scalarField r(mag(patch().Cf() - centre_));

	// Get flow rate at the instant
	scalar curFlowRate = flowRate_(this->db().time().timeOutputValue());

    if (scaledFlowRate_)
    {
        // Updated volumetric flow-rate
        curFlowRate *= pow(2.0*sqrt(radiusSqr)/avgDiameter_, powerCoeff_);
    }

	// Compute mean velocity
	scalar avgU = -curFlowRate/(gSumArea_ + SMALL);

	vectorField n = patch().nf();

	const surfaceScalarField& phi =
		db().lookupObject<surfaceScalarField>(phiName_);

	// Evaluate boundary condition
    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // Volumetric flow-rate
		operator==(n*2.0*avgU*(1.0 - (r*r)/(radiusSqr)));
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        // Mass flow-rate
		operator==(n*2.0*(avgU/rhop)*(1.0 - (r*r)/(radiusSqr)));
    }
    else
    {
        FatalErrorIn
        (
            "flowRateInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << nl << exit(FatalError);
    }
#else

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one());
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            updateValues(rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoInlet_);
        }
    }

#endif

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::paraboloidFlowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

#if defined(FOAMEXTEND)
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }

	os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
	os.writeKeyword("scaledFlowRate") << scaledFlowRate_ << token::END_STATEMENT << nl;
	os.writeKeyword("avgDiameter") << avgDiameter_ << token::END_STATEMENT << nl;
	os.writeKeyword("powerCoeff") << powerCoeff_ << token::END_STATEMENT << nl;

	flowRate_.write(os);
#else
    flowRate_->writeData(os);
    if (!volumetric_)
    {
        os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
        os.writeEntryIfDifferent<scalar>("rhoInlet", -VGREAT, rhoInlet_);
    }

    os.writeEntry("centre", centre_);
    os.writeEntry("scaledFlowRate", scaledFlowRate_);
    os.writeEntry("avgDiameter", avgDiameter_);
    os.writeEntry("powerCoeff", powerCoeff_);
#endif

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       paraboloidFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
