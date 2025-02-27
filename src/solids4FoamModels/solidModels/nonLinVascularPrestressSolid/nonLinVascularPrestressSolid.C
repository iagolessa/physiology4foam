/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "nonLinVascularPrestressSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinVascularPrestressSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, nonLinVascularPrestressSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, nonLinVascularPrestressSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void nonLinVascularPrestressSolid::predict()
{
    Info<< "Linear predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinVascularPrestressSolid::nonLinVascularPrestressSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    S0_
    (
        IOobject
        (
            "S0",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    zeroDisplacement_
    (
        IOobject
        (
            "zeroDisplacement",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zeroDisplacement",
            dimensionSet(0,1,0,0,0,0,0),
            vector::zero
        )
    )
{
    DisRequired();

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAMESIORFOUNDATION
            mesh().ddtScheme("ddt(" + D().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinVascularPrestressSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    if (predictor_)
    {
        predict();
    }

    int iCorr = 0;
#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector> solverPerfD;
    SolverPerformance<vector>::debug = 0;
#else
    lduSolverPerformance solverPerfD;
    blockLduMatrix::debug = 0;
#endif

    Info<< "Setting the displacement displacement to zero" << endl;

    // Force displacement to zero at the beginning of each time-step
    D() = zeroDisplacement_;

	// Update gradient of displacement
	mechanical().grad(zeroDisplacement_, gradD());

	// Update gradient of displacement increment
	gradDD() = gradD() - gradD().oldTime();

	// Total deformation gradient
	F_ = I + gradD().T();

	// Inverse of the deformation gradient
	Finv_ = inv(F_);

	// Jacobian of the deformation gradient
	J_ = det(F_);

	mechanical().correct(sigma());

    Info<< "Max. D " << max(D()).value() << " "
        << "Min. D " << min(D()).value() << endl;

    Info<< "Max. sigma " << max(mag(sigma())).value() << " "
        << "Min. sigma " << min(mag(sigma())).value() << endl;

    Info<< "F " << max(F_).value() << endl;

	// Define total stress
	volSymmTensorField totalSigma = sigma();

    Info<< "2nd PK stress before update:" << endl 
		<< "Max. S " << max(mag(S0_)).value() << " "
        << "Min. S " << min(mag(S0_)).value() << endl;

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();
		S0_.storePrevIter();
		
		// Update sigma to account for new S0
		totalSigma = sigma() + (1/J_)*symm(F_ & S0_ & F_.T());

        // Momentum equation total displacement total Lagrangian form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(J_*Finv_ & totalSigma, "div(sigma)")
          + rho()*g()
          + stabilisation().stabilisation(DD(), gradDD(), impK_)
        );

        // Under-relax the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Fixed or adaptive field under-relaxation
        relaxField(D(), iCorr);

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    while
    (
       !converged
        (
            iCorr,
#ifdef OPENFOAMESIORFOUNDATION
            mag(solverPerfD.initialResidual()),
            max
            (
                solverPerfD.nIterations()[0],
                max
                (
                    solverPerfD.nIterations()[1],
                    solverPerfD.nIterations()[2]
                )
            ),
#else
            solverPerfD.initialResidual(),
            solverPerfD.nIterations(),
#endif
            D()
        ) && ++iCorr < nCorr()
    );

    Info<< "Max. D " << max(mag(D())).value() << " "
        << "Min. D " << min(mag(D())).value() << endl;

    Info<< "2nd PK stress before update:" << endl 
		<< "Max. S " << max(mag(S0_)).value() << " "
        << "Min. S " << min(mag(S0_)).value() << endl;

    Info<< "F " << max(F_).value() << endl;

	// Update second Piola-Kirchhoff tensor
    S0_ += J_*symm(Finv_ & sigma() & Finv_.T());
	S0_.relax();

    Info<< "2nd PK stress after update:" << endl 
		<< "Max. S " << max(mag(S0_)).value() << " "
        << "Min. S " << min(mag(S0_)).value() << endl;

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

#ifndef OPENFOAMESIORFOUNDATION
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> nonLinVascularPrestressSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    vectorField nCurrent = Finv.T() & n;
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
