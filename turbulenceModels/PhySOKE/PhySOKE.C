/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "PhySOKE.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PhySOKE, 0);
addToRunTimeSelectionTable(RASModel, PhySOKE, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void PhySOKE::correctNut()
{
    correctNonlinearStress(fvc::grad(U_));
}

void PhySOKE::correctNonlinearStress
(
    const volTensorField& gradU
)
{
    volVectorField U = this->U();
    volTensorField R(-skew(gradU));
    volSymmTensorField S(symm(gradU));

    volVectorField gradp = fvc::grad(this->mesh_.lookupObject<volScalarField>("p"));

    volScalarField d = wallDist(this->mesh_).y(); // 读取壁面距离场

    volScalarField q1 = 0.5*(sqr(tr(gradU)) - tr(((gradU) & (gradU))));
    volScalarField q2 = k_;
    volScalarField q3 = min(((sqrt(k_)*d)/(50*this->nu())), scalar(2));
    volScalarField q4 = U & gradp;
    volScalarField q5 = k_/epsilon_;
    volScalarField q6 = sqrt(gradp & gradp); // grad(p)分量影响收敛性，需要检查
    volScalarField q7 = mag((U*U) && gradU);
    volScalarField q8 = U & fvc::grad(k_);
    volScalarField q9 = mag(this->R());

    nut_ = 0.09*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    forAll(this->mesh_.C(), cell){       
        this->bBot_[cell].component(symmTensor::XX) = (q2[cell]*sqr(q3[cell])/sqr(1+2*q3[cell]));
        this->bBot_[cell].component(symmTensor::XY) = q5[cell]*q8[cell]/(3*q3[cell]+1);
        this->bBot_[cell].component(symmTensor::YY) = -q2[cell]/(2+q3[cell]);
        this->bBot_[cell].component(symmTensor::ZZ) = (q2[cell]*sqr(q3[cell])*sqr(q3[cell]-2));
    }

    this->aBot_ = this->alpha_*this->rho_*(this->bBot_);

    this->q2_ = q2;
    this->q3_ = q3;
    this->q5_ = q5;
    this->q8_ = q8;

    nonlinearStress_ = this->bBot_;
    this->Rnew = this->R() + this->bBot_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PhySOKE::PhySOKE
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    nonlinearEddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),


    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            coeffDict_,
            1.44
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.92
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    bBot_ // 各向异性雷诺应力的非线性部分(无量纲)
    (
        IOobject
        (
            "bBot",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("bBot", dimensionSet(0, 2, -2, 0, 0, 0, 0), Zero)
    ),
    aBot_ 
    (
        IOobject
        (
            "aBot",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("aBot", dimensionSet(0, 2, -2, 0, 0, 0, 0), Zero)
    ),
    q2_
    (
        IOobject
        (
            "q2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("q2", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0)
    ),
    q3_
    (
        IOobject
        (
            "q3",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("q3", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    ),
    q5_
    (
        IOobject
        (
            "q5",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("q5", dimensionSet(0, 0, 1, 0, 0, 0, 0), 0)
    ),
    q6_
    (
        IOobject
        (
            "q6",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("q6", dimensionSet(0, 1, -2, 0, 0, 0, 0), 0)
    ),
    q8_
    (
        IOobject
        (
            "q8",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("q8", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0)
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Rnew
    (
        IOobject
        (
            "Rnew",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("Rnew", dimensionSet(0, 2, -2, 0, 0, 0, 0), Zero)
    )

{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool PhySOKE::read()
{
    if (nonlinearEddyViscosity<incompressible::RASModel>::read())
    {
        Ceps1_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void PhySOKE::correct()
{
    if (!turbulence_)
    {
        return;
    }
    
    nonlinearEddyViscosity<incompressible::RASModel>::correct();

    tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    volScalarField G
    (
        GName(),
        (nut_*twoSymm(gradU) - nonlinearStress_) && gradU
    );


    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        Ceps1_*G*epsilon_/k_
       - fvm::Sp(Ceps2_*epsilon_/k_, epsilon_)
    );

    epsEqn.ref().relax();
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
