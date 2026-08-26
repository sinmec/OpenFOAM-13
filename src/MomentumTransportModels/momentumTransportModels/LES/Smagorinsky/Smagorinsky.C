/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "Smagorinsky.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> Smagorinsky<BasicMomentumTransportModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*this->Ck_*this->delta()*(dev(D) && D));

    return volScalarField::New
    (
        this->groupName("k"),
        sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
    );
}


template<class BasicMomentumTransportModel>
void Smagorinsky<BasicMomentumTransportModel>::correctNut()
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    // Viscosidade turbulenta do Smagorinsky
    volScalarField nuTurb = this->Ck_*this->delta()*sqrt(k);

    // Fração volumétrica da fase contínua (H2O) - passada pelo solver
    const volScalarField& alphaLiq = this->alpha_;

    // Limitadores para evitar divisão por zero e sigFpe
    volScalarField alphaLim = min(max(alphaLiq, scalar(0.05)), scalar(1.0));

    // Viscosidade molecular da fase contínua
    const volScalarField& nu = this->nu();

    // Correlação: nu' = nu * [(1/alpha)^2 - 1]
    volScalarField nuParticle = nu * (pow(alphaLim, -2) - 1.0);

    // nut_ armazena turbulenta + partículas
    // A classe base calcula nuEff() = nu() + nut_ automaticamente
    this->nut_ = nuTurb + nuParticle;

    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Smagorinsky<BasicMomentumTransportModel>::Smagorinsky
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    LESeddyViscosity<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    )
{
    // FORÇA a leitura dos coeficientes do dicionário
    this->read();
    
    Info<< ">>> Smagorinsky loaded with Ck = " << this->Ck_
        << ", Ce = " << this->Ce_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool Smagorinsky<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        // Tenta ler Ck e Ce do sub-dicionário SmagorinskyCoeffs
        const dictionary& coeffs = this->coeffDict();
        
        coeffs.readIfPresent("Ck", this->Ck_);
        coeffs.readIfPresent("Ce", this->Ce_);

        Info<< ">>> Smagorinsky read() updated Ck = " << this->Ck_
            << ", Ce = " << this->Ce_ << endl;

        return true;
    }

    return false;
}


template<class BasicMomentumTransportModel>
void Smagorinsky<BasicMomentumTransportModel>::correct()
{
    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //