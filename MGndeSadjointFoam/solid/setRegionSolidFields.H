const fvMesh& mesh = solidRegions[i]; //

//solidThermo& thermo = thermoSolid[i]; /thermoSolid

   

//tmp<volScalarField> trho = thermo.rho(); //te thermsolidThermrh
//const volScalarField& rho = trho(); /trh

//

//volScalarField& e = thermo.he(); //Given the enthalpy of the solid region


const Foam::fvModels& fvModels = fvModelsSolid[i];
Foam::fvConstraints& fvConstraints = fvConstraintsSolid[i];

#include "checkRadiationModel.H"

solidNoLoopControl& pimple = pimples.solid(i);
