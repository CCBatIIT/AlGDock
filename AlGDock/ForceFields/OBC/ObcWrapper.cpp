#include "ObcParameters.h"
#include "ReferenceObc.h"
#include "ObcWrapper.h"
#include <iostream>

extern "C" {

ObcParameters* newObcParameters(int numAtoms,
                                double strength,
                                const double* charges,
                                const double* atomicRadii,
                                const double* scaleFactors) {
  
  std::vector<double> charges_v(charges, charges + numAtoms);
  std::vector<double> atomicRadii_v(atomicRadii, atomicRadii + numAtoms);
  std::vector<double> scaleFactors_v(scaleFactors, scaleFactors + numAtoms);
  
  ObcParameters* obcParameters = new ObcParameters(numAtoms,
    ObcParameters::ObcTypeII);
  
  obcParameters->setStrength(strength);
    
  obcParameters->setPartialCharges(charges_v);
  obcParameters->setAtomicRadii(atomicRadii_v);
  obcParameters->setScaledRadiusFactors(scaleFactors_v);

  obcParameters->setSolventDielectric(static_cast<double>(78.5));
  obcParameters->setSoluteDielectric(static_cast<double>(1.0));
  obcParameters->setPi4Asolv(4*M_PI*2.25936);
  obcParameters->setUseCutoff(static_cast<double>(1.5));
  
  return obcParameters;
}

void setNumberOfAtoms(ObcParameters* obcParameters, int numAtoms) {
  obcParameters->setNumberOfAtoms(numAtoms);
}

void obcParameterReport(ObcParameters* obcParameters) {
  const int numberOfAtoms = obcParameters->getNumberOfAtoms();
  const double dielectricOffset = obcParameters->getDielectricOffset();
  const double cutoffDistance = obcParameters->getCutoffDistance();
  const double soluteDielectric = obcParameters->getSoluteDielectric();
  const double solventDielectric = obcParameters->getSolventDielectric();
  const double strength = obcParameters->getStrength();
  std::cout << "numberOfAtoms = " << numberOfAtoms << std::endl;
  std::cout << "dielectricOffset = " << dielectricOffset << std::endl;
  std::cout << "cutoffDistance = " << cutoffDistance << std::endl;
  std::cout << "soluteDielectric = " << soluteDielectric << std::endl;
  std::cout << "solventDielectric = " << solventDielectric << std::endl;
  std::cout << "strength = " << strength << std::endl;
}

ReferenceObc* newReferenceObc(ObcParameters* obcParameters) {
  
  ReferenceObc *obc = new ReferenceObc(obcParameters);
  obc->setIncludeAceApproximation(true);

  return obc;
}

double computeBornEnergy(ReferenceObc* self, ObcParameters* obcParameters, const double* Igrid,vector3* atomCoordinates) {
  return self->computeBornEnergy(obcParameters, atomCoordinates, obcParameters->getPartialCharges(), Igrid);
}

double computeBornEnergyForces(ReferenceObc* self, ObcParameters* obcParameters, const double* Igrid, vector3* atomCoordinates, vector3* forces) {
  return self->computeBornEnergyForces(obcParameters, atomCoordinates, obcParameters->getPartialCharges(), Igrid, forces);
}

void deleteReferenceObc(ReferenceObc* self) {
//  delete self->getObcParameters();
  delete self;
}

} // extern "C"
