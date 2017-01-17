#ifndef __OBCWRAPPER_H
#define __OBCWRAPPER_H

#ifdef __cplusplus

extern "C" {
#endif

typedef struct ObcParameters ObcParameters;
typedef struct ReferenceObc ReferenceObc;
// typedef double vector3[3];

ObcParameters* newObcParameters(int numParticles,
                                double strength,
                                const double* charges,
                                const double* atomicRadii,
                                const double* scaleFactors);

void setNumberOfAtoms(ObcParameters* obcParameters, int numAtoms);

void obcParameterReport(ObcParameters* obcParameters);
  
ReferenceObc* newReferenceObc(ObcParameters* obcParameters);

double computeBornEnergy(ReferenceObc* self,
                         ObcParameters* obcParameters,
                         vector3* atomCoordinates);

double computeBornEnergyForces(ReferenceObc* self,
                               ObcParameters* obcParameters,
                               vector3* atomCoordinates,
                               vector3* forces);

void deleteReferenceObc(ReferenceObc* self);

#ifdef __cplusplus
}
#endif
#endif
