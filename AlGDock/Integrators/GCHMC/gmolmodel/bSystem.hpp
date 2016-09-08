#ifndef BSYSTEM_H_
#define BSYSTEM_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include <thread>
#include <array>
//#include <random>
#include <math.h>

#include <Eigen/Dense>
//#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
//using Eigen::MatrixXd;

#include "Simbody.h"
#include "Molmodel.h"

// From MMTK:
#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include <boost/timer.hpp>

//#include "/home/lspirido/Installers/armadillo-6.400.3/include/armadillo.hpp"

/*
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif
#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif
#ifndef DEBUG_WRITEDATA
#define DEBUG_WRITEDATA
#endif
*/
/*
#ifndef DEBUG_WRITEALLPDBS
#define DEBUG_WRITEALLPDBS
#endif
#ifndef DEBUG_WRITEPDBS
#define DEBUG_WRITEPDBS
#endif
*/
/*
#ifndef DEBUG_WRITEPFPDB
#define DEBUG_WRITEPFPDB
#endif

#ifndef DEBUG_CONF
#define DEBUG_CONF
#endif
#ifndef DEBUG_SPECIFIC
#define DEBUG_SPECIFIC
#endif
*/
/*
#ifndef DEBUG_ENERGY
#define DEBUG_ENERGY
#endif
*/
/*
#ifndef DEBUG_TIME
#define DEBUG_TIME
#endif
*/

#include "bMoleculeReader.hpp"
#include "bAddParams.hpp"
#include "server.hpp"
#include "bMainResidue.hpp"
#include "bArgParser.hpp"

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix);

void writePdb(SimTK::PdbStructure pdb, const char *FN);

class SymSystem;

//==============================================================================
//                           CLASS GridForce
//==============================================================================
/**
 * External Custom Forces Class. It sends the cartesian coordinates calculated
 * by Compound and applies the forces passed by MMTK through PyArrayObjects 
 **/
class GridForce : public SimTK::Force::Custom::Implementation {
 public:
  SimTK::CompoundSystem *compoundSystem;
  TARGET_TYPE **indexMap;
  TARGET_TYPE *PrmToAx_po;
  TARGET_TYPE *MMTkToPrm_po;
  TARGET_TYPE **coords;
  TARGET_TYPE **vels;
  TARGET_TYPE **grads;
  int *fassno;
  int *flag;
  // From MMTK:
  PyFFEvaluatorObject *evaluator;
  energy_data *p_energy_po;
  PyArrayObject *configuration;
  PyUniverseSpecObject *universe_spec;
  TARGET_TYPE *shm;
  SymSystem *Caller;

  GridForce(SimTK::CompoundSystem *compoundSystem, SimTK::SimbodyMatterSubsystem& matter
            , TARGET_TYPE **indexMap, TARGET_TYPE *PrmToAx_po, TARGET_TYPE *MMTkToPrm_po
            , TARGET_TYPE **coords, TARGET_TYPE **vels, TARGET_TYPE **grads
            , int *fassno
            // From MMTK:
            , PyFFEvaluatorObject *evaluator
            , energy_data *p_energy_po
            , PyArrayObject *configuration
            , PyUniverseSpecObject *universe_spec
            , TARGET_TYPE *shm
            , SymSystem *Caller
            );

  void calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
    SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const;

  SimTK::Real calcPotentialEnergy(const SimTK::State& state) const;

  bool dependsOnlyOnPositions() const;

 private:
  SimTK::SimbodyMatterSubsystem& matter;
};

/////////////////////////////
//////// MIDVV INTEGRATOR ///
/////////////////////////////
class MidVVIntegratorRep;
//#include "/home/lspirido/Installers/simbody/simbody-Simbody-3.0/SimTKmath/Integrators/src/AbstractIntegratorRep.h"
#include "SimTKmath/Integrators/src/AbstractIntegratorRep.h"
#include "MidVVIntegrator.hpp"
/////////////////////////////
//////// MIDVV INTEGRATOR ///
/////////////////////////////

//==============================================================================
//                           CLASS SymSystem
//==============================================================================
/**
 *  Symbols System Class. Manages Molmodel-gMolmodel-MMTK communication.
 **/
class SymSystem{
 public:
  SimTK::CompoundSystem *system;
  SimTK::SimbodyMatterSubsystem *matter;
  SimTK::GeneralForceSubsystem *forces;
  SimTK::Force::Custom *ExtForce;
  SimTK::DecorationSubsystem *decorations;
  SimTK::DuMMForceFieldSubsystem *forceField;
  bMoleculeReader *mr;  // local
  bMainResidue *lig1;  // local
  SimTK::Visualizer *viz;
  #ifdef NOSETHERMOS
  SimTK::NoseHooverThermostat *thermo;
  #endif
  #ifdef VELSTHERMOS
  SimTK::VelocityRescalingThermostat *vthermo;
  #endif
  SimTK::Real startT;
  //RungeKuttaMersonIntegrator *integ;
  //VerletIntegrator *integ;
  MidVVIntegrator *integ;
  SimTK::TimeStepper *ts;
  TARGET_TYPE *PrmToAx_po;
  TARGET_TYPE *MMTkToPrm_po;
  string mol2F, rbF, gaffF, frcmodF, ictdF;

  TARGET_TYPE **coords;
  TARGET_TYPE **vels;
  TARGET_TYPE **inivels;
  TARGET_TYPE **indexMap;
  int **_indexMap;
  TARGET_TYPE **grads;
  int arrays_cut;

  int *passno;
  int *vassno;
  int *fassno;
  int *sassno;
  TARGET_TYPE *sysTimestep;
  TARGET_TYPE **QVector;
  SimTK::Transform *TVector;
  int **mbxTreeMat;    // tree representing the bonding
  SimTK::Real *branchMassVec; // branch masses self body included


  // From MMTK:
  PyFFEvaluatorObject *evaluator;
  energy_data *p_energy_po;
  PyArrayObject *configuration;
  PyUniverseSpecObject *universe_spec;
  TARGET_TYPE *shm;

  double **sysRetConfsPois;
  double *sysRetPotEsPoi;
  double *sysAccs;

  long unsigned int *pyseed;
  int massMatNumOpt;
  int metroFixmanOpt;
  
  SymSystem(
    string mol2F, string rbF, string gaffF, string frcmodF,
    string ictdF, TARGET_TYPE *PrmToAx_po, TARGET_TYPE *MMTkToPrm_po,
    // From MMTK:
    PyFFEvaluatorObject *evaluator,
    energy_data *p_energy_po,
    PyArrayObject *configuration,
    PyUniverseSpecObject *universe_spec,
    TARGET_TYPE *shm
  );

  void InitSimulation(
    TARGET_TYPE **coords,
    TARGET_TYPE **vels,
    TARGET_TYPE **inivels,
    TARGET_TYPE **indexMap,
    TARGET_TYPE **grads,
    TARGET_TYPE extTimestep,
    bool first_time
  );

  // Manages the TimeStepper actions
  void Advance(int nosteps);

  ~SymSystem(); // destructor
};

#endif /*BSYSTEM_H_*/
