//@file This is included in bSysSystem.cpp
class MidVVIntegratorRep;
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
//#include "/home/lspirido/Installers/simbody/simbody-Simbody-3.0/SimTKmath/Integrators/src/AbstractIntegratorRep.h"
#include "SimTKmath/Integrators/src/AbstractIntegratorRep.h"

//==============================================================================
//                           CLASS MidVVIntegrator
//==============================================================================
/**
 * Hamiltonian Monte Carlo Integrator Class with Verlet type integration between steps.
 **/
class MidVVIntegratorRep : public SimTK::AbstractIntegratorRep {
 public:
  TARGET_TYPE *PrmToAx_po;
  TARGET_TYPE *MMTkToPrm_po;
  TARGET_TYPE **coords;
  TARGET_TYPE **vels;
  TARGET_TYPE **grads;
  TARGET_TYPE *shm;
  int arrays_cut;
  unsigned int natms;
  long int trial;
  unsigned long int noSteps;
  int massMatNumOpt;
  int metroFixmanOpt;
  unsigned long int stepsPerTrial;
  unsigned long int ntrials;
  unsigned long int step;
  unsigned long int totStepsInCall;
  SimTK::CompoundSystem *compoundSystem;
  SymSystem *Caller;
  const SimTK::Compound& c;
  const SimTK::SimbodyMatterSubsystem& matter;

  TARGET_TYPE **QVector;
  SimTK::Transform *TVector;
  int *beginFlag;
  int *step0Flag;
  int *metroFlag;

  TARGET_TYPE Tb; // Bath temperature
  SimTK::Real kTb;
  SimTK::Real sqrtkTb;
  SimTK::Real invsqrtkTb;

  // Convenient intermediate data
  SimTK::Matrix Mi; // Mass matrix at the beg of MD
  SimTK::Matrix Mj; // Mass matrix at the end of MD
  SimTK::Real detsqrtMi; // det(sqrt(Mass matrix at the beg of MD)
  SimTK::Real detsqrtMj; // det(sqrt(Mass matrix at the beg of MD)
  SimTK::Real detMBATi; // det(BAT Mass matrix at the beg of MD)
  SimTK::Real detMBATj; // det(BAT Mass matrix at the end of MD)
  SimTK::Real detMi; // det(Mass matrix at the beg of MD)
  SimTK::Real detMj; // det(Mass matrix at the end of MD)
  SimTK::Real numDetMi; // det(Mass matrix at the beg of MD) numerically computed
  SimTK::Real numDetMj; // det(Mass matrix at the end of MD) numerically computed
  SimTK::Real numDetSqrtMi; // det(Mass matrix at the beg of MD) numerically computed
  SimTK::Real numDetSqrtMj; // det(Mass matrix at the end of MD) numerically computed
  SimTK::Real UFixi; // Fixman potential at the beg of MD
  SimTK::Real UFixj; // Fixman potential at the beg of MD
  SimTK::Real detMBATct; // ct part of the det(BAT Mass matrix at the beg of MD)

  TARGET_TYPE ke;
  TARGET_TYPE pe;
  TARGET_TYPE ini_pe;
  TARGET_TYPE timeToReach;

  SimTK::Random::Uniform unirand;
  SimTK::Random::Gaussian gaurand;
  boost::random::mt19937 rng;
  boost::random::uniform_real_distribution<> boostRealRand;

  SimTK::Real starttime;
  int trouble;

  std::map<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex> mbx2aIx;
  std::map<SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex> aIx2mbx;

 public:
  // SimTK::Compound -> Poss, Vels, Accs, Data
  // shm -> coords, vels, grads, ShmData
  void printData(const SimTK::Compound&, SimTK::State&);
  void printPoss(const SimTK::Compound&, SimTK::State&);
  void printVels(const SimTK::Compound&, SimTK::State&);
  void printAccs(const SimTK::Compound&, SimTK::State&);
  void printForces(const SimTK::Compound&, SimTK::State&);
  void printForcesNorms(const SimTK::Compound&, SimTK::State&);

  // Convenient intermediate functions
  SimTK::Real getUncorrAnglesSinesSq(const SimTK::Compound&, SimTK::State&, std::vector<int>, SimTK::Real *, int *); // Get uncorrelated angles sines squared
  SimTK::Real calcDetMBAT(const SimTK::Compound&, SimTK::State&); // Compute BAT mass matrix determinant (Jain et al, 2013)
  SimTK::Real calcDetMBATct(const SimTK::Compound&, SimTK::State&); // Compute BAT mass matrix ct part determinant (Jain et al, 2013)
  SimTK::Real calcDetM(const SimTK::Compound&, SimTK::State&); // Compute mass matrix determinant (Jain et al, 2013)
  SimTK::Real calcMAndDetM(const SimTK::Compound&, SimTK::State&); // Compute mass matrix and its determinant
  SimTK::Real calcUFix(const SimTK::Compound&, SimTK::State&); // Compute Fixman potential
  SimTK::Real calcUFixNum(const SimTK::Compound&, SimTK::State&); // Compute Fixman potential numerically
  SimTK::Real calcTotLinearVel(SimTK::State&);
  SimTK::Real calcKEFromAtoms(const SimTK::Compound&, SimTK::State&);
  SimTK::Real calcPEFromMMTK(void);

  // Flag operations
  void raiseBeginFlag(void){*this->beginFlag = 1;}
  void  dropBeginFlag(void){*this->beginFlag = 0;}
  int    getBeginFlag(void){return *this->beginFlag;}

  void raiseStep0Flag(void){*this->step0Flag = 1;}
  void  dropStep0Flag(void){*this->step0Flag = 0;}
  int    getStep0Flag(void){return *this->step0Flag;}

  void raiseMetroFlag(void){*this->metroFlag = 1;}
  void  dropMetroFlag(void){*this->metroFlag = 0;}
  int    getMetroFlag(void){return *this->metroFlag;}

  // Initialize velocities functions
  void initializeVelsFromRandU(const SimTK::Compound&, SimTK::State&, SimTK::Real);    // From random generalized vels
  void initializeVelsFromVRT(const SimTK::Compound&, SimTK::State&, SimTK::Real);    // From Velocity Rescaling Thermostate
  void initializeVelsFromAtoms(const SimTK::Compound&, SimTK::State&, SimTK::Real);  // From atoms inertia moments
  void initializeVelsFromMobods(const SimTK::Compound&, SimTK::State&, SimTK::Real); // From Mobilized bodies spatial inertia
  void initializeVelsFromJacobian(const SimTK::Compound&, SimTK::State&, SimTK::Real); // From Mobilized bodies spatial inertia
  void initializeVelsFromConst(const SimTK::Compound&, SimTK::State&, SimTK::SpatialVec *); // From const vector of vels

  // Assign conformations
  void setTVector(SimTK::State&); // Compound -> TVector
  void assignConfFromTVector(SimTK::State&); // XVector[] -> Compound
  void assignConfAndTVectorFromShm0(SimTK::State&); // shm -> Compound
  void assignConfAndTVectorFromShm1(SimTK::State&); // shm -> Compound; use matchDefault
  void assignConfAndTVectorFromShm2(SimTK::State&); // shm -> Compound; use Assembler
  SimTK::State& assignConfAndTVectorFromShm3(SimTK::State&); // shm -> Compound; use matchDefault and Assembler
  SimTK::State& assignConfAndTVectorFromShm3Opt(SimTK::State&); // shm -> Compound; use matchDefault and Assembler
  void writePossVelsToShm(const SimTK::Compound&, SimTK::State&); // Compound -> shm
  void writePossToShm(const SimTK::Compound&, SimTK::State&); // Compound -> shm
  void setMMTKConf(const SimTK::Compound&, SimTK::State&); // Compound -> MMTK

  int isTrouble(void);
  void setTrouble(void);
  void resetTrouble(void);

  void setTb(TARGET_TYPE);
  TARGET_TYPE getTb(void);
  void setKe(TARGET_TYPE);
  TARGET_TYPE getKe(void);
  void setPe(TARGET_TYPE);
  TARGET_TYPE getPe(void);
  void setIniPe(TARGET_TYPE);
  TARGET_TYPE getIniPe(void);

  void setTimeToReach(TARGET_TYPE);
  TARGET_TYPE getTimeToReach(void);

  void setTrial(long int);
  long int getTrial(void);
  void incrTrial(void);

  unsigned long int getNoSteps(void);
  void setNoSteps(unsigned long int);

  int getMassMatNumOpt(void);
  void setMassMatNumOpt(int);

  int getMetroFixmanOpt(void);
  void setMetroFixmanOpt(int);

  unsigned long int getStepsPerTrial(void);
  void setStepsPerTrial(unsigned long int);

  unsigned long int getNtrials(void);
  void setNtrials(unsigned long int);

  void resetStep(void);
  unsigned long int getStep(void);
  void incrStep(void);
  void incrStep(unsigned int);

  void resetTotStepsInCall(void);
  unsigned long int getTotStepsInCall(void);
  void incrTotStepsInCall(void);
  void incrTotStepsInCall(unsigned int);

  void metropolis(const SimTK::Compound&, SimTK::State&); // Metropolis function
  void try_finalize(const SimTK::Compound&, SimTK::State&, int, int); 

 public:
    MidVVIntegratorRep(SimTK::Integrator* handle, const SimTK::System& sys
                       , TARGET_TYPE *PrmToAx_po
                       , TARGET_TYPE *MMTkToPrm_po
                       , SimTK::CompoundSystem *compoundSystem
                       , SymSystem *Caller
                      );
 protected:
     bool attemptDAEStep
       (SimTK::Real t1, SimTK::Vector& yErrEst, int& errOrder, int& numIterations);
};


class SimTK_SIMMATH_EXPORT MidVVIntegrator : public SimTK::Integrator {
 public:

  void printData(const SimTK::Compound&, SimTK::State&);
  void printPoss(const SimTK::Compound&, SimTK::State&);
  void printVels(const SimTK::Compound&, SimTK::State&);
  void printAccs(const SimTK::Compound&, SimTK::State&);
  void printForces(const SimTK::Compound&, SimTK::State&); 
  void printForcesNorms(const SimTK::Compound&, SimTK::State&); 

  SimTK::Real getUncorrAnglesSinesSq(const SimTK::Compound&, SimTK::State&, std::vector<int>, SimTK::Real *, int *); // Get uncorrelated angles sines squared
  SimTK::Real calcDetMBATct(const SimTK::Compound&, SimTK::State&); // Compute BAT constant part of the mass matrix determinant (Jain et al, 2013)
  SimTK::Real calcDetMBAT(const SimTK::Compound&, SimTK::State&); // Compute BAT mass matrix determinant (Jain et al, 2013)
  SimTK::Real calcDetM(const SimTK::Compound&, SimTK::State&); // Compute mass matrix determinant (Jain et al, 2013)
  SimTK::Real calcMAndDetM(const SimTK::Compound&, SimTK::State&); // Compute mass matrix and its determinant
  SimTK::Real calcUFix(const SimTK::Compound&, SimTK::State&); // Compute Fixman potential
  SimTK::Real calcUFixNum(const SimTK::Compound&, SimTK::State&); // Compute Fixman potential numerically
  SimTK::Real calcTotLinearVel(SimTK::State&);
  SimTK::Real calcKEFromAtoms(const SimTK::Compound&, SimTK::State&);
  SimTK::Real calcPEFromMMTK(void);

  void raiseBeginFlag(void){((MidVVIntegratorRep *)rep)->raiseBeginFlag();}
  void  dropBeginFlag(void){((MidVVIntegratorRep *)rep)->dropBeginFlag();}
  int    getBeginFlag(void){return ((MidVVIntegratorRep *)rep)->getBeginFlag();}

  void raiseStep0Flag(void){((MidVVIntegratorRep *)rep)->raiseStep0Flag();}
  void  dropStep0Flag(void){((MidVVIntegratorRep *)rep)->dropStep0Flag();}
  int    getStep0Flag(void){return ((MidVVIntegratorRep *)rep)->getStep0Flag();}

  void raiseMetroFlag(void){((MidVVIntegratorRep *)rep)->raiseMetroFlag();}
  void  dropMetroFlag(void){((MidVVIntegratorRep *)rep)->dropMetroFlag();}
  int    getMetroFlag(void){return ((MidVVIntegratorRep *)rep)->getMetroFlag();}

  // Initialize velocities functions
  void initializeVelsFromRandU(const SimTK::Compound&, SimTK::State&, SimTK::Real);    // From random generalized vels
  void initializeVelsFromVRT(const SimTK::Compound&, SimTK::State&, SimTK::Real); // From Velocity Rescaling Thermostate
  void initializeVelsFromAtoms(const SimTK::Compound&, SimTK::State&, SimTK::Real); // From atoms inertia moments
  void initializeVelsFromMobods(const SimTK::Compound&, SimTK::State&, SimTK::Real); // From Mobilized bodies spatial inertia
  void initializeVelsFromJacobian(const SimTK::Compound&, SimTK::State&, SimTK::Real); // From Mobilized bodies spatial inertia
  void initializeVelsFromConst(const SimTK::Compound&, SimTK::State&, SimTK::SpatialVec *); // From constat vector of vels

  // Assign conformations
  void setTVector(SimTK::State&); // Compound -> TVector
  void assignConfFromTVector(SimTK::State&); // TVector[] -> Compound
  void assignConfAndTVectorFromShm0(SimTK::State&); // shm -> Compound
  void assignConfAndTVectorFromShm1(SimTK::State&); // shm -> Compound; use matchDefault
  void assignConfAndTVectorFromShm2(SimTK::State&); // shm -> Compound; use Assembler
  SimTK::State& assignConfAndTVectorFromShm3(SimTK::State&); // shm -> Compound; use matchDefault and Assembler
  SimTK::State& assignConfAndTVectorFromShm3Opt(SimTK::State&); // shm -> Compound; use matchDefault and Assembler
  void writePossVelsToShm(const SimTK::Compound&, SimTK::State&); // Compound -> shm
  void writePossToShm(const SimTK::Compound&, SimTK::State&); // Compound -> shm
  void setMMTKConf(const SimTK::Compound&, SimTK::State&); // Compound -> MMTK

  int isTrouble(void);
  void setTrouble(void);
  void resetTrouble(void);

  void setTb(TARGET_TYPE);
  TARGET_TYPE getTb(void);
  void setKe(TARGET_TYPE);
  TARGET_TYPE getKe(void);
  void setPe(TARGET_TYPE);
  TARGET_TYPE getPe(void);
  void setIniPe(TARGET_TYPE);
  TARGET_TYPE getIniPe(void);

  void setTimeToReach(TARGET_TYPE);
  TARGET_TYPE getTimeToReach(void);

  void setTrial(long int);
  long int getTrial(void);
  void incrTrial(void);

  unsigned long int getNoSteps(void);
  void setNoSteps(unsigned long int);

  int getMassMatNumOpt(void);
  void setMassMatNumOpt(int);

  int getMetroFixmanOpt(void);
  void setMetroFixmanOpt(int);

  unsigned long int getStepsPerTrial(void);
  void setStepsPerTrial(unsigned long int);

  unsigned long int getNtrials(void);
  void setNtrials(unsigned long int);

  void resetStep(void);
  unsigned long int getStep(void);
  void incrStep(void);
  void incrStep(unsigned int);

  void resetTotStepsInCall(void);
  unsigned long int getTotStepsInCall(void);
  void incrTotStepsInCall(void);
  void incrTotStepsInCall(unsigned int);

  void metropolis(const SimTK::Compound&, SimTK::State&); // Metropolis function

 public:
  /**
   * Create a VerletIntegrator for integrating a System with variable size steps.
   **/
  explicit MidVVIntegrator(const SimTK::System& sys
                           , TARGET_TYPE *PrmToAx_po
                           , TARGET_TYPE *MMTkToPrm_po
                           , SimTK::CompoundSystem *compoundSystem
                           , SymSystem *Caller
                          );
  /**
   * Create a VerletIntegrator for integrating a System with fixed size steps.
   **/
  MidVVIntegrator(const SimTK::System& sys, SimTK::Real stepSize
                  , TARGET_TYPE *PrmToAx_po
                  , TARGET_TYPE *MMTkToPrm_po
                  , SimTK::CompoundSystem *compoundSystem
                  , SymSystem *Caller
                  );
  ~MidVVIntegrator();
};



