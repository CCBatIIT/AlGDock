#ifndef BMAINRESIDUE_H_
#define BMAINRESIDUE_H_

#include "bMoleculeReader.hpp"
#include "bgeneral.hpp"
#include "server.hpp"

/*
#ifndef MAIN_RESIDUE_DEBUG_SPECIFIC
#define MAIN_RESIDUE_DEBUG_SPECIFIC 1
#endif
#ifndef MAIN_RESIDUE_DEBUG_LEVEL01
#define MAIN_RESIDUE_DEBUG_LEVEL01
#endif

#ifndef MAIN_RESIDUE_DEBUG_LEVEL02
#define MAIN_RESIDUE_DEBUG_LEVEL02
#endif
*/
//using namespace SimTK;

void mol_StructureChainsBuild (MolStructure *, int);

//==============================================================================
//                           CLASS MainResidue
//==============================================================================
/**
 * Main Residue Class. It represents the main compound.
 **/
class bMainResidue : public SimTK::Compound{
public:

  MolAtom *bMolAtomList;
  MolStructure *struc;
  MolModel *model;
  bool hasBuiltSystem;
  unsigned int natms;
  bSpecificAtom *bAtomList;
  //std::vector<bBond> bonds; // RESTORE
  unsigned int nbnds; // EU
  bBond *bonds; // EU
  std::string ictdF;
  TARGET_TYPE *PrmToAx_po;
  TARGET_TYPE *MMTkToPrm_po;

  bMainResidue(
    SimTK::DuMMForceFieldSubsystem &dumm,
    unsigned int natms,
    bSpecificAtom *bAtomList,
    unsigned int nbnds,
    //std::vector<bBond> bonds, // RESTORE
    bBond *bonds, // EU
    TARGET_TYPE **coords,
    TARGET_TYPE **indexMap,
    TARGET_TYPE *PrmToAx_po,
    TARGET_TYPE *MMTkToPrm_po,
    bool first_time=true,
    std::string ictdF="IC"
  );

  ~bMainResidue();
};

#endif //BMAINRESIDUE_H_
