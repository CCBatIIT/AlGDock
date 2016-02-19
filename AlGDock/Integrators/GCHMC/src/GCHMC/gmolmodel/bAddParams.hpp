#ifndef __BADDPARAMS__
#define __BADDPARAMS__

#include <iostream>
#include <fstream>
#include <sstream>
#include "Simbody.h"
#include "Molmodel.h"
#include "bMoleculeReader.hpp"

/*
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif

#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif
*/

// Conveninet functions (Pointers have to be allocated first).
bool Type2atomType(char *Type, char *atomType, int ATOMTYPE_MAX_LEN);
bool Type2atomType(string Type, char *atomType, int ATOMTYPE_MAX_LEN);

//==============================================================================
//                           FUNCTION AddGaffParams
//==============================================================================
/**
 * Add gaff Parameters Function. It adds the parameters read from an Amber parameters
 * file to the DuMM force field.
 **/
void bAddGaffParams(
  SimTK::DuMMForceFieldSubsystem& dumm,
  const char *filename,
  int natms,
  bSpecificAtom *bAtomList,
  vector<bBond> bonds,
  string frcmodfn
);

#endif  //__BADDPARAMS__

