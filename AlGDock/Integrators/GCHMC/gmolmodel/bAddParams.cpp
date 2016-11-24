#include "bAddParams.hpp"


using namespace SimTK;
//using namespace std;


SimTK::Real bDihedral(SimTK::Vec3 pos0, SimTK::Vec3 pos1, SimTK::Vec3 pos2, SimTK::Vec3 pos3){
  SimTK::Vec3 diffs[3];
  SimTK::Vec3 normals[2];
  double dots[2];
  double psin, pcos, dih;

  diffs[0] = pos1 - pos0;
  diffs[1] = pos2 - pos1;
  diffs[2] = pos3 - pos2;

  normals[0] = diffs[0] % diffs[1];
  normals[1] = diffs[1] % diffs[2];
  dots[0] = (normals[0][0] * diffs[2][0]) + (normals[0][1] * diffs[2][1]) + (normals[0][2] * diffs[2][2]);
  dots[1] = (diffs[1][0] * diffs[1][0]) + (diffs[1][1] * diffs[1][1]) + (diffs[1][2] * diffs[1][2]);
  psin = dots[0] * std::sqrt(dots[1]);
  pcos = (normals[0][0] * normals[1][0]) + (normals[0][1] * normals[1][1]) + (normals[0][2] * normals[1][2]);

  return atan2(psin, pcos);
}


/***********************************
 * ================================
 *    ADD NEW PARAMETERS FUNCTION 
 * ================================
 ***********************************/

bool Type2atomType(char *Type, char *atomType, int ATOMTYPE_MAX_LEN){
  if((Type == NULL) || (atomType == NULL))
    return false;
  bZeroCharArray(atomType, ATOMTYPE_MAX_LEN);
  atomType[0] = Type[5];
  if(strlen(Type) >= 7){
    atomType[1] = Type[6];
  }
  else{
    atomType[1] = '\0';
  }
  return true;
}

bool Type2atomType(string Type, char *atomType, int ATOMTYPE_MAX_LEN){
  if(atomType == NULL)
    return false;
  bZeroCharArray(atomType, ATOMTYPE_MAX_LEN);
  atomType[0] = Type.at(5);
  if(Type.length() >= 7){
    atomType[1] = Type.at(6);
  }
  else{
    atomType[1] = '\0';
  }
  return true;
}

void bAddGaffParams(
    DuMMForceFieldSubsystem& dumm,
    const char *filename,
    int natms,
    bSpecificAtom *bAtomList,
    //vector<bBond> bonds, // RESTORE
    unsigned int nbnds,
    bBond *bonds, // EU
    string frcmodfn
  ){

  #ifdef DEBUG_PARAMS_LEVEL01
  std::cout<<"bAddGaffParams() BEGIN"<<std::endl<<std::flush;
  #endif

  /*Atomic masses*/
  std::map<string, float> am;
  /* //Molmodel masses
  am["x"]  = 16.00;    //dummy atom
  am["z"]  = 1.000;    //dummy atom
  am["c"]  = 12.01;
  am["h"]  = 1.008;
  am["f"]  = 19.00;
  am["cl"] = 35.45;
  am["br"] = 79.90;
  am["i"]  = 126.9;
  am["n"]  = 14.01;
  am["o"]  = 16.00;
  am["p"]  = 30.97;
  am["s"]  = 32.06;
  */
  // MMTK masses * //
  am["x"]  = 16.00;    //dummy atom
  am["z"]  = 1.000;    //dummy atom
  am["c"]  = am["C"] = 12.00;
  am["h"]  = am["H"] =  1.00;
  am["f"]  = am["F"]  = 19.00;
  am["cl"] = am["Cl"] = 35.45;
  am["br"] = am["Br"] = 80.00;
  am["i"]  = am["I"]  = 127.0;
  am["n"]  = am["N"]  = 14.00;
  am["o"]  = am["O"]  = 16.00;
  am["p"]  = am["P"]  = 31.00;
  am["s"]  = am["S"]  = 32.00;
  am["b2"]  = 14.00; // EU8/8

  // * Valences (TODO read from file)* //
  std::map<string, float> val;
  val["x"]  = 2;      // dummy atom
  val["z"]  = 1;      // dummy atom

  val["N3"] = 4;
  val["H"] = 1;
  val["CT"] = 4;
  val["HP"] = 1;
  val["HC"] = 1;
  val["C"] = 3;
  val["O"] = 1;
  val["N"] = 3;
  val["H1"] = 1;
  val["O2"] = 1;

  val["c"]  = 3;
  val["c1"] = 2;
  val["c2"] = 3;
  val["c3"] = 4;
  val["ca"] = 3;
  val["cp"] = 3;
  val["cq"] = 3;
  val["cc"] = 3;
  val["cd"] = 3;
  val["ce"] = 3;
  val["cf"] = 3;
  val["cg"] = 2;
  val["ch"] = 2;
  val["cx"] = 4;
  val["cy"] = 4;
  val["cu"] = 3;
  val["cv"] = 3;
  val["cz"] = 3;
  val["h1"] = 1;
  val["h2"] = 1;
  val["h3"] = 1;
  val["h4"] = 1;
  val["h5"] = 1;
  val["ha"] = 1;
  val["hc"] = 1;
  val["hn"] = 1;
  val["ho"] = 1;
  val["hp"] = 1;
  val["hs"] = 1;
  val["hw"] = 1;
  val["hx"] = 1;
  val["f"] =  1;
  val["cl"] = 1;
  val["br"] = 1;
  val["i"] =  1;
  val["n"] =  2;
  val["n1"] = 1;
  val["n2"] = 2;
  val["n3"] = 3;
  val["n4"] = 4;
  val["na"] = 2;
  val["nb"] = 2;
  val["nc"] = 2;
  val["nd"] = 2;
  val["ne"] = 2;
  val["nf"] = 2;
  val["nh"] = 3;
  val["no"] = 4;
  val["o"] =  1;
  val["o2"] = 1;
  val["oh"] = 2;
  val["os"] = 2;
  val["ow"] = 2;
  val["p2"] = 4;
  val["p3"] = 4;
  val["p4"] = 4;
  val["p5"] = 4;
  val["pb"] = 3;
  val["pc"] = 3;
  val["pd"] = 3;
  val["pe"] = 3;
  val["pf"] = 3;
  val["px"] = 4;
  val["py"] = 4;
  val["s"] =  1;
  val["s2"] = 2;
  val["s4"] = 3;
  val["s6"] = 4;
  val["sh"] = 4;
  val["ss"] = 4;
  val["sx"] = 3;
  val["sy"] = 4;
  val["cl"] = 1;
  val["b2"] = 2; // EU8/8
  
  /*Valence layer electrons*/
  std::map<string, float> vle;
  vle["x"]  = 8;  //dummy atom
  vle["z"]  = 1;  //dummy atom
  vle["c"]  = vle["C"] = 6;
  vle["h"]  = vle["H"] = 300; // for MMTK mass
  vle["f"]  = vle["F"] = 9;
  vle["cl"] = vle["Cl"] = 17;
  vle["br"] = vle["Br"] = 35;
  vle["i"]  = vle["I"] = 53;
  vle["n"]  = vle["N"] = 7;
  vle["o"]  = vle["O"] = 8;
  vle["p"]  = vle["P"] = 15;
  vle["s"]  = vle["S"] = 16;
  vle["b2"] = 2; // EU8/8
  vle["mg"] = vle["Mg"] = vle["MG"] = 12;
  vle["ep"] = 1; // Amber Extra Point Type

  /*Van der Waals radii - read from gaff MOD4 section*/
  std::map<string, float> vdw;
  std::map<string, int> Type2Ix;

  /*To check for duplicates in gaff*/
  std::vector<bBond> bondsSoFar;
  std::vector<bBond>::iterator bSFit1;
  std::vector<bBond>::iterator bSFit2;

  std::vector<intriad> anglesSoFar;
  std::vector<intriad>::iterator aSFit1;
  std::vector<intriad>::iterator aSFit2;

  

  int aIx;
  FILE *fpo;    // gaff file pointer
  FILE *frcmod;  // frcmod file pointer
  int atomTypeFlag = 0;
  int bondFlag = 0;
  int angleFlag = 0;
  int torsFlag = 0;
  //int imprFlag = 0;
  string line;

  int LINE_MAX_LEN = 500;
  int ATOMTYPE_MAX_LEN = 9;
  int PARAMDEF_MAX_LEN = 500;

  char line_c[LINE_MAX_LEN];
  bZeroCharArray(line_c, LINE_MAX_LEN);

  char atomType1[ATOMTYPE_MAX_LEN];
  char atomType2[ATOMTYPE_MAX_LEN];
  char atomType3[ATOMTYPE_MAX_LEN];
  char atomType4[ATOMTYPE_MAX_LEN];
  bZeroCharArray(atomType1, ATOMTYPE_MAX_LEN);
  bZeroCharArray(atomType2, ATOMTYPE_MAX_LEN);
  bZeroCharArray(atomType3, ATOMTYPE_MAX_LEN);
  bZeroCharArray(atomType4, ATOMTYPE_MAX_LEN);

  char paramDef[PARAMDEF_MAX_LEN];
  bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
 
  string buff1; string buff2;
  string buff3; string buff4;
  string buff;
  string buffvle;
  float vdwValue;
  float atomicMass;
  float well;
  float k1;    // force ct
  float equil1;  // equil value
  float k2;    // periodicity at tors
  float equil2;  // phase at tors
  float k3;
  float equil3;
  float perio1=0, perio2=0, perio3=0;
  int cnt;
  
  vector<float> torsK;
  vector<float> torsEquil;
  vector<float> torsPhase;
  //float phase1;

  /*
  //Levels for the dummies
  //Z        l3
  //|       /
  //X<-O<-l2-l3
  //   |    \
  //   H     l3
  //natms = Z
  //natms - 1 = X
  */
  set<string> l1types;  //level 1 = the hydrogen
  set<string>::iterator l1it;
  set<string> l2types;  //first molecule atom linked to H
  set<string>::iterator l2it;
  set<string> l3types;
  set<string>::iterator l3it;

  int l1no;
  vector<int> l2no; //many atoms on level 3
  vector<int> l3no; //many atoms on level 3
  for(l1no=0; l1no<natms; l1no++){
    if(toupper(bAtomList[l1no].elem) == 'O'){  //target atom
      l1types.insert(bAtomList[l1no].fftype);
      break;
    }
  }

  //2nd level
  for(unsigned int k=0; k<nbnds; k++){
    if(l1no == bonds[k].i){
      if(strcmp(bAtomList[bonds[k].j-1].fftype, "gaff_x")){
        l2no.push_back(bonds[k].j);
        l2types.insert(bAtomList[bonds[k].j-1].fftype);
      }
    }
    else if(l1no == bonds[k].j){
      if(strcmp(bAtomList[bonds[k].i-1].fftype, "gaff_x")){
        l2no.push_back(bonds[k].i);
        l2types.insert(bAtomList[bonds[k].i-1].fftype);
      }
    }
  }
  ///////////////////////////////////

  //3rd level
  for(unsigned int k=0; k<nbnds; k++){
    for(unsigned int i=0; i<l2no.size(); i++){
        if(l2no[i] == bonds[k].i){
          if(bonds[k].j != l1no){
            l3no.push_back(bonds[k].j);
            l3types.insert(bAtomList[bonds[k].j-1].fftype);
          }
        }
        else if(l2no[i] == bonds[k].j){
          if(bonds[k].i != l1no){
            l3no.push_back(bonds[k].i);
            l3types.insert(bAtomList[bonds[k].i-1].fftype);
          }
        }
    }
  }
  ///////////////////////////////////


  frcmod = fopen(frcmodfn.c_str(), "r");
  fpo = fopen(filename, "r");  // for gaff

  /*Jump to MOD4 section*/
  while(fgets(line_c, 200, fpo)){
    line = line_c;
    if(line.at(0) == 'M'){
      if(line.length() >= 4){
        if((line.at(1) == 'O') && (line.at(2) == 'D') &&
          (line.at(3) == '4')){
          std::cout<<"Found MOD4"<<std::endl<<std::flush;
          break;
        }
      }
    }
  }

  /*READ VDW VALUES*/
  while(fgets(line_c, 200, fpo)){
    line = line_c;
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      break;
    }
    sscanf(line_c, "%s%f", &atomType1, &vdwValue);
    buff1 = atomType1;
    vdw.insert (std::pair<string, float>(buff1, vdwValue));
  }
  vdw.insert(std::pair<string, float>("x", 1.0));    // dummy atom
  vdw.insert(std::pair<string, float>("z", 1.0));    // dummy atom
  vdw.insert(std::pair<string, float>("b2", 1.0));    // EU8/8

  /*Skip the first line*/
  rewind(fpo);
  fgets(line_c, 200, fpo);
  atomTypeFlag = 1;

  /*READ ATOM TYPES*/
  /*First deal with the dummy atoms*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams ATOM TYPES defineAtomClass - dummy atoms"<<std::endl<<std::flush;
  #endif

  aIx = dumm.getNextUnusedAtomClassIndex();
  Type2Ix.insert (std::pair<string, float>("gaff_b2", aIx)); // EU8/8
  dumm.defineAtomClass(
    (DuMM::AtomClassIndex)aIx,
    "gaff_b2",
    vle["b2"],
    val["b2"],
    vdw["b2"]/10,
    0.1  //random
  );

  aIx = dumm.getNextUnusedAtomClassIndex();
  Type2Ix.insert (std::pair<string, float>("gaff_x", aIx));
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"defineAtomClass gaff_x "<<" vle "<<vle["x"]<<" val "<<val["x"]<<" vdw "
    <<vdw["x"]/10<<" well "<<std::endl<<std::flush;
  #endif
  dumm.defineAtomClass(
    (DuMM::AtomClassIndex)aIx,
    "gaff_x",
    vle["x"],
    val["x"],
    vdw["x"]/10,
    0.1  //random
  );

  aIx = dumm.getNextUnusedAtomClassIndex();
  Type2Ix.insert (std::pair<string, float>("gaff_z", aIx));
  dumm.defineAtomClass(
    (DuMM::AtomClassIndex)aIx,
    "gaff_z",
    vle["z"],
    val["z"],
    vdw["z"]/10,
    0.1  //random
  );

  aIx = dumm.getNextUnusedAtomClassIndex();
  Type2Ix.insert (std::pair<string, float>("gaff_2C", aIx));
  dumm.defineAtomClass(
    (DuMM::AtomClassIndex)aIx,
    "gaff_2C",
    vle["c"],
    val["c3"],
    (1.9080)/10,
    0.1094  // de la c3
  );
  ////////////////////////////////////
  ////////////////////////////////////
  
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams ATOM TYPES defineAtomClass from gaff MOD4"<<std::endl<<std::flush;
  #endif

  while(fgets(line_c, 200, fpo)){
    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<line_c<<std::endl<<std::flush;
    #endif
    line = line_c;
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      atomTypeFlag = 0;
      bondFlag = 1;
      break;
    }
    aIx = dumm.getNextUnusedAtomClassIndex();
    sscanf(line_c, "%s%f%f", &atomType1, &atomicMass, &well);
    buff1 = "gaff_";
    buff1 += atomType1;
    if((tolower(atomType1[0]) == 'c') && (tolower(atomType1[1]) == 'l')){
      buffvle = "cl";
    }
    else if ((tolower(atomType1[0]) == 'b') && (tolower(atomType1[1]) == 'r')){
      buffvle = "br";
    }
    else if ((tolower(atomType1[0]) == 'm') && (tolower(atomType1[1]) == 'g')){
      buffvle = "mg";
    }
    else if ((tolower(atomType1[0]) == 'e') && (tolower(atomType1[1]) == 'p')){
      buffvle = "ep";
    }
    else{
      buffvle = atomType1[0]; 
      std::transform(buffvle.begin(), buffvle.end(), buffvle.begin(), ::tolower); // LS
    }
    Type2Ix.insert (std::pair<string, float>(buff1, aIx));
    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"defineAtomClass "<<buff1.c_str()<<" buffvle "<<buffvle<<" vle "<<vle[buffvle]<<" val "<<val[atomType1]<<" vdw "
      <<vdw[buff1.substr(5,buff1.length()-5)]/10<<" well "<<well<<" aIx "<<aIx<<std::endl<<std::flush;
    #endif
    dumm.defineAtomClass(
      (DuMM::AtomClassIndex)aIx,
      buff1.c_str(),
      vle[buffvle],
      val[atomType1],
      vdw[buff1.substr(5,buff1.length()-5)]/10,
      well
    );
  }

  /*READ BOND PARAMS*/
  /*First deal with the dummy atoms*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams BONDS defineBondStretch dummy atoms"<<std::endl<<std::flush;
  #endif

  dumm.defineBondStretch(
    (DuMM::AtomClassIndex)(Type2Ix["gaff_z"]),
    (DuMM::AtomClassIndex)(Type2Ix["gaff_x"]),
    0.1,  //k1
    0.1    //equil1
  );

  for(l1it = l1types.begin(); l1it != l1types.end(); ++l1it){
    dumm.defineBondStretch(
      (DuMM::AtomClassIndex)(Type2Ix["gaff_x"]),
      (DuMM::AtomClassIndex)(Type2Ix[(*l1it).c_str()]),
      0.1,  //k1
      0.1    //equil1
    );
  }

  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams BONDS defineBondStretch gaff"<<std::endl<<std::flush;
  #endif

  fgets(line_c, 200, fpo);
  while(fgets(line_c, 200, fpo)){
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      bondFlag = 0;
      angleFlag = 1;
      break;
    }
    line = line_c;
    bSubstr(atomType1, line_c, 0, 2);
    if(atomType1[1] == ' '){atomType1[1] = '\0';}
    else{atomType1[2] = '\0';}
    bSubstr(atomType2, line_c, 3, 2);
    if(atomType2[1] == ' '){atomType2[1] = '\0';}
    else{atomType2[2] = '\0';}
    buff1 = "gaff_"; buff1 += atomType1;
    buff2 = "gaff_"; buff2 += atomType2;
    bSubstr(paramDef, line_c, 7,5);
    k1 = atof(paramDef);
    bSubstr(paramDef, line_c, 16,6);
    equil1 = atof(paramDef);

    // Check if bond is unique
    bBond thisBond(Type2Ix[buff1], Type2Ix[buff2]);
    bSFit1 = bondsSoFar.begin();
    bSFit2 = bondsSoFar.end();
    int found = 0;
    while(bSFit1 != bSFit2){
      if(bSFit1->isTheSameAs(&thisBond)){
        found = 1;
        break;
      }
      ++bSFit1;
    }

    //Add bond param
    equil1 = equil1/10.0;  // From angstroms to nanometers
    k1 = 418.4 * k1;    // From kcal to kJ
    if(!found){
      dumm.defineBondStretch(
        (DuMM::AtomClassIndex)(Type2Ix[buff1]),
        (DuMM::AtomClassIndex)(Type2Ix[buff2]),
        int(k1),
        equil1
        );
      bondsSoFar.push_back(thisBond);
    }
  }

  /*READ ANGLE PARAMS*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams ANGLES defineBondBend dummy atoms"<<std::endl<<std::flush;
  #endif

  /*Deal with dummy atoms first*/
  for(l1it = l1types.begin(); l1it != l1types.end(); ++l1it){
    dumm.defineBondBend(
      (DuMM::AtomClassIndex)(Type2Ix["gaff_z"]),
      (DuMM::AtomClassIndex)(Type2Ix["gaff_x"]),
      (DuMM::AtomClassIndex)(Type2Ix[(*l1it).c_str()]),
      0.1,
      90.0
    );
  }

  for(l2it = l2types.begin(); l2it != l2types.end(); ++l2it){
    for(l1it = l1types.begin(); l1it != l1types.end(); ++l1it){
      dumm.defineBondBend(
        (DuMM::AtomClassIndex)(Type2Ix["gaff_x"]),
        (DuMM::AtomClassIndex)(Type2Ix[(*l1it).c_str()]),
        (DuMM::AtomClassIndex)(Type2Ix[(*l2it).c_str()]),
        0.1,
        90.0
      );
    }
  }


  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams BONDS defineBondStretch frcmod"<<std::endl<<std::flush;
  #endif

  /*Read bonds from frcmod NEW */
  while(fgets(line_c, 200, frcmod)){
    line = line_c;
    buff1 = line.substr(0,4);
    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"bAddParams: ini "<<line_c<<std::endl;
    std::cout<<"bAddParams: buf "<<buff1<<std::endl;
    #endif
    if(buff1 == "BOND"){
      break;
    }
  }

  buff1 = "";
  while(fgets(line_c, 200, frcmod)){
    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"bAddParams: "<<line_c<<std::endl;
    #endif
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      break;
    }
    line = line_c;
    bSubstr(atomType1, line_c, 0, 2);
    if(atomType1[1] == ' '){atomType1[1] = '\0';}
    else{atomType1[2] = '\0';}
    bSubstr(atomType2, line_c, 3, 2);
    if(atomType2[1] == ' '){atomType2[1] = '\0';}
    else{atomType2[2] = '\0';}
    buff1 = "gaff_"; buff1 += atomType1;
    buff2 = "gaff_"; buff2 += atomType2;
    bSubstr(paramDef, line_c, 7, 6);
    k1 = atof(paramDef);
    bSubstr(paramDef, line_c, 16,5);
    equil1 = atof(paramDef);

    // Check if bond is unique
    bBond thisBond(Type2Ix[buff1], Type2Ix[buff2]);
    bSFit1 = bondsSoFar.begin();
    bSFit2 = bondsSoFar.end();
    int found = 0;
    while(bSFit1 != bSFit2){
      if(bSFit1->isTheSameAs(&thisBond)){
        found = 1;
        break;
      }
      ++bSFit1;
    }

    //Add bond param
    equil1 = equil1/10.0;  // From angstroms to nanometers
    k1 = 4.184 * k1;    // From kcal to kJ // EU8/8
    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"bAddParams: BOND "<<buff1<<' '<<buff2<<' '<<k1/418.4<<' '<<equil1*10<<std::endl;
    #endif
    if(!found){
      dumm.defineBondStretch(
        (DuMM::AtomClassIndex)(Type2Ix[buff1]),
        (DuMM::AtomClassIndex)(Type2Ix[buff2]),
        int(k1),
        equil1
        );
      bondsSoFar.push_back(thisBond);
      #ifdef DEBUG_PARAMS_LEVEL02
      std::cout<<"bAddParams: BondStretch "<<Type2Ix[buff1]<<' '<<Type2Ix[buff2]<<" added "<<k1<<' '<<equil1<<std::endl;
      #endif
    }
  }

  /*Read angles from frcmod*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams ANGLES defineBondBend frcmod"<<std::endl<<std::flush;
  #endif

  while(fgets(line_c, 200, frcmod)){
    line = line_c;
    buff1 = line.substr(0,5);
    if(buff1 == "ANGLE"){
      break;
    }
  }

  buff1 = "";
  while(fgets(line_c, 200, frcmod)){
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      break;
    }
    line = line_c;
    bSubstr(atomType1, line_c, 0, 2);
    if(atomType1[1] == ' '){atomType1[1] = '\0';}
    else{atomType1[2] = '\0';}
    bSubstr(atomType2, line_c, 3, 2);
    if(atomType2[1] == ' '){atomType2[1] = '\0';}
    else{atomType2[2] = '\0';}
    bSubstr(atomType3, line_c, 6, 2);
    if(atomType3[1] == ' '){atomType3[1] = '\0';}
    else{atomType3[2] = '\0';}
    buff1 = "gaff_"; buff1 += atomType1;
    buff2 = "gaff_"; buff2 += atomType2;
    buff3 = "gaff_"; buff3 += atomType3;
    bSubstr(paramDef, line_c, 10,7);
    k1 = atof(paramDef);
    bSubstr(paramDef, line_c, 22,7);
    equil1 = atof(paramDef);

    // Check if angle is unique
    intriad thisAngle(Type2Ix[buff1], Type2Ix[buff2], Type2Ix[buff3]);
    aSFit1 = anglesSoFar.begin();
    aSFit2 = anglesSoFar.end();
    int found = 0;
    while(aSFit1 != aSFit2){
      if(aSFit1->isTheSameAs(&thisAngle)){
        found = 1;
        break;
      }
      ++aSFit1;
    }

    //Add angle param
    k1 = 4.184 * k1;    // From kcal to kJ
    if(!found){
      #ifdef DEBUG_PARAMS_LEVEL02
      std::cout<<"bAddParams: Angle Bond Bend "<<buff1<<' '<<buff2<<' '<<buff3<<" "<<k1<<' '<<equil1<<std::endl;
      std::cout<<"bAddParams: Angle Bond Bend "<<Type2Ix[buff1]
        <<' '<<Type2Ix[buff2]<<' '<<Type2Ix[buff3]<<" "<<k1<<' '<<equil1<<std::endl;
      #endif
      dumm.defineBondBend(
        (DuMM::AtomClassIndex)(Type2Ix[buff1]),
        DuMM::AtomClassIndex(Type2Ix[buff2]),
        DuMM::AtomClassIndex(Type2Ix[buff3]),
        k1, 
        equil1
      );
      anglesSoFar.push_back(thisAngle);
    }
  }

  /*Read angles from gaff*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams ANGLES defineBondBend gaff"<<std::endl<<std::flush;
  #endif

  while(fgets(line_c, 200, fpo)){
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      angleFlag = 0;
      torsFlag = 1;
      break;
    }
    line = line_c;
    bSubstr(atomType1, line_c, 0, 2);
    if(atomType1[1] == ' '){atomType1[1] = '\0';}
    else{atomType1[2] = '\0';}
    bSubstr(atomType2, line_c, 3, 2);
    if(atomType2[1] == ' '){atomType2[1] = '\0';}
    else{atomType2[2] = '\0';}
    bSubstr(atomType3, line_c, 6, 2);
    if(atomType3[1] == ' '){atomType3[1] = '\0';}
    else{atomType3[2] = '\0';}
    buff1 = "gaff_"; buff1 += atomType1;
    buff2 = "gaff_"; buff2 += atomType2;
    buff3 = "gaff_"; buff3 += atomType3;
    bSubstr(paramDef, line_c, 10,7);
    k1 = atof(paramDef);
    bSubstr(paramDef, line_c, 22,7);
    equil1 = atof(paramDef);

    // Check if angle is unique
    intriad thisAngle(Type2Ix[buff1], Type2Ix[buff2], Type2Ix[buff3]);
    aSFit1 = anglesSoFar.begin();
    aSFit2 = anglesSoFar.end();
    int found = 0;
    while(aSFit1 != aSFit2){
      if(aSFit1->isTheSameAs(&thisAngle)){
        found = 1;
        break;
      }
      ++aSFit1;
    }

    //Add angle param

    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"bAddParams hop angles "<<buff1<<" "<<buff2<<" "<<buff3<<std::endl<<std::flush;
    std::cout<<"bAddParams hop angles "<<Type2Ix[buff1]<<" "<<Type2Ix[buff2]<<" "<<Type2Ix[buff3]<<std::endl<<std::flush;
    #endif

    k1 = 4.184 * k1;    // From kcal to kJ
    if(!found){
      dumm.defineBondBend(
        (DuMM::AtomClassIndex)(Type2Ix[buff1]),
        DuMM::AtomClassIndex(Type2Ix[buff2]),
        DuMM::AtomClassIndex(Type2Ix[buff3]),
        k1, 
        equil1
      );
      anglesSoFar.push_back(thisAngle);
    }
  }

  /*READ SPECIFIC TORSION PARAMS*/
  std::vector<string> spTaT1;
  std::vector<string>::iterator spit1;
  std::vector<string> spTaT2;
  std::vector<string>::iterator spit2;
  std::vector<string> spTaT3;
  std::vector<string>::iterator spit3;
  std::vector<string> spTaT4;
  std::vector<string>::iterator spit4;

  std::vector<string> uspTaT1;
  std::vector<string>::iterator uspit1;
  std::vector<string> uspTaT2;
  std::vector<string>::iterator uspit2;
  std::vector<string> uspTaT3;
  std::vector<string>::iterator uspit3;
  std::vector<string> uspTaT4;
  std::vector<string>::iterator uspit4;
  std::vector<float> uspTk1;
  std::vector<float> uspTequil1;
  char *atomType = new char[ATOMTYPE_MAX_LEN];


  /*Deal with dummy atoms*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams DIHEDRALS defineBondTorsion dummy atoms"<<std::endl<<std::flush;
  #endif

  for(l2it = l2types.begin(); l2it != l2types.end(); ++l2it){
    for(l1it = l1types.begin(); l1it != l1types.end(); ++l1it){
      dumm.defineBondTorsion(
        (DuMM::AtomClassIndex)(Type2Ix["gaff_z"]),
        (DuMM::AtomClassIndex)(Type2Ix["gaff_x"]),
        (DuMM::AtomClassIndex)(Type2Ix[(*l1it).c_str()]),
        (DuMM::AtomClassIndex)(Type2Ix[(*l2it).c_str()]),
        1,
        0.1,
        0.0
      );
      Type2atomType("gaff_z", atomType1, ATOMTYPE_MAX_LEN);
      Type2atomType("gaff_x", atomType2, ATOMTYPE_MAX_LEN);
      Type2atomType((*l1it).c_str(), atomType3, ATOMTYPE_MAX_LEN);
      Type2atomType((*l2it).c_str(), atomType4, ATOMTYPE_MAX_LEN);
      spTaT1.push_back(atomType1);
      spTaT2.push_back(atomType2);
      spTaT3.push_back(atomType3);
      spTaT4.push_back(atomType4);
    }
  }

  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams hop 1"<<std::endl<<std::flush;
  #endif


  for(l3it = l3types.begin(); l3it != l3types.end(); ++l3it){
    for(l2it = l2types.begin(); l2it != l2types.end(); ++l2it){
      if(strcmp((*l2it).c_str(), "gaff_ho")){  // DIRTY fix
        for(l1it = l1types.begin(); l1it != l1types.end(); ++l1it){
          dumm.defineBondTorsion(
            (DuMM::AtomClassIndex)(Type2Ix["gaff_x"]),
            (DuMM::AtomClassIndex)(Type2Ix[(*l1it).c_str()]),
            (DuMM::AtomClassIndex)(Type2Ix[(*l2it).c_str()]),
            (DuMM::AtomClassIndex)(Type2Ix[(*l3it).c_str()]),
            1,
            0.1,
            0.0
          );
          Type2atomType("gaff_x", atomType1, ATOMTYPE_MAX_LEN);
          Type2atomType((*l1it).c_str(), atomType2, ATOMTYPE_MAX_LEN);
          Type2atomType((*l2it).c_str(), atomType3, ATOMTYPE_MAX_LEN);
          Type2atomType((*l3it).c_str(), atomType4, ATOMTYPE_MAX_LEN);
          spTaT1.push_back(atomType1);
          spTaT2.push_back(atomType2);
          spTaT3.push_back(atomType3);
          spTaT4.push_back(atomType4);
        }
      }
    }
  }
  /////////////////////////
  
  /*Read torsions from frcmod*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams DIHEDRALS defineBondTorsion frcmod"<<std::endl<<std::flush;
  #endif

  string prevQuad = "";
  buff1 = "";
  atomType1[0] = '\0'; atomType1[1] = '\0';
  atomType2[0] = '\0'; atomType2[1] = '\0';
  atomType3[0] = '\0'; atomType3[1] = '\0';
  atomType4[0] = '\0'; atomType4[1] = '\0';

  while(fgets(line_c, 200, frcmod)){
    line = line_c;
    buff1 = line.substr(0,4);
    if(buff1 == "DIHE"){
      break;
    }
  }

  while(fgets(line_c, 200, frcmod)){
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      torsFlag = 0;
      break;
    }
    line = line_c;

    if(prevQuad == line.substr(0,11)){
      ++cnt;
      prevQuad = line.substr(0,11);
    }
    else{
      if((prevQuad.substr(0,1) != "X") && (atomType1[0] != '\0')){ // Specific Quads
        spTaT1.push_back(atomType1);
        spTaT2.push_back(atomType2);
        spTaT3.push_back(atomType3);
        spTaT4.push_back(atomType4);
    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<line;
    std::cout<<"bAddParams frcmod dihe "<<buff1<<" "<<buff2<<" "<<buff3<<" "<<buff4<<std::endl<<std::flush;
    std::cout<<"bAddParams frcmod dihe "<<Type2Ix[buff1]<<" "<<Type2Ix[buff2]<<" "<<Type2Ix[buff3]<<" "<<Type2Ix[buff4]<<std::endl<<std::flush;
    #endif

        if(cnt == 1){
          dumm.defineBondTorsion(
            (DuMM::AtomClassIndex)(Type2Ix[buff1]),
            (DuMM::AtomClassIndex)(Type2Ix[buff2]),
            (DuMM::AtomClassIndex)(Type2Ix[buff3]),
            (DuMM::AtomClassIndex)(Type2Ix[buff4]),
            1, k1, equil1
            );
        }
        else if(cnt == 2){
          dumm.defineBondTorsion(
            (DuMM::AtomClassIndex)(Type2Ix[buff1]),
            (DuMM::AtomClassIndex)(Type2Ix[buff2]),
            (DuMM::AtomClassIndex)(Type2Ix[buff3]),
            (DuMM::AtomClassIndex)(Type2Ix[buff4]),
            1, k1, equil1,
            2, k2, equil2
            );
        }
        else if(cnt >= 3){
          dumm.defineBondTorsion(
            (DuMM::AtomClassIndex)(Type2Ix[buff1]),
            (DuMM::AtomClassIndex)(Type2Ix[buff2]),
            (DuMM::AtomClassIndex)(Type2Ix[buff3]),
            (DuMM::AtomClassIndex)(Type2Ix[buff4]),
            1, k1, equil1,
            2, k2, equil2,
            3, k3, equil3
            );
        }
      } // if(line.at(0) != 'X')
      cnt = 1;
      prevQuad = line.substr(0,11);
    }

    // Read atom types
    bSubstr(atomType1, line_c, 0, 2);
    if(atomType1[1] == ' '){atomType1[1] = '\0';}
    else{atomType1[2] = '\0';}
    bSubstr(atomType2, line_c, 3, 2);
    if(atomType2[1] == ' '){atomType2[1] = '\0';}
    else{atomType2[2] = '\0';}
    bSubstr(atomType3, line_c, 6, 2);
    if(atomType3[1] == ' '){atomType3[1] = '\0';}
    else{atomType3[2] = '\0';}
    bSubstr(atomType4, line_c, 9, 2);
    if(atomType4[1] == ' '){atomType4[1] = '\0';}
    else{atomType4[2] = '\0';}
    buff1 = "gaff_"; buff1 += atomType1;
    buff2 = "gaff_"; buff2 += atomType2;
    buff3 = "gaff_"; buff3 += atomType3;
    buff4 = "gaff_"; buff4 += atomType4;

    if(cnt == 1){
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 14,1); perio1 = atof(paramDef);
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 18,5); k1 = atof(paramDef);
      k1 = k1/perio1;
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 30,6); equil1 = atof(paramDef);
    }
    else if(cnt == 2){
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 14,1); perio2 = atof(paramDef);
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 18,5); k2 = atof(paramDef);
      k2 = k2/perio2;
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 30,6); equil2 = atof(paramDef);
    }
    else if(cnt == 3){
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 14,1); perio3 = atof(paramDef);
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 18,5); k3 = atof(paramDef);
      k3 = k3/perio3;
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 30,6); equil3 = atof(paramDef);
    }
    equil1 = ANG_360_TO_180(equil1);
    equil2 = ANG_360_TO_180(equil2);
    equil3 = ANG_360_TO_180(equil3);
    
    bZeroCharArray(line_c, LINE_MAX_LEN);
  }

  /* Final insertion from FRCMOD */
    if(atomType1[0] != '\0'){ // Specific Quads
      spTaT1.push_back(atomType1);
      spTaT2.push_back(atomType2);
      spTaT3.push_back(atomType3);
      spTaT4.push_back(atomType4);
    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<line;
    std::cout<<"bAddParams frcmod dihe cnt "<<cnt<<' '<<buff1<<" "<<buff2<<" "<<buff3<<" "<<buff4<<std::endl<<std::flush;
    std::cout<<"bAddParams frcmod dihe "<<Type2Ix[buff1]<<" "<<Type2Ix[buff2]<<" "<<Type2Ix[buff3]<<" "<<Type2Ix[buff4]<<std::endl<<std::flush;
    #endif
      if(cnt == 1){
        dumm.defineBondTorsion(
          (DuMM::AtomClassIndex)(Type2Ix[buff1]),
          (DuMM::AtomClassIndex)(Type2Ix[buff2]),
          (DuMM::AtomClassIndex)(Type2Ix[buff3]),
          (DuMM::AtomClassIndex)(Type2Ix[buff4]),
          1, k1, equil1
          );
      }
      else if(cnt == 2){
        dumm.defineBondTorsion(
          (DuMM::AtomClassIndex)(Type2Ix[buff1]),
          (DuMM::AtomClassIndex)(Type2Ix[buff2]),
          (DuMM::AtomClassIndex)(Type2Ix[buff3]),
          (DuMM::AtomClassIndex)(Type2Ix[buff4]),
          1, k1, equil1,
          2, k2, equil2
          );
      }
      else if(cnt >= 3){
        dumm.defineBondTorsion(
          (DuMM::AtomClassIndex)(Type2Ix[buff1]),
          (DuMM::AtomClassIndex)(Type2Ix[buff2]),
          (DuMM::AtomClassIndex)(Type2Ix[buff3]),
          (DuMM::AtomClassIndex)(Type2Ix[buff4]),
          1, k1, equil1,
          2, k2, equil2,
          3, k3, equil3
          );
      }
    } // if(atomType[0] != '\0')

  /*Read from GAFF*/
  #ifdef DEBUG_PARAMS_LEVEL02
  std::cout<<"bAddParams DIHEDRALS defineBondTorsion gaff"<<std::endl<<std::flush;
  #endif

  prevQuad = "";
  atomType1[0] = '\0'; atomType1[1] = '\0';
  atomType2[0] = '\0'; atomType2[1] = '\0';
  atomType3[0] = '\0'; atomType3[1] = '\0';
  atomType4[0] = '\0'; atomType4[1] = '\0';
  while(fgets(line_c, 200, fpo)){
    if((line_c[0] == '\n') || (line_c[0] == '\r')){
      torsFlag = 0;
      break;
    }
    line = line_c;

    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"bAddParams hop 3 "<<buff1<<" "<<buff2<<" "<<buff3<<" "<<buff4<<std::endl<<std::flush;
    #endif

    if(prevQuad == line.substr(0,11)){
      ++cnt;
      prevQuad = line.substr(0,11);
    }
    else{
      if((prevQuad.substr(0,1) != "X") && (atomType1[0] != '\0')){ // Specific Quads
        spTaT1.push_back(atomType1);
        spTaT2.push_back(atomType2);
        spTaT3.push_back(atomType3);
        spTaT4.push_back(atomType4);
        if(cnt == 1){
          dumm.defineBondTorsion(
            (DuMM::AtomClassIndex)(Type2Ix[buff1]),
            (DuMM::AtomClassIndex)(Type2Ix[buff2]),
            (DuMM::AtomClassIndex)(Type2Ix[buff3]),
            (DuMM::AtomClassIndex)(Type2Ix[buff4]),
            1, k1, equil1
            );
        }
        else if(cnt == 2){
          dumm.defineBondTorsion(
            (DuMM::AtomClassIndex)(Type2Ix[buff1]),
            (DuMM::AtomClassIndex)(Type2Ix[buff2]),
            (DuMM::AtomClassIndex)(Type2Ix[buff3]),
            (DuMM::AtomClassIndex)(Type2Ix[buff4]),
            1, k1, equil1,
            2, k2, equil2
            );
        }
        else if(cnt >= 3){
          dumm.defineBondTorsion(
            (DuMM::AtomClassIndex)(Type2Ix[buff1]),
            (DuMM::AtomClassIndex)(Type2Ix[buff2]),
            (DuMM::AtomClassIndex)(Type2Ix[buff3]),
            (DuMM::AtomClassIndex)(Type2Ix[buff4]),
            1, k1, equil1,
            2, k2, equil2,
            3, k3, equil3
            );
        }
      } // if(line.at(0) != 'X')
      else if(atomType1[0] != '\0'){ // Unspecific Quads (line.at(0) == 'X')
        uspTaT1.push_back(atomType1);
        uspTaT2.push_back(atomType2);
        uspTaT3.push_back(atomType3);
        uspTaT4.push_back(atomType4);
        uspTk1.push_back(k1);
        uspTequil1.push_back(equil1);
      }
      cnt = 1;
      prevQuad = line.substr(0,11);
    }

    // Read atom types
    bSubstr(atomType1, line_c, 0, 2);
    if(atomType1[1] == ' '){atomType1[1] = '\0';}
    else{atomType1[2] = '\0';}
    bSubstr(atomType2, line_c, 3, 2);
    if(atomType2[1] == ' '){atomType2[1] = '\0';}
    else{atomType2[2] = '\0';}
    bSubstr(atomType3, line_c, 6, 2);
    if(atomType3[1] == ' '){atomType3[1] = '\0';}
    else{atomType3[2] = '\0';}
    bSubstr(atomType4, line_c, 9, 2);
    if(atomType4[1] == ' '){atomType4[1] = '\0';}
    else{atomType4[2] = '\0';}
    buff1 = "gaff_"; buff1 += atomType1;
    buff2 = "gaff_"; buff2 += atomType2;
    buff3 = "gaff_"; buff3 += atomType3;
    buff4 = "gaff_"; buff4 += atomType4;

    if(cnt == 1){
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 14,1); perio1 = atof(paramDef);
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 18,5); k1 = atof(paramDef);
      k1 = k1/perio1;
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 30,6); equil1 = atof(paramDef);
    }
    else if(cnt == 2){
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 14,1); perio2 = atof(paramDef);
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 18,5); k2 = atof(paramDef);
      k2 = k2/perio2;
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 30,6); equil2 = atof(paramDef);
    }
    else if(cnt == 3){
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 14,1); perio3 = atof(paramDef);
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 18,5); k3 = atof(paramDef);
      k3 = k3/perio3;
      bZeroCharArray(paramDef, PARAMDEF_MAX_LEN);
      bSubstr(paramDef, line_c, 30,6); equil3 = atof(paramDef);
    }
    equil1 = ANG_360_TO_180(equil1);
    equil2 = ANG_360_TO_180(equil2);
    equil3 = ANG_360_TO_180(equil3);
    
    bZeroCharArray(line_c, LINE_MAX_LEN);
  }

  /* Final insertion from GAFF */
    if(atomType1[0] != '\0'){ // Specific Quads
      spTaT1.push_back(atomType1);
      spTaT2.push_back(atomType2);
      spTaT3.push_back(atomType3);
      spTaT4.push_back(atomType4);
      if(cnt == 1){
        dumm.defineBondTorsion(
          (DuMM::AtomClassIndex)(Type2Ix[buff1]),
          (DuMM::AtomClassIndex)(Type2Ix[buff2]),
          (DuMM::AtomClassIndex)(Type2Ix[buff3]),
          (DuMM::AtomClassIndex)(Type2Ix[buff4]),
          1, k1, equil1
          );
      }
      else if(cnt == 2){
        dumm.defineBondTorsion(
          (DuMM::AtomClassIndex)(Type2Ix[buff1]),
          (DuMM::AtomClassIndex)(Type2Ix[buff2]),
          (DuMM::AtomClassIndex)(Type2Ix[buff3]),
          (DuMM::AtomClassIndex)(Type2Ix[buff4]),
          1, k1, equil1,
          2, k2, equil2
          );
      }
      else if(cnt >= 3){
        dumm.defineBondTorsion(
          (DuMM::AtomClassIndex)(Type2Ix[buff1]),
          (DuMM::AtomClassIndex)(Type2Ix[buff2]),
          (DuMM::AtomClassIndex)(Type2Ix[buff3]),
          (DuMM::AtomClassIndex)(Type2Ix[buff4]),
          1, k1, equil1,
          2, k2, equil2,
          3, k3, equil3
          );
      }
    } // if(atomType[0] != '\0')

    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"bAddParams hop 4 "<<buff1<<" "<<buff2<<" "<<buff3<<" "<<buff4<<std::endl<<std::flush;
    #endif


  /*Add parameters for unspecific quads (from unspecific params)*/
  //vector<bBond> cen; // RESTORE
  //vector<bBond> l; // RESTORE
  //vector<bBond> r; // RESTORE
  //vector<bBond>::iterator it; // RESTORE
  bBond *cen, *l, *r;

  cen = bonds;
  l   = bonds;
  r   = bonds;
  unsigned int ceni=0, li=0, ri=0;
  int CORRECT = 1;  // boolean like
  int INCORRECT = 0;  // l or r need to be swapped
  int lpos, rpos;
  int found_in_specifics = 0;
  int found_in_unspecifics = 0;
  set<string> toBeAdd;
  set<string>::iterator toBeAddIt;
  pair<std::set<string>::iterator, bool> toBeAddRet1;
  pair<std::set<string>::iterator, bool> toBeAddRet2;

  for(ceni=0; ceni<nbnds; ceni++){
    for(li=0; li<nbnds; li++){
      if( ((cen[ceni].i == l[li].i) && (cen[ceni].j != l[li].j)) ||
        ((cen[ceni].i == l[li].j) && (cen[ceni].j != l[li].i)) ){
        //std::cout<<cen[ceni].getString()<<'|'<<l[li].getString()<<std::endl;
        for(ri=0; ri<nbnds; ri++){
          if( ((cen[ceni].j == r[ri].i) && (cen[ceni].i != r[ri].j)) ||
            ((cen[ceni].j == r[ri].j) && (cen[ceni].i != r[ri].i)) ){ // l-cen-r is a quad
            if(cen[ceni].i == l[li].j){lpos=CORRECT;}
            else{lpos=INCORRECT; l[li].swap();}
            if(cen[ceni].j == r[ri].i){rpos=CORRECT;}
            else{rpos=INCORRECT; r[ri].swap();}
            buff1 = bAtomList[l[li].i-1].fftype[5];
            if(bAtomList[l[li].i-1].fftype[6]){buff1 += bAtomList[l[li].i-1].fftype[6];}
            buff2 = bAtomList[cen[ceni].i-1].fftype[5];
            if(bAtomList[cen[ceni].i-1].fftype[6]){buff2 += bAtomList[cen[ceni].i-1].fftype[6];}
            buff3 = bAtomList[cen[ceni].j-1].fftype[5];
            if(bAtomList[cen[ceni].j-1].fftype[6]){buff3 += bAtomList[cen[ceni].j-1].fftype[6];}
            buff4 = bAtomList[r[ri].j-1].fftype[5];
            if(bAtomList[r[ri].j-1].fftype[6]){buff4 += bAtomList[r[ri].j-1].fftype[6];}
            found_in_specifics = 0;
            found_in_unspecifics = 0;
            for(unsigned int s=0; s<spTaT1.size(); s++){
              if( ((spTaT1[s] == buff1) && (spTaT2[s] == buff2) &&
                 (spTaT3[s] == buff3) && (spTaT4[s] == buff4)) ||
                ((spTaT1[s] == buff4) && (spTaT2[s] == buff3) &&
                 (spTaT3[s] == buff2) && (spTaT4[s] == buff1))){
                found_in_specifics = 1;
                break;
              }
            }
            if(!found_in_specifics){
              for(unsigned int us=0; us<uspTaT1.size(); us++){
                if( ((uspTaT2[us] == buff2) && (uspTaT3[us] == buff3)) ||
                  ((uspTaT2[us] == buff3) && (uspTaT3[us] == buff2))){
                  found_in_unspecifics = 1;
                  buff  = buff1; buff += '-';
                  buff += buff2; buff += '-';
                  buff += buff3; buff += '-';
                  buff += buff4;
                  toBeAddRet1 = toBeAdd.insert(buff);
                  buff  = buff4; buff += '-';
                  buff += buff3; buff += '-';
                  buff += buff2; buff += '-';
                  buff += buff1;
                  toBeAddRet2 = toBeAdd.insert(buff);
                  if( ((toBeAddRet1.second == true) && //non-palindromic
                     (toBeAddRet2.second == true) &&
                     ((buff1 != buff4) || (buff2 != buff3))) ||
                    ((toBeAddRet1.second == true) && // palindromic
                     (toBeAddRet2.second == false) &&
                     (buff1 == buff4) && (buff2 == buff3))
                    ){

    #ifdef DEBUG_PARAMS_LEVEL02
    std::cout<<"bAddParams hop 5 "<<bAtomList[l[li].i-1].fftype<<" "<<bAtomList[cen[ceni].i-1].fftype
      <<" "<<bAtomList[cen[ceni].j-1].fftype<<" "<<bAtomList[r[ri].j-1].fftype<<std::endl<<std::flush;
    std::cout<<Type2Ix[bAtomList[l[li].i-1].fftype]<<' '<<Type2Ix[bAtomList[cen[ceni].i-1].fftype]<<' '
      <<Type2Ix[bAtomList[cen[ceni].j-1].fftype]<<' '<<Type2Ix[bAtomList[r[ri].j-1].fftype]<<std::endl;
    #endif

                    dumm.defineBondTorsion(
                    (DuMM::AtomClassIndex)(Type2Ix[bAtomList[l[li].i-1].fftype]),
                    (DuMM::AtomClassIndex)(Type2Ix[bAtomList[cen[ceni].i-1].fftype]),
                    (DuMM::AtomClassIndex)(Type2Ix[bAtomList[cen[ceni].j-1].fftype]),
                    (DuMM::AtomClassIndex)(Type2Ix[bAtomList[r[ri].j-1].fftype]),
                    1,
                    uspTk1[us],
                    uspTequil1[us]
                    );
                  }
                  else{
                    ;// TODO
                  }
                  break;
                }
              }
              if(!found_in_specifics && !found_in_unspecifics){
                ; // TODO
              }
            }

          }
        }
      }
    }
  }

  fclose(fpo);
  fclose(frcmod);    

  #ifdef DEBUG_PARAMS_LEVEL01
  std::cout<<"bAddGaffParams() END"<<std::endl<<std::flush;
  #endif
}  


