#include "bMainResidue.hpp"
using namespace SimTK;
////////////////////////////
////// BMAINRESIDUE ////////
////////////////////////////
/*Any kind of molecule*/
  bMainResidue::bMainResidue(
    DuMMForceFieldSubsystem &dumm,
    unsigned int natms,
    bSpecificAtom *bAtomList,
    std::vector<bBond> bonds,
    TARGET_TYPE **coords,
    TARGET_TYPE **indexMap,
    TARGET_TYPE *PrmToAx_po,
    TARGET_TYPE *MMTkToPrm_po,
    bool first_time,
    std::string ictdF
)
{
    this->natms = natms;
    this->bAtomList = bAtomList;
    this->bonds = bonds;
    this->ictdF = ictdF;
    this->PrmToAx_po = PrmToAx_po;
    this->MMTkToPrm_po = MMTkToPrm_po;

    if (bAtomList == NULL){
      std::cout<<"bMainResidue: NULL bAtomList"<<std::endl<<std::fflush;
      exit(1);
    }
  
    int noDummies = 0;
    std::stringstream sbuff;
    std::stringstream otsbuff;  // other sbuff
    std::string buff[bonds.size()];
    std::string otbuff;        // other buff
    int bi[bonds.size()], bj[bonds.size()];
    unsigned int k;

     /* =================================================
      *    BUILD AN UNLINKED MolModel NESTED IN THIS CLASS 
      * =================================================*/
    bMolAtomList = new MolAtom[natms];
    mol_MolModelCreate ("MainModel", &model);

    for(k=0; k<natms; k++){
      bAtomAssign(&bMolAtomList[k], &bAtomList[k]);
      bMolAtomList[k].id = k;
      bMolAtomList[k].insertion_code = ' ';
      bMolAtomList[k].het = 1;
      bMolAtomList[k].chain_id = 'A';

      mol_MolModelCurrStrucGet (model, &struc);
      mol_StructureAtomAdd (struc, MOL_FALSE, &bMolAtomList[k]);
    }

    mol_StructureChainsBuild(struc, 1);
    struc[0].chain_names = new char*[1];
    struc[0].chain_names[0] = &bMolAtomList[0].chain_id;

     /* ========================================
      *    BUILD THIS CLASS' MOLECULE TOPOLOGY
      * ========================================*/

    dumm.setBondStretchGlobalScaleFactor(0);
    dumm.setBondBendGlobalScaleFactor(0);
    dumm.setBondTorsionGlobalScaleFactor(0);
    dumm.setAmberImproperTorsionGlobalScaleFactor(0);
    dumm.setVdw12ScaleFactor(0);
    dumm.setVdw13ScaleFactor(0);
    dumm.setVdw14ScaleFactor(0);
    dumm.setVdw15ScaleFactor(0);
    dumm.setVdwGlobalScaleFactor(0);
    dumm.setCoulomb12ScaleFactor(0);
    dumm.setCoulomb13ScaleFactor(0);
    dumm.setCoulomb14ScaleFactor(0);
    dumm.setCoulomb15ScaleFactor(0);
    dumm.setCoulombGlobalScaleFactor(0);
    dumm.setGbsaGlobalScaleFactor(0);

    setCompoundName("bMainResidue");

    /*First atom*/
    int found = 0;
    int inter = -1;
    unsigned int fst = 0;  // First bond to be used
    std::vector<int> pushed;
    std::vector<int> crbonds;
    std::vector<int>::iterator pushedIt1;
    std::vector<int>::iterator pushedIt2;

    /* Set up the connectivity */
    for (k=0; k < bonds.size(); k++){
      bi[k] = bonds[k].i - 1;
      bj[k] = bonds[k].j - 1;
    }
    
    /*Set base atom and link second one if dummies not present*/
    if(noDummies == 0){
      for (k=0; k < bonds.size(); k++){
        if((!bonds[k].isRingClosing()) && (bAtomList[bi[k]].freebonds != 1)){
          fst = k;
          break;
        }
      }
    }
    else{    // if dummies
      fst = bonds.size() - 2;
    }
    
    setBaseAtom( *(bAtomList[bi[fst]].bAtomType) );
    setAtomBiotype(bAtomList[bi[fst]].name, "mainRes", bAtomList[bi[fst]].biotype);
    convertInboardBondCenterToOutboard();

    sbuff.str("");
    sbuff<<bAtomList[bi[fst]].name<<"/bond"<<1;
    buff[fst] = sbuff.str();

    bondAtom(*bAtomList[bj[fst]].bAtomType, buff[fst].c_str());
    setAtomBiotype(bAtomList[bj[fst]].name, "mainRes", bAtomList[bj[fst]].biotype);
    bAtomList[bi[fst]].freebonds = 2;    // linked to Ground and bj (++ in loop)
 
    pushed.push_back(bi[fst]);
    pushed.push_back(bj[fst]);
    
    /* Rearrange and connect */
    while(pushed.size() < 2*bonds.size()){
      for(unsigned int m=0; m<bonds.size(); m++){
        if(!bonds[m].isRingClosing()){
          pushedIt1 = find(pushed.begin(), pushed.end(), bj[m]);
          pushedIt2 = find(pushed.begin(), pushed.end(), bi[m]);
          found = 0;
          if((pushedIt1 == pushed.end()) && (pushedIt2 != pushed.end())){
          // bj not found, bi found
            found = 1;
          }
          if((pushedIt1 != pushed.end()) && (pushedIt2 == pushed.end())){
          // bj found, bi not found => swap
            found = 1;
            inter = bi[m]; bi[m] = bj[m]; bj[m] = inter; // swap
          }
          if(found == 1){
            if(m != fst){
              sbuff.str("");
              sbuff<<bAtomList[bi[m]].name<<"/bond"<<bAtomList[bi[m]].freebonds;
              buff[m] = sbuff.str();
              bondAtom(*bAtomList[bj[m]].bAtomType,buff[m].c_str());
              setAtomBiotype(bAtomList[bj[m]].name, "mainRes", bAtomList[bj[m]].biotype);
              if(bi[m] == bi[fst]){
                ++bAtomList[bi[m]].freebonds; // The first atom
              }
              else{
                --bAtomList[bi[m]].freebonds;
              }
            }
            pushed.push_back(bi[m]);
            pushed.push_back(bj[m]);
            break;
          }
        }
      }
      if(found == 0){
        break;
      }
    }
    /*Close the rings*/
      for(unsigned int m=0; m<bonds.size(); m++){
        if(bonds[m].isRingClosing()){
          found = 0;
          pushedIt1 = find(crbonds.begin(), crbonds.end(), m);
          if(pushedIt1 == crbonds.end()){
            found = 1;
          }
          if(found == 1){
              sbuff.str("");
              sbuff<<bAtomList[bi[m]].name
                <<"/bond"<<bAtomList[bi[m]].freebonds;
              buff[m] = sbuff.str();

              otsbuff.str("");
              otsbuff<<bAtomList[bj[m]].name              
                <<"/bond"<<bAtomList[bj[m]].freebonds;
              otbuff = otsbuff.str();
              
              std::cout<<bj[m]+1<<" "<<bi[m]+1<<std::endl
                <<buff[m].c_str()<<' '<<otbuff.c_str()<<std::endl;
              
              addRingClosingBond(
                buff[m].c_str(),
                otbuff.c_str(),
                0.14,
                109*Deg2Rad,
                BondMobility::Rigid
                );
              setAtomBiotype(bAtomList[bi[m]].name, "mainRes", bAtomList[bi[m]].biotype);
              setAtomBiotype(bAtomList[bj[m]].name, "mainRes", bAtomList[bj[m]].biotype);
              if(bi[m] == bi[fst]){
                ++bAtomList[bi[m]].freebonds; // The first atom
              }
              else{
                --bAtomList[bi[m]].freebonds;
              }
          }
          crbonds.push_back(m);
        }
      }

  /*
   Create charged atom types
   Must be called AFTER first mainRes is declared,
   so Biotypes and atom classes will be defined
  */
    std::string abuff;
    for(k=0; k<natms; k++){
      abuff =  "mainRes ";
      abuff += bAtomList[k].biotype;
      dumm.defineChargedAtomType(
        dumm.getNextUnusedChargedAtomTypeIndex(),
        abuff.c_str(),
        dumm.getAtomClassIndex(bAtomList[k].fftype), //Amber
        bAtomList[k].charge
      );
      dumm.setBiotypeChargedAtomType( 
        dumm.getChargedAtomTypeIndex(abuff.c_str()),
        Biotype::get("mainRes", bAtomList[k].biotype).getIndex()
      );

    }

    /* Assign AtomIndex values to atoms in bAtomList[] */
    unsigned int ix = 0; // MINE change type
    unsigned int jx = 0; // MINE change type
    std::string cname, myname;
    for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
      for(ix = 0; ix<natms; ix++){
        cname = getAtomName(aIx);
        myname = bAtomList[ix].name;
        if(cname == myname){
          bAtomList[ix].atomIndex = aIx;
          break;
        }
      }
    }

    /* Assign BondIndex values to bonds in bonds[] */
    int inumber, jnumber;
    Compound::AtomIndex iIx, jIx;
    for (unsigned int r=0 ; r<getNumBonds(); r++){
      iIx = getBondAtomIndex(Compound::BondIndex(r), 0);
      jIx = getBondAtomIndex(Compound::BondIndex(r), 1);

      for(ix = 0; ix<natms; ix++){
        if(bAtomList[ix].atomIndex == iIx){
          inumber = bAtomList[ix].number;
        }
      }
      for(jx = 0; jx<natms; jx++){
        if(bAtomList[jx].atomIndex == jIx){
          jnumber = bAtomList[jx].number;
        }
      }

      for(unsigned int m=0; m<bonds.size(); m++){
        if(((bonds[m].i == inumber) && (bonds[m].j == jnumber)) ||
           ((bonds[m].i == jnumber) && (bonds[m].j == inumber))){
          bonds[m].setBondIndex(Compound::BondIndex(r));
        }
      }
    }


    /* ============= */
    /* Fill indexMap */
    /* ============= */
   for(ix = 0; ix<natms; ix++){
      //indexMap[bAtomList[ix].atomIndex][0] = bAtomList[ix].atomIndex;
      //indexMap[bAtomList[ix].atomIndex][1] = ix;
      indexMap[ix][0] = ix;
      indexMap[ix][1] = bAtomList[ix].atomIndex;
    }

    for(unsigned int i=0; i<natms; i++){
      for(unsigned int k=0; k<natms; k++){
        if(indexMap[k][1] == i){
          PrmToAx_po[i] = TARGET_TYPE(k);
        }
      }
    }  

    for(unsigned int i=0; i<natms; i++){
      for(unsigned int k=0; k<natms; k++){
        if(indexMap[k][2] == i){
          MMTkToPrm_po[i] = TARGET_TYPE(k);
        }
      }
    }

    /* Create atomTargets from passed coords array*/
    std::map<AtomIndex, Vec3> atomTargets;  
    ix = 0;
    if(first_time == true){ // Take coordinates from memory
      for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
       Vec3 v(coords[ (int)(indexMap[ix][2]) ][0],
              coords[ (int)(indexMap[ix][2]) ][1],
              coords[ (int)(indexMap[ix][2]) ][2]);
        atomTargets.insert(pair<AtomIndex, Vec3>
          (bAtomList[ix].atomIndex, v)
        );
        ix++;
      }
    }
    else{ // Take coordinates from
      int ixi, prmtopi_from_SimTK, prmtopi_from_MMTK;
      for(ix = 0; ix < getNumAtoms(); ++ix){
        ixi                 = indexMap[ix][0];
        prmtopi_from_SimTK  = indexMap[ix][1];
        prmtopi_from_MMTK   = indexMap[ix][2];
        Vec3 v(coords[prmtopi_from_MMTK][0]/10, coords[prmtopi_from_MMTK][1]/10, coords[prmtopi_from_MMTK][2]/10);

        atomTargets.insert(pair<AtomIndex, Vec3>
          (bAtomList[ix].atomIndex, v)
        );
      }
    }

    matchDefaultTopLevelTransform(atomTargets);
    matchDefaultConfiguration(atomTargets, Match_Exact, true, 150.0); //Compound::Match_Idealized

    PdbStructure  pdb(*this);

    for(unsigned int i=0; i<natms; i++){
      const PdbAtom& P = pdb.getAtom(String(bAtomList[i].name), PdbResidueId(1), String(" "));
      std::string s(P.getName());
      const Vec3& PC = P.getCoordinates();
      bAtomList[i].x = PC[0]*10;
      bAtomList[i].y = PC[1]*10;
      bAtomList[i].z = PC[2]*10;
    }

    /* Set all the bonds to ... */
    if(ictdF=="IC"){
      for (int r=0 ; r<getNumBonds(); r++){
        setBondMobility(BondMobility::Free, Compound::BondIndex(r));
      }
    }
    else if(ictdF=="TD"){
      for (int r=0 ; r<getNumBonds(); r++){
        setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
      }
    }
    else{
      fprintf(stderr, "Dynamics type unknown\n");
      exit(1);
    }

    /* Set the Rigidity, rigid RIGID */
    /*
    for(unsigned int m=0; m<bonds.size(); m++){
      if(bonds[m].isInRing()){
        setBondMobility(BondMobility::Rigid, bonds[m].getBondIndex());
        std::cout<<"Bond "<<m<<" Rigidized"<<std::endl;
      }
    }
    */

    std::ostringstream sstream;
    sstream<<"pdbs/sb_"<<"ini"<<".pdb";
    std::string ofilename = sstream.str();
    std::cout<<"Writing pdb file: "<<ofilename<<std::endl;
    std::filebuf fb;
    fb.open(ofilename.c_str(), std::ios::out);
    std::ostream os(&fb);
    pdb.write(os);
    fb.close();

  }

  bMainResidue::~bMainResidue(){
  }
  
////////////////////////////
////// END BMAINRESIDUE ////
////////////////////////////

