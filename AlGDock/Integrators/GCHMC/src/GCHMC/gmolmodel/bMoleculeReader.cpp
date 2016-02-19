#include "bMoleculeReader.hpp"

#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif

#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif

//using namespace std;
using namespace SimTK;

/****
 *  Trivalent Atom Class with tetrahedral geometry.
 *  Bond centers are named "bond1", "bond2", and "bond3"
 ****/
TrivalentAtomTetra::TrivalentAtomTetra(
        const Compound::AtomName& atomName,   ///< name for new atom
        const Element& element   ///< chemical element for new atom
) : Compound::SingleAtom(atomName, element)
{
  static const Angle TetrahedralAngle = 109.47 * Deg2Rad;
  // BondCenter1 dihedral will be relative to BondCenter2
  addFirstBondCenter( "bond1", atomName );
  // bond centers 2 and 3 dihedrals relative to bond center 1
  //addSecondBondCenter( "bond2", atomName,  TetrahedralAngle);
  addLeftHandedBondCenter( "bond2", atomName, TetrahedralAngle, TetrahedralAngle );
  addRightHandedBondCenter( "bond3", atomName, TetrahedralAngle, TetrahedralAngle );
  // Choice of inboard bond may differ from bond priority - user may change this
  setInboardBondCenter("bond1");
  setCompoundName("TrivalentAtomTetra"); // overridden
}


/****
 * bPDBReader
 ****/
bPDBReader::bPDBReader()
{
    ;
}
bPDBReader::~bPDBReader()
{
    ;
}


/****
 * bSpecificAtom
 ****/
bSpecificAtom::bSpecificAtom(){
    nbonds = 0;
    freebonds = 0;
    number = 0;
    bZeroStr(name);
    bZeroStr(mol2name);
    bZeroStr(fftype);
    bZeroStr(biotype);
    x = -999;
    y = -999;
    z = -999;
  }

bSpecificAtom::~bSpecificAtom(){;}

/********************
 *     FUNCTIONS
 * ******************/
int bAtomAssign(MolAtom *dest, const bSpecificAtom *src)
{
  dest->name = new char[5];
  strncpy(dest->name, src->name, 4);

  dest->num = src->number;

  switch(src->elem){
    case 'C': dest->type = MOL_ATOM_ELEMENT_CARBON; break;
    case 'H': dest->type = MOL_ATOM_ELEMENT_HYDROGEN; break;
    case 'N': dest->type = MOL_ATOM_ELEMENT_NITROGEN; break;
    case 'O': dest->type = MOL_ATOM_ELEMENT_OXYGEN; break;
    case 'S': dest->type = MOL_ATOM_ELEMENT_SULFUR; break;
    case 'P': dest->type = MOL_ATOM_ELEMENT_PHOSPHORUS; break;
    default: dest->type = MOL_ATOM_ELEMENT_UNKNOWN; break;
  }

  dest->pos[0] = src->x;
  dest->pos[1] = src->y;
  dest->pos[2] = src->z;

  return 0;
}

/********************************/

/****
 * intpair
 ****/
intpair::intpair(){
    i = 0; j = 0;
  }
intpair::intpair(int inI, int inJ){
    this->i = inI;
    this->j = inJ;
}
intpair::~intpair(){}

bool intpair::operator==(const intpair *other){
  return (
    ((this->i == other->i) && (this->j == other->j)) ||
    ((this->i == other->j) && (this->j == other->i))
  );
}

bool intpair::operator!=(const intpair *other){
  return (
    ((this->i != other->i) || (this->j != other->j)) &&
    ((this->i != other->j) || (this->j != other->i))
  );
}

bool intpair::isTheSameAs(const intpair *other){
  return (
    ((this->i == other->i) && (this->j == other->j)) ||
    ((this->i == other->j) && (this->j == other->i))
  );
}

void intpair::swap(void){
  int inter;
  inter = this->i;
  this->i = this->j;
  this->j = inter;
}

void intpair::dump(void){
  std::cout<<i<<' '<<j<<std::endl;
}

std::string intpair::getString(void){
  std::stringstream ret;
  ret<<i<<' '<<j;
  return ret.str();
}

/********************************/


/****
 * bBond
 ****/
bBond::bBond(void) : intpair(){
  inring = 0;
  rigid = 0;
  ring_closing = 0; // later
  ring_no = 0; // later
}

bBond::bBond(int a, int b) : intpair(a, b){
  inring = 0;
  rigid = 0;
  ring_closing = 0;
  ring_no = 0; // later
}

bBond::~bBond(void){;}

bool bBond::isInRing(void){
  if(this->inring == 1)
    return true;
  else
    return false;
}

bool bBond::isRingClosing(void){
  if(this->ring_closing == 1)
    return true;
  else
    return false;
}

bool bBond::isRigid(void){
  if(this->rigid == 1)
    return true;
  else
    return false;
}

int bBond::ringNo(void){
  return this->ring_no;
}

void bBond::setInRing(void){
  this->inring = 1;
}

void bBond::setAsRingClosing(void){
  this->ring_closing = 1;
}

void bBond::setAsRigid(void){
  this->rigid = 1;
}

void bBond::setRingNo(int rn){
  this->ring_no = rn;
}

Compound::BondIndex bBond::getBondIndex(void){
  return bondIndex;
}
  
void bBond::setBondIndex(Compound::BondIndex otherIx){
  this->bondIndex = otherIx;
}



/********************************/

/****
 * intriad
 ****/
intriad::intriad(){
  i=0; j=0; k=0;
}
intriad::intriad(int inI, int inJ, int inK){
  this->i = inI;
  this->j = inJ;
  this->k = inK;
}
intriad::~intriad(){}
//Sorry for crowding
bool intriad::operator==(const intriad *other){
  return (
    ((this->i==other->i)&&(this->j==other->j)&&(this->k==other->k))||
    ((this->i==other->j)&&(this->j==other->i)&&(this->k==other->k))||
    ((this->i==other->k)&&(this->j==other->i)&&(this->k==other->j))||
    ((this->i==other->k)&&(this->j==other->j)&&(this->k==other->i))||
    ((this->i==other->i)&&(this->j==other->k)&&(this->k==other->j))||
    ((this->i==other->j)&&(this->j==other->k)&&(this->k==other->i))
  );
}
bool intriad::isTheSameAs(const intriad *other){
  return (
    ((this->i==other->i)&&(this->j==other->j)&&(this->k==other->k))||
    ((this->i==other->k)&&(this->j==other->j)&&(this->k==other->i))
  );
}
void intriad::dump(void){
  std::cout<<i<<' '<<j<<' '<<k<<std::endl;
}
std::string intriad::getString(void){
  std::stringstream ret;
  ret<<i<<' '<<j<<' '<<k;
  return ret.str();
}

/********************************/

/****
 * bMoleculeReader
 ****/
bMoleculeReader::bMoleculeReader(DuMMForceFieldSubsystem& dumm,
        const char *filename,
        const char *filetype,
        const char *rbfilename){
  #ifdef DEBUG_LEVEL02
  std::cout<<"bMoleculeReader::bMoleculeReader() BEGIN"<<std::endl<<std::flush;
  #endif

  FILE *fpo;
  fpo = fopen(filename, "r");
  /*rewind(fpo);*/ /*Doesn't work for Windows files*/
  MAX_LINE_LENGTH = 10000;
  char *buff = new char[80];
  //bZeroStr(buff);
  bZeroCharArray(buff, 80);
  natms = 0;
  //int nbonds = 0;
  bBond buffij;
  unsigned int i = 0;
  unsigned int j = 0;
  char line_c[MAX_LINE_LENGTH];
  std::string line;
  char elem = 'x';
  unsigned int lno = 0;
  int noDummies = 0;

  /*+++++++++ SDF type ++++++++++*/
  if(strstr(filetype, "sdf")){
    ; // TODO
  }

  /*+++++++++ MOL2 type ++++++++++*/
  else if(strstr(filetype, "mol2")){
    bZeroStr(buff);
    natms = 0;
    i = 0; j = 0;
    char elem = 'x';
    unsigned int lno = 0;
 
    /*Read number of atoms*/
    while(fgets(line_c, MAX_LINE_LENGTH, fpo)){
      #ifdef DEBUG_LEVEL02
      std::cout<<"line: "<<line_c<<std::flush;
      #endif
      ++lno;
      if(lno == 3){
        bSubstr(buff, line_c, 0, 5);
        natms = atoi(buff);
        #ifdef DEBUG_LEVEL02
        std::cout<<"natms "<<natms<<" noDummies "<<noDummies<<std::endl;
        #endif
        bZeroStr(buff);

        //bSubstr(buff, line_c, 6, 6); // EU
        //nbnds = atoi(buff); // EU
        //#ifdef DEBUG_LEVEL02 // EU
        //std::cout<<"nbnds: "<<nbnds<<std::endl<<std::flush; // EU
        //#endif // EU
        //bZeroStr(buff); // EU

        break;
      }
    }
    bAtomList = new bSpecificAtom[natms + noDummies]; /*Alloc - 1 for dummy*/

    #ifdef DEBUG_LEVEL02
    if(bAtomList == NULL){
      std::cout<<"bMoleculeReader: bAtomList is NULL"<<std::endl<<std::flush;
    }
    else{
      std::cout<<"bMoleculeReader: bAtomList allocated"<<std::endl<<std::flush;
    }
    #endif

    /*Jump to the atom section*/
    while(fgets(line_c, MAX_LINE_LENGTH, fpo)){
      if(line_c[0] == '@'){
        break;
      }
    }

    /*Read position, element, fftype and charge*/
    lno = 0;
    while(fgets(line_c, MAX_LINE_LENGTH, fpo) && (lno < natms)){
      ++lno;
      if(line_c != NULL){
        std::cout<<"bMoleculeReader::bMoleculeReader() line_c is not NULL and has "<<strlen(line_c)<<" chars"<<std::endl<<std::flush;
      }
      //line = line_c; // RESTORE
      bZeroStr(bAtomList[lno-1].name);
      bZeroStr(bAtomList[lno-1].fftype);
      //elem = line.at(8); // RESTORE
      elem = line_c[8]; // EU
      sprintf(buff, "%c%d", elem, lno); /*This is not the name from mol2*/
      strncpy(bAtomList[lno-1].name, buff, 4);
      bAtomList[lno-1].number = lno;
      bAtomList[lno-1].elem = elem;
      bZeroStr(buff);
      bSubstr(buff, line_c, 47,2);
      sprintf(bAtomList[lno-1].fftype, "gaff_%s",buff);
      if(bAtomList[lno-1].fftype[6] == ' '){bAtomList[lno-1].fftype[6] = '\0';}
      bZeroStr(buff);
      bSubstr(buff, line_c, 67,9);
      bAtomList[lno-1].charge = atof(buff);

      bZeroStr(buff);
      bSubstr(buff, line_c, 8,4);
      strncpy(bAtomList[lno-1].mol2name, buff, 4);

      bZeroStr(buff);
      bSubstr(buff, line_c, 17,9);
      bAtomList[lno-1].x = atof(buff);
      bZeroStr(buff);
      bSubstr(buff, line_c, 27,9);
      bAtomList[lno-1].y = atof(buff);
      bZeroStr(buff);
      bSubstr(buff, line_c, 37,9);
      bAtomList[lno-1].z = atof(buff);
    }
    #ifdef DEBUG_LEVEL02
    std::cout<<"bMoleculeReader::bMoleculeReader 2"<<std::endl<<std::flush;
    #endif


    /*READ BONDS*/
    fgets(line_c, MAX_LINE_LENGTH, fpo); /*Skip the @<TRIPOS>BOND line*/
    lno = 0;
    do{
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() bond line= "<<line_c<<std::flush;
      #endif
      ++lno;
      //line = line_c; // RESTORE
      sscanf(line_c, "%d%d%d", &i, &buffij.i, &buffij.j);
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() bond line sscanf "<<std::endl<<std::flush;
      std::cout<<"bMoleculeReader::bMoleculeReader() i, buffij.i buffij.j "
        <<i<<" "<<buffij.i<<" "<<buffij.j<<std::endl<<std::flush;
      #endif
      bonds.push_back(buffij);
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() bonds push "<<std::endl<<std::flush;
      #endif
    }
    while(fgets(line_c, MAX_LINE_LENGTH, fpo) && (line_c[0] != '@'));

    
    /*Assign atoms nbonds and freebonds*/
    for(i=0; i<natms;i++){
      bAtomList[i].nbonds = 0;
      for(j=0; j<bonds.size(); j++){
        if((bAtomList[i].number == bonds[j].i) ||\
           (bAtomList[i].number == bonds[j].j)){
          ++bAtomList[i].nbonds;
          ++bAtomList[i].freebonds;
        }
      }
    }

    #ifdef DEBUG_LEVEL02
    std::cout<<"bMoleculeReader::bMoleculeReader() nbonds & freebonds assigned "<<line_c<<std::endl<<std::fflush;
    #endif


    /*TODO Develop TrivalentAtomTetra for adding custom angles*/
    // Every *valentAtom is derived from SingleAtom in turn derived from Compound with one atom (AtomIndex 0)
    for(i=0; i<natms+noDummies; i++){
      if(bAtomList[i].nbonds == 1){
        if(toupper(bAtomList[i].elem) == 'H'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(1, "Hydrogen", "H", 1.008));
            UnivalentAtom(bAtomList[i].name, Element::Hydrogen());
        }
        else if(toupper(bAtomList[i].name[0]) == 'C'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(17, "Chlorine", "Cl", 35.45));
            UnivalentAtom(bAtomList[i].name, Element::Chlorine());
        }
        else if(toupper(bAtomList[i].name[0]) == 'O'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(8, "Oxygen", "O", 16.00));
            UnivalentAtom(bAtomList[i].name, Element::Oxygen());
        }
        else if(toupper(bAtomList[i].name[0]) == 'F'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(9, "Fluorine", "F", 19.00));
            UnivalentAtom(bAtomList[i].name, Element::Fluorine());
        }
        else if(toupper(bAtomList[i].name[0]) == 'B'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(35, "Bromine", "Br", 79.90));
            UnivalentAtom(bAtomList[i].name, Element::Bromine());
        }
        else if(toupper(bAtomList[i].name[0]) == 'I'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(53, "Iodine", "I", 126.9));
            UnivalentAtom(bAtomList[i].name, Element::Iodine());
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.1112);
      }
      else if (bAtomList[i].nbonds == 2){
        if(toupper(bAtomList[i].elem) == 'H'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name, Element(1, "Hydrogen", "H", 1.008));
            BivalentAtom(bAtomList[i].name, Element::Hydrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'C'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", 12.01));
            BivalentAtom(bAtomList[i].name,  Element::Carbon());
        }
        else if(toupper(bAtomList[i].elem) == 'O'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", 16.00),
            BivalentAtom(bAtomList[i].name,  Element::Oxygen(),
            109.47*Deg2Rad);
        }
        else if(toupper(bAtomList[i].elem) == 'N'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", 14.01));
            BivalentAtom(bAtomList[i].name,  Element::Nitrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'S'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", 32.06),
            BivalentAtom(bAtomList[i].name,  Element::Sulfur(),
            109.47*Deg2Rad);
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
      }
      else if (bAtomList[i].nbonds == 3){
        if(toupper(bAtomList[i].elem) == 'C'){
          bAtomList[i].bAtomType = new
            //TrivalentAtom(bAtomList[i].name, Element(6, "Carbon", "C", 12.01),
            TrivalentAtom(bAtomList[i].name, Element::Carbon(),
              120*Deg2Rad, 120*Deg2Rad
            );
        }
        else if(toupper(bAtomList[i].elem) == 'O'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(8, "Oxygen", "O", 16.00));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Oxygen());
        }
        else if(toupper(bAtomList[i].elem) == 'N'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(7, "Nitrogen", "N", 14.01));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Nitrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'S'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(16, "Sulfur", "S", 32.06));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Sulfur());
        }
        else if(toupper(bAtomList[i].elem) == 'P'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(15, "Phosphorus", "P", 30.97));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Phosphorus());
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
      }
      else if (bAtomList[i].nbonds == 4){
        if(toupper(bAtomList[i].elem) == 'C'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", 12.01));
            QuadrivalentAtom(bAtomList[i].name,  Element::Carbon());
        }
        else if(toupper(bAtomList[i].elem) == 'O'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", 16.00));
            QuadrivalentAtom(bAtomList[i].name,  Element::Oxygen());
        }
        else if(toupper(bAtomList[i].elem) == 'N'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", 14.01));
            QuadrivalentAtom(bAtomList[i].name,  Element::Nitrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'S'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", 32.06));
            QuadrivalentAtom(bAtomList[i].name,  Element::Sulfur());
        }
        else if(toupper(bAtomList[i].elem) == 'P'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(15, "Phosphorus", "P", 30.97));
            QuadrivalentAtom(bAtomList[i].name,  Element::Phosphorus());
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
      }

      bZeroStr(bAtomList[i].biotype);
      sprintf(bAtomList[i].biotype, "%s_%s", \
        bAtomList[i].name, bAtomList[i].fftype);
    }

    if (bAtomList == NULL){
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader constructor exit: NULL bAtomList\n";fflush(stdout);
      #endif
    }
    else{
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader constructor: bAtomList not NULL\n";fflush(stdout);
      #endif
    }
    
  }

  /* Just checking *////////
  #ifdef DEBUG_LEVEL02
  std::cout<<"Just checking\n";
  for(i=0; i<natms;i++){
    printf(" -- name(%s) mol2name(%s) number(%d) elem(%c) fftype(%s) biotype(%s) charge(%f)\n", 
      bAtomList[i].name, bAtomList[i].mol2name, bAtomList[i].number, 
      bAtomList[i].elem, bAtomList[i].fftype,
      bAtomList[i].biotype, bAtomList[i].charge);
    fflush(stdout);
  }
  for(i=0; i<bonds.size(); i++){
    std::cout<<"bond: "<<bonds[i].i<<" "<<bonds[i].j<<std::endl;
    fflush(stdout);
  }
  #endif
  ///////////////////////////


  fclose(fpo);
  #ifdef DEBUG_LEVEL02
  std::cout<<"File closed succesfully\n";fflush(stdout);
  #endif
  //////////////////////////////////////////


  /*Now read rigid bodies specifications*/
  FILE *rfpo;
  rfpo = fopen(rbfilename, "r");
  if(rfpo == NULL){
    printf("Usage:\n<program> -mol2 <mol2_file> -rb <rb_file> -gaff  gaff.dat -frcmod <frcmod_file>\n");
    printf("rb_file not provided. Exiting...\n");
    exit(1);
  }

  #ifdef DEBUG_LEVEL02
  std::cout<<"\nREAD RIGID BODIES BEGIN\n"; fflush(stdout); 
  #endif

  std::string sbuff;
  std::vector<int> ring;
  char curr_char;
  int bond_found = 0;
  unsigned int boi;
  int par_open = 0;
  int ring_no = -1;
  int ring_closing = 0;

  while(fgets(line_c, MAX_LINE_LENGTH, rfpo)){
    line = line_c; //RESTORE
    if((sbuff = line.substr(0,5)) == "rings"){ // RESTORE
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() sbuff"<<sbuff<<std::endl;
      #endif
      sbuff = "";
      //for(i=0; i<line.size(); i++){ // RESTORE
      for(i=0; i<strlen(line_c); i++){ // EU
        //curr_char = line.at(i); // RESTORE
        curr_char = line_c[i]; // RESTORE
        if(curr_char == ']'){
          ++ring_no;
          ring.push_back(atoi(sbuff.c_str()));
          bond_found = 0;
          for(boi=0; boi<bonds.size(); boi++){
            for(unsigned int ri=0; ri<ring.size(); ri++){
              for(unsigned int rj=0; rj<ring.size(), rj<ri; rj++){
                if( ((ring[ri] == bonds[boi].i) && (ring[rj] == bonds[boi].j)) ||
                  ((ring[ri] == bonds[boi].j) && (ring[rj] == bonds[boi].i))){
                  bonds[boi].setInRing();
                  bonds[boi].setRingNo(ring_no);
                  if(ring_closing == 0){
                    bonds[boi].setAsRingClosing();
                  }
                  #ifdef DEBUG_LEVEL02
                  std::cout<<"RING BOND: "<<bonds[boi].i
                    <<' '<<bonds[boi].j<<' '
                    <<bonds[boi].ringNo()<<' '
                    <<bonds[boi].isRingClosing()<<std::endl;
                  #endif
                  ++ring_closing;
                  bond_found = 1;
                  break;
                }
              }
            }
          }
          ring.clear();
          ring_closing = 0;
        }
        else if(curr_char == ','){
          ring.push_back(atoi(sbuff.c_str()));
          sbuff = "";
        }
        else if((curr_char >= '0') && (curr_char <='9')){
          sbuff += curr_char;
        }
      }
    }

    else if((sbuff = line.substr(0,17)) == "non_ring_pi_bonds"){ // RESTORE
      #ifdef DEBUG_LEVEL02
      std::cout<<std::endl<<sbuff<<std::endl<<std::flush;
      #endif
      sbuff = "";
      for(i=0; i<line.size(); i++){
        //curr_char = line.at(i); // RESTORE
        curr_char = line_c[i]; // RESTORE
        if(curr_char == ')'){
          par_open = 0;
          ring.push_back(atoi(sbuff.c_str()));
          // Set the bond as rigid
          for(boi = 0; boi<bonds.size(); boi++){
            if( (bonds[boi].i == ring[0]) && (bonds[boi].j == ring[1]) ||
              (bonds[boi].i == ring[1]) && (bonds[boi].j == ring[0])){
              bonds[boi].setAsRigid();
              #ifdef DEBUG_LEVEL02
              std::cout<<"RIGID BOND: "<<bonds[boi].i
                <<' '<<bonds[boi].j<<std::endl;
              fflush(stdout);
              #endif
              ring.clear();
              sbuff = "";
              break;
            }
          }
        }
        else if((curr_char == ',') && (par_open == 1)){
          ring.push_back(atoi(sbuff.c_str()));
          sbuff = "";
        }
        else if((curr_char >= '0') && (curr_char <='9')){
          sbuff += curr_char;
        }
        else if(curr_char == '('){
          par_open = 1;
        }
      }
    }
  }

  #ifdef DEBUG_LEVEL01
  std::cout<<"\nREAD RIGID BODIES END\n\n"; fflush(stdout); 
  #endif

  fclose(rfpo);
}

bMoleculeReader::~bMoleculeReader(){;}




