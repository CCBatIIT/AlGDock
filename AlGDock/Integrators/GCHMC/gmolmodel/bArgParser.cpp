#include "bArgParser.hpp"
////////////////////////////
////// ARGPARSER ///////////
////////////////////////////

bArgParser::bArgParser(int argc, const char **argv){ // constructor

  MAX_NO_OPT = argc/2;
  std::cout<<"MAX_NO_OPT "<<MAX_NO_OPT<<std::endl;
  option = new int[MAX_NO_OPT];
  for(int i=0; i<MAX_NO_OPT; i++){
    option[i] = 0;
  }

  // * Set the argument strings * //
  for(int a=0; a<argc; a++){
    if((option[AR_LIGDIR] == 1) && (argv[a][0] != '-')){
      mol2F = argv[a];
      mol2F += "/ligand.mol2";
      rbF = argv[a];
      rbF += "/ligand.rb";
      frcmodF = argv[a];
      frcmodF += "/ligand.frcmod";
    }
    else if((option[AR_GAFF]   == 1) && (argv[a][0] != '-')){
      gaffF = argv[a];
    }
    else if((option[AR_ICTD] == 1) && (argv[a][0] != '-')){
      ictdF = argv[a];
    }
  
    // * Set the flags * //
    if(strcmp(argv[a], "-ligdir") == 0){
      for(int i=0; i<MAX_NO_OPT; i++){
        option[i] = 0;
      }
      option[AR_LIGDIR] = 1;
    }
    if(strcmp(argv[a], "-gaff") == 0){
      for(int i=0; i<MAX_NO_OPT; i++){
        option[i] = 0;
      }
      option[AR_GAFF] = 1;  
    }
    if(strcmp(argv[a], "-ictd") == 0){
      for(int i=0; i<MAX_NO_OPT; i++){
        option[i] = 0;
      }
      option[AR_ICTD] = 1;  
    }
  }
  printf("bArgParser END\n"); fflush(stdout);

} // constructor
 
void bArgParser::Print(void){
  std::cout<<"bArgParser arguments are:"<<std::endl;
  std::cout<<"mol2 file: "<<mol2F<<std::endl;
  std::cout<<"rb file:   "<<rbF<<std::endl;
  std::cout<<"gaff file: "<<gaffF<<std::endl;
  std::cout<<"frcmod file: "<<frcmodF<<std::endl;
  std::cout<<"dynamics type : "<<ictdF<<std::endl;
  //std::cout<<"memid : "<<args[AR_MEMID]<<std::endl;
}


bArgParser::~bArgParser(){} // destructor

////////////////////////////
////// END ARGPARSER ///////
////////////////////////////

