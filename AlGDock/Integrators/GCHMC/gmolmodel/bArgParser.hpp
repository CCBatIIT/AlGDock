#ifndef BARGPARSER_H_
#define BARGPARSER_H
#include <stdio.h>
#include <iostream>
#include <string>
#include <string.h>

//==============================================================================
//                           CLASS ArgParser
//==============================================================================
/**
 * Argument Parser Class. It is used by the GCHMCIntegrator class in simmain.
 * It holds the pathnames of the files needed to construct the compound: a mol2
 * file, a rigid body specification file, the Amber gaff file and the type of 
 * dynamics that it supposed to be done.
 **/
class bArgParser{
 public:
  int MAX_NO_OPT;
  int *option;
  std::string mol2F, rbF, gaffF, frcmodF, ictdF; 
  enum{
    AR_LIGDIR,
    AR_GAFF,
    AR_ICTD
  };

  bArgParser(int argc, const char **argv); // constructor
  ~bArgParser(); // destructor
  void Print(void);
};

#endif //BARGPARSER_H_
