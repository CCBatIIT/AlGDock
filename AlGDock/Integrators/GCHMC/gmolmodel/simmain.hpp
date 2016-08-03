#ifndef SIMMAIN_HPP
#define SIMMAIN_HPP
//#include <boost/python/module.hpp>
//#include <boost/python/def.hpp>
#include <boost/python.hpp>
using namespace boost::python;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <sys/types.h>
#include <sys/ipc.h>

#include "Simbody.h"
#include "Molmodel.h"
#include "bMoleculeReader.hpp"
#include "bAddParams.hpp"
#include "bSystem.hpp"
#include "server.hpp"

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
//#include "MMTK/trajectory.h"

#include <numpy/arrayobject.h>


#endif // SIMMAIN_HPP


