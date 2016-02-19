#!/bin/bash

# -std=gnu++0x este pentru multithreading

SIMTK_INCL_FLAGS="-I/home/lspirido/Installers/simbody/SimTK-3.0/include -I/home/lspirido/Installers/molmodel/molmodel-3.0.0/src"
MMTK_INCL_FLAGS="-I/export/apps/canopy/1.5.0/Canopy_64bit/User/include/python2.7/ -I/export/apps/canopy/1.5.0/appdata/canopy-1.5.0.2717.rh5-x86_64/include/python2.7/  -I/opt/rocks/include/python2.6"
ARMA_INCL_FLAGS="-I/home/lspirido/Installers/armadillo-6.400.3/include/"
EIGEN_INCL_FLAGS="-I/export/apps/amber/12/AmberTools/src/mtkpp/src/eigen3b2"
MMTK_LIBS_DIRS="/export/apps/canopy/1.5.0/Canopy_64bit/User/lib/python2.7/site-packages/MMTK/linux2/"
SIMTK_RDIRS_FLAGS="-Wl,-rpath,/home/lspirido/Installers/simbody/SimTK-3.0/lib -Wl,-rpath,simbody"
SIMTK_DIR="./simbody/"
#NUMPY_DIR="./NUMPY/"
LIBGC="./"
SIMTK_LIBS_FLAGS="-lSimTKsimbody -lSimTKmath -lSimTKcommon"
MOLMODEL_LIBS_FLAGS="-lSimTKmolmodel"
MMTK_CXX_OPT="-pthread -fno-strict-aliasing -DNDEBUG -g -O2 -DLIBM_HAS_ERFC -D_LONG64_ -DEXTENDED_TYPES -DUSE_NETCDF_H_FROM_SCIENTIFIC=1 -DNUMPY=1 -O3 -ffast-math -fomit-frame-pointer -keep-inline-functions"
SIMTK_DEPS_FLAGS="-llapack -lblas -lpthread -ldl -lrt"
CXX_OBJ="g++ -Wall -fPIC -g -c"

rm libgmolmodel.so*
rm objs_compilation.out
echo 0 > objs_compilation.out

if [[ gmolmodel/bgeneral.cpp -nt gmolmodel/bgeneral.o ]] || [[ gmolmodel/bgeneral.hpp -nt gmolmodel/bgeneral.o ]]; then
  echo; echo "$CXX_OBJ gmolmodel/bgeneral.cpp $SIMTK_INCL_FLAGS"
  $CXX_OBJ gmolmodel/bgeneral.cpp $SIMTK_INCL_FLAGS -o gmolmodel/bgeneral.o
  echo $? >> objs_compilation.out
fi
if [[ gmolmodel/bArgParser.cpp -nt gmolmodel/bArgParser.o ]] || [[ gmolmodel/bArgParser.hpp -nt gmolmodel/bArgParser.o ]]; then
  echo; echo "$CXX_OBJ gmolmodel/bArgParser.cpp $SIMTK_INCL_FLAGS"
  $CXX_OBJ gmolmodel/bArgParser.cpp $SIMTK_INCL_FLAGS -o gmolmodel/bArgParser.o
  echo $? >> objs_compilation.out
fi
if [[ gmolmodel/bMoleculeReader.cpp -nt gmolmodel/bMoleculeReader.o ]] || [[ gmolmodel/bMoleculeReader.hpp -nt gmolmodel/bMoleculeReader.o ]]; then
  echo; echo "$CXX_OBJ gmolmodel/bMoleculeReader.cpp $SIMTK_INCL_FLAGS"
  $CXX_OBJ gmolmodel/bMoleculeReader.cpp $SIMTK_INCL_FLAGS -o gmolmodel/bMoleculeReader.o
  echo $? >> objs_compilation.out
fi
if [[ gmolmodel/bAddParams.cpp -nt gmolmodel/bAddParams.o ]] || [[ gmolmodel/bAddParams.hpp -nt gmolmodel/bAddParams.o ]]; then
  echo; echo "$CXX_OBJ gmolmodel/bAddParams.cpp $SIMTK_INCL_FLAGS"
  $CXX_OBJ gmolmodel/bAddParams.cpp $SIMTK_INCL_FLAGS -o gmolmodel/bAddParams.o
  echo $? >> objs_compilation.out
fi
if [[ gmolmodel/bMainResidue.cpp -nt gmolmodel/bMainResidue.o ]] || [[ gmolmodel/bMainResidue.hpp -nt gmolmodel/bMainResidue.o ]]; then
  echo; echo "$CXX_OBJ gmolmodel/bMainResidue.cpp $SIMTK_INCL_FLAGS"
  $CXX_OBJ gmolmodel/bMainResidue.cpp $SIMTK_INCL_FLAGS -o gmolmodel/bMainResidue.o
  echo $? >> objs_compilation.out
fi
if [[ gmolmodel/bSystem.cpp -nt gmolmodel/bSystem.o ]] || [[ gmolmodel/bSystem.hpp -nt gmolmodel/bSystem.o ]] || [[ gmolmodel/MidVVIntegrator.cpp -nt gmolmodel/bSystem.o ]] || [[ gmolmodel/MidVVIntegrator.hpp -nt gmolmodel/bSystem.o ]]; then
  echo; echo "g++ -std=gnu++0x -g -Wall -fPIC -c gmolmodel/bSystem.cpp $SIMTK_INCL_FLAGS  $MMTK_INCL_FLAGS $EIGEN_INCL_FLAGS -I/export/apps/canopy/1.5.0/Canopy_64bit/User/lib/python2.7/site-packages/numpy/core/include/ -l_compiled_base -llapack -o gmolmodel/bSystem.o"
  g++ -std=gnu++0x -g -Wall -fPIC -c gmolmodel/bSystem.cpp -DNUMPY=1 $SIMTK_INCL_FLAGS $MMTK_INCL_FLAGS $EIGEN_INCL_FLAGS -I/export/apps/canopy/1.5.0/Canopy_64bit/User/lib/python2.7/site-packages/numpy/core/include/ -l_compiled_base -llapack -o gmolmodel/bSystem.o
  echo $? >> objs_compilation.out
fi

awk '{if($1 != 0) exit(1);}' objs_compilation.out
if [ $? -ne 0 ] ; then
        echo "Gmolmodel objects compilation failed"
        exit 1;
fi


echo; echo "g++ -std=gnu++0x -fPIC -shared -Wl,-soname,libgmolmodel.so.1 -o libgmolmodel.so.1.0  gmolmodel/bMoleculeReader.o  gmolmodel/bAddParams.o  gmolmodel/bArgParser.o  gmolmodel/bMainResidue.o   gmolmodel/bSystem.o   gmolmodel/bgeneral.o  -I./ -L$LIBGC  -g  -Wl,-R./ -Wl,-Bstatic $SIMTK_INCL_FLAGS  -Wl,-Bdynamic  $SIMTK_RDIRS_FLAGS -L$SIMTK_DIR -L$LIBGC $SIMTK_LIBS_FLAGS $MOLMODEL_LIBS_FLAGS $SIMTK_DEPS_FLAGS"
g++ -std=gnu++0x -g -fPIC -shared -Wl,-soname,libgmolmodel.so.1.0 -o gmolmodel/libgmolmodel.so.1.0  gmolmodel/bMoleculeReader.o  gmolmodel/bAddParams.o  gmolmodel/bArgParser.o  gmolmodel/bMainResidue.o  gmolmodel/bSystem.o   gmolmodel/bgeneral.o -I./ -LLIBGC  -g  -Wl,-R./ -Wl,-Bstatic $SIMTK_INCL_FLAGS  -Wl,-Bdynamic  $SIMTK_RDIRS_FLAGS -L$SIMTK_DIR -L$LIBGC $SIMTK_LIBS_FLAGS $MOLMODEL_LIBS_FLAGS $SIMTK_DEPS_FLAGS
echo $? | awk '{if($1==0){ print "Gmolmodel libraries compilation succesful!";}else print "Gmolmodel libraries compilation failed";}'

cd gmolmodel/
ln -sf libgmolmodel.so.1.0 libgmolmodel.so.1
ln -sf libgmolmodel.so.1 libgmolmodel.so
cp libgmolmodel.so.* ../
cd ../



