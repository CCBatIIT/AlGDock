
INCL_FLAGS[$CURR]="-I./"
INCL_FLAGS[$GMOL]="-Igmolmodel/"
INCL_FLAGS[$SIMB]="-Isimbody/"
INCL_FLAGS[$SIMTK]="-I/home/lspirido/Installers/simbody/SimTK-3.0/include -I/home/lspirido/Installers/molmodel/molmodel-3.0.0/src"
INCL_FLAGS[$MMTK]="-I/export/apps/canopy/1.5.0/Canopy_64bit/User/include/python2.7/ -I/export/apps/canopy/1.5.0/appdata/canopy-1.5.0.2717.rh5-x86_64/include/python2.7/  -I/opt/rocks/include/python2.6"
INCL_FLAGS[$NUMPY]="-I/export/apps/canopy/1.5.0/Canopy_64bit/User/lib/python2.7/site-packages/numpy/core/include/ -I/opt/rocks/include/python2.6/"
INCL_FLAGS[$EIGEN]="-I/export/apps/amber/12/AmberTools/src/mtkpp/src/eigen3b2"

RDIRS_FLAGS[$CURR]="-Wl,-rpath,./"
RDIRS_FLAGS[$GMOL]="-Wl,-rpath,gmolmodel"
RDIRS_FLAGS[$SIMB]="-Wl,-rpath,simbody"
RDIRS_FLAGS[$SIMTK]="-Wl,-rpath,/home/lspirido/Installers/simbody/SimTK-3.0/lib"
RDIRS_FLAGS[$CANOPY]="-Wl,-rpath,/export/apps/canopy/1.5.0/Canopy_64bit/User/lib/python2.7/site-packages/MMTK/linux2"
RDIRS_FLAGS[$MMTK]="-Wl,-rpath,MMTK/linux2/"
RDIRS_FLAGS[$NUMPY]="-Wl,-rpath,NUMPY"
RDIRS_FLAGS[$GCHMC]="-Wl,-rpath,GCHMC/"

LDIRS_FLAGS[$CURR]="-L./"
LDIRS_FLAGS[$GMOL]="-Lgmolmodel/"
LDIRS_FLAGS[$SIMB]="-Lsimbody/"
LDIRS_FLAGS[$MMTK]="-L/export/apps/canopy/1.5.0/Canopy_64bit/User/lib/python2.7/site-packages/MMTK/linux2/"
LDIRS_FLAGS[$NUMPY]="-LNUMPY/"
LDIRS_FLAGS[$GCHMC]="-LGCHMC/"

LIBS_FLAGS[$GMOL]="-lgmolmodel"
LIBS_FLAGS[$SIMTK]="-lSimTKsimbody -lSimTKmath -lSimTKcommon"
LIBS_FLAGS[$NUMPY]="-l_compiled_base"
LIBS_FLAGS[$MOLMO]="-lSimTKmolmodel"
LIBS_FLAGS[$BOOST]="-lboost_python"
##############################################################

CXX="g++"
MMTK_COMP_OPT="-pthread -fno-strict-aliasing -DNDEBUG -g -O2 -DLIBM_HAS_ERFC -D_LONG64_ -DEXTENDED_TYPES -DUSE_NETCDF_H_FROM_SCIENTIFIC=1 -DNUMPY=1 -O3 -ffast-math -fomit-frame-pointer -keep-inline-functions"
MMTK_CXX_OPT="-pthread -fno-strict-aliasing -DNDEBUG -g -O2 -DLIBM_HAS_ERFC -D_LONG64_ -DEXTENDED_TYPES -DUSE_NETCDF_H_FROM_SCIENTIFIC=1 -DNUMPY=1 -O3 -ffast-math -fomit-frame-pointer -keep-inline-functions"
LIBGC="./"
SIMTK_DEPS_FLAGS="-llapack -lblas -lpthread -ldl -lrt"
REC="-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"
EXTRA_OPT="-Xlinker -export-dynamic -std=gnu++0x"
DEBUG_OPT="-Wall"
##############################################################

/home/lspirido/Installers/boost_1_55_0/tools/build/v2/engine/bin.linuxx86_64/bjam cxxflags="$DEBUG_OPT $EXTRA_OPT $MMTK_COMP_OPT ${RDIRS_FLAGS[*]} ${LDIRS_FLAGS[*]} ${LIBS_FLAGS[*]}  ${DEPS_FLAGS[*]}"; >log 2>&1;   echo "================="; tail -5 log; echo "################"; egrep "error|Error|failed|fault" log;



