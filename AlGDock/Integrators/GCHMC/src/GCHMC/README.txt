to build:
  export BOOST_ROOT=/home/lspirido/Installers/boost_1_55_0
  export BOOST_BUILD_PATH=/home/lspirido/Installers/boost_1_55_0
  export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH 

./compile_gmolmodel.bash  >log 2>&1; egrep "error|Error|failed|fault|===" log | head
rm GCHMC.so
./bjam.bash


python simulator09_GCHMC.py --mol_name=2butanol >out 2>&1

