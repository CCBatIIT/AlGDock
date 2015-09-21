#!/bin/bash

export WORK_DIR=`pwd`

# Make sure the latest AlGDock is installed
cd ..
python setup.py install

# Clean the old distribution
cd $WORK_DIR
rm -rf dist build

# Run pyinstaller
pyinstaller BindingPMF.spec

# Compress the final results
cd dist
tar czf algdock.tar.gz AlGDock/

# Copy the compressed file to the OSG grid and to CCB
# scp algdock.tar.gz dminh@ccb.tbc.iit.edu:/share/apps/algdock/
scp algdock.tar.gz daveminh@login.osgconnect.net:~/public/
cp algdock.tar.gz /share/apps/algdock/
