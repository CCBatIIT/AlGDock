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

# Copy the compressed file different locations
echo CCB cluster
cp algdock.tar.gz /share/apps/algdock/
rm -rf /share/apps/algdock/AlGDock
tar -xvf /share/apps/algdock/algdock.tar.gz -C /share/apps/algdock

echo OSG Connect cluster
scp algdock.tar.gz daveminh@login.osgconnect.net:~/public/

echo Bridges cluster on XSEDE
scp algdock.tar.gz dminh@bridges.psc.xsede.org:~/software/AlGDock/algdock.tar.gz
ssh dminh@bridges.psc.xsede.org 'tar -xvf ~/software/AlGDock/algdock.tar.gz -C ~/software/AlGDock/'
ssh dminh@bridges.psc.xsede.org 'rm -rf ~/software/AlGDock/pyinstaller'
ssh dminh@bridges.psc.xsede.org 'mv ~/software/AlGDock/AlGDock ~/software/AlGDock/pyinstaller'

echo Comet cluster on XSEDE
scp algdock.tar.gz daveminh@comet.sdsc.xsede.org:~/software/AlGDock/algdock.tar.gz
ssh daveminh@comet.sdsc.xsede.org 'tar -xvf ~/software/AlGDock/algdock.tar.gz -C ~/software/AlGDock/'
ssh daveminh@comet.sdsc.xsede.org 'rm -rf ~/software/AlGDock/pyinstaller'
ssh daveminh@comet.sdsc.xsede.org 'mv ~/software/AlGDock/AlGDock ~/software/AlGDock/pyinstaller'

echo LLNL machines
scp algdock.tar.gz minh1@cab.llnl.gov:/g/g19/minh1/software/AlGDock/algdock.tar.gz
echo Excute the following in /g/g19/minh1/software/AlGDock/
echo tar -xvf /g/g19/minh1/software/AlGDock/algdock.tar.gz -C /g/g19/minh1/software/AlGDock/
echo rm -rf /g/g19/minh1/software/AlGDock/pyinstaller
echo mv /g/g19/minh1/software/AlGDock/AlGDock /g/g19/minh1/software/AlGDock/pyinstaller

# echo Remember to go to '/Users/dminh/Google Drive/Software' and execute:
# echo scp dminh@ccb.tbc.iit.edu:/share/apps/algdock/algdock.tar.gz .
