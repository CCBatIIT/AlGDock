# The example is 1of6

import AlGDock.BindingPMF
import os, shutil, glob
import time

execfile('test_python.py')

for (dock_cycle,cores) in zip(range(3,9),[32,16,8,4,2,1]):
  start_time = time.time()
  self = AlGDock.BindingPMF.BPMF(\
    dir_dock='dock', dir_cool='cool', \
    dock_repX_cycles=dock_cycle, \
    cores=cores, \
    run_type='all')
  del self
  print 'Time for one docking cycle with %d cores: %f s'%(\
    cores, time.time()-start_time)
