import AlGDock.BindingPMF
import os, shutil, glob

self = AlGDock.BindingPMF.BPMF(dir_dock='dock', dir_cool='cool')
self._replica_exchange('dock')
