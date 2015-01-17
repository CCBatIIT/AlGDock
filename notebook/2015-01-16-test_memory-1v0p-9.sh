#!/bin/bash

python -m memory_profiler $ALGDOCK --run_type memory_test

# This is code for the memory_test option
#      for k in range(len(self.dock_protocol)):
#        self._set_universe_evaluator(self.dock_protocol[k])

# Loading the grids takes ~212 MiB. The netcdf files take 69 MiB.
#  1546     97.0 MiB      0.0 MiB         fflist.append(self._forceFields['site'])
#  1547    309.6 MiB    212.6 MiB       for scalable in self._scalables:
#  1597     84.7 MiB      0.3 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597     96.7 MiB      0.0 MiB         self.universe, self.universe._forcefield, None, None, None, None)
# Each evaluator with 4 terms takes 43.7 MiB (MM, site, sLJr, ELE)
#  1597    352.9 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    396.6 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    440.3 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    484.0 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    527.7 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    571.4 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    615.1 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    658.8 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    702.5 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    746.2 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    789.9 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    833.6 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    877.3 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    921.0 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597    964.7 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1008.4 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1052.1 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1095.8 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1139.5 MiB     43.7 MiB         self.universe, self.universe._forcefield, None, None, None, None)
# Each evaluator with 7 terms takes ~125 MiB (MM, site, sLJr, sLJa, LJa, LJr, ELE)
#  1597   1336.2 MiB    121.4 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1462.4 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1588.6 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1714.8 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1841.0 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   1967.2 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2093.5 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2219.7 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2345.9 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2472.2 MiB    126.3 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2598.4 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2724.6 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2850.8 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
#  1597   2977.1 MiB    126.2 MiB         self.universe, self.universe._forcefield, None, None, None, None)
# The evaluator with 5 terms takes 82.6 MiB (MM, site, LJa, LJr, ELE)
#  1597   3059.7 MiB     82.6 MiB         self.universe, self.universe._forcefield, None, None, None, None)

# Even with a short cycle of replica exchange, there is no more than 3192.6 MiB.
# Executing a few cycles does not increase it.
# How does the process use ~9000 MiB on the open science grid?

