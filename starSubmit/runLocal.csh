#!/bin/tcsh

starver SL16d

# -- Test on single file:
cons && root -l -b -q StRoot/macros/loadSharedBesLibraries.C StRoot/macros/runPicoBesNetParticleMaker.C++g'("/project/projectdirs/starprod/picodsts/Run14/AuAu/15GeV/Pico16b/P14ii/064/15064001/st_physics_adc_15064001_raw_2000004.picoDst.root", "outputTest", 2, 0, 1)'

# -- Test on file list
#cons && root -l -b -q StRoot/macros/loadSharedBesLibraries.C StRoot/macros/runPicoBesNetParticleMaker.C++g'("lists/test.list", "outputTest", 2, 0, 1)'





