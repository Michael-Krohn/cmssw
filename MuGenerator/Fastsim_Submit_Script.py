#!/usr/bin/python2

import sys

#eval `scramv1 runtime -sh`
#Needs original file, location of outputfiles, number of events per file, number of files

if sys.argv[1].isdigit:
    JobNumber = int(sys.argv[1])
else:
    sys.exit("Command line argument not a number")
if len(sys.argv) < 1:
    sys.exit("Not enough arguments")

JobNumber = JobNumber+2
from subprocess import call

call(["/local/cms/user/revering/dphoton/CMSSW_10_2_11_patch1/src/DarkPhoton/MuGenerator/generate.sh", str(JobNumber), "ZMM_14TeV_TuneCUETP8M1_cfi_py_GEN_SIM_RECOBEFMIX_DIGI_L1_DIGI2RAW_L1Reco_RECO.py"])


