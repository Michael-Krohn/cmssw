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

#JobNumber = JobNumber+301
from subprocess import call

call(["/data/cmszfs1/user/revering/dphoton/DBremAdd/CMSSW_10_2_6/src/DarkPhoton/MuGenerator/ZMM14TeV_SignalGen.sh", str(JobNumber), "ZMM14TeV_pythia8_GEN_cfg.py"])


