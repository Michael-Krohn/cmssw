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

JobNumber = JobNumber+1011
from subprocess import call

call(["/local/cms/user/revering/dphoton/CMSSW_10_0_2/src/DarkPhoton/MuGenerator/generate.sh", str(JobNumber), "ZMM14TeV_pythia8_RAW_cfg.py", "ZMM14TeV_pythia8_RECO_cfg.py"])


