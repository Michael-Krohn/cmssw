 #!/usr/bin/env python

import argparse
import os
import ROOT
import selectionvars_tool as tool

parser = argparse.ArgumentParser(description="Create stack plots for the mu+X skim")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd()+"/plots/")
parser.add_argument("--t", "--title", dest="label", help="Data label", default="Mu+X Data Skim")
arg = parser.parse_args()

ROOT.gROOT.SetBatch()
outDir = arg.outDir
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
if not os.path.isdir(outDir):
     print "Output directory " + outDir + " does not exist!"
     quit()

mcFiles = []
#mcFIle object goes Skim File,Results File, Label, Cross Section, Fill Color
#Histograms are stacked in the order they are added here
mcFiles.append(tool.histfile(ROOT.TFile("../ZmmDBrem.root","READ"), "ZmmDBrem, mAp=1.0 GeV",2))
mcFiles.append(tool.histfile(ROOT.TFile("../ZmmSim.root","READ"), "DYJets MC",1))
mcFiles.append(tool.histfile(ROOT.TFile("../ZmmDBrem_0p1.root","READ"), "ZmmDBrem, mAp=0.1 GeV",3))

info = tool.stackInfo(mcFiles)
info.plotAll("allEvents","Private MC, All Events",outDir)


