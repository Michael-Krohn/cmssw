 #!/usr/bin/env python

import argparse
from time import strftime
import os
import math
import ROOT
from array import array   
from normalize_comp import plotcomp

parser = argparse.ArgumentParser(description="Parse Root File fit eta peaks")
parser.add_argument("rootfile", help="name of files in directory")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd())
arg = parser.parse_args()


outDir = arg.outDir
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
infileType = "."+arg.rootfile.split(".")[-1]

if not os.path.isdir(outDir):
   print "Output directory does not exist!"
   quit()

if not os.path.isfile(arg.rootfile):
   print "Input file does not exist!"
   quit()

if infileType != ".root":
   print "Input file is of incorrect type \"%s\"!"%(infileType)
   quit()

rfile = ROOT.TFile(arg.rootfile,"READ")
for i in range(7):
   MuonHist = rfile.Get(("demo/Layer{}Eta").format(i+7))
   for j in range(MuonHist.GetNbinsY()):
      etahist = ROOT.TH1F()
      
      
   plotcomp(MuonHist, RMuonHist, "Energy (GeV)", "Events", "BeforeHit")

rfile.Close()
