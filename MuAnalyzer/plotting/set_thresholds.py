 #!/usr/bin/env python

import argparse
from time import strftime
import os
import math
import ROOT
from array import array   
from plotcomp import plotcomp
from calc_probs import calc_probs
from plotone import plotone

parser = argparse.ArgumentParser(description="Parse Root Files to make plots of muon hit / regular hit with arbitrary thresholds")
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
probs = []
op1miss = []
op1misserr = []
op2miss = []
op2misserr = []
for i in range(7):
   MuonHist = rfile.Get("demo/Layer{}Spectra".format(i+1))
   RMuonHist = rfile.Get("demo/RLayer{}Spectra".format(i+1))
   IntHist = MuonHist.Clone("h3")
   RIntHist = RMuonHist.Clone("h3")
   nbins = IntHist.GetNbinsX()
   TotalMuonInt = MuonHist.Integral(1,nbins+1)
   RTotalMuonInt = RMuonHist.Integral(1,nbins+1)
   MuonPoints, RMuonPoints = array('d'), array('d')
   for j in range(nbins):
      Mint = MuonHist.Integral(j+1,nbins+1)
      Rint = RMuonHist.Integral(j+1,nbins+1)
      Mval = Mint/TotalMuonInt
      IntHist.SetBinContent(j+1,1-Mval)
      MuonPoints.append(Mval)
      if(IntHist.GetBinLowEdge(j+1)==0.2):
         print IntHist.GetBinContent(j+1)
         op2miss.append(IntHist.GetBinContent(j+1))
	 op2misserr.append(math.sqrt(TotalMuonInt-Mint)/TotalMuonInt)
      if(IntHist.GetBinLowEdge(j+1)==0.1):
         op1miss.append(IntHist.GetBinContent(j+1))
	 op1misserr.append(math.sqrt(TotalMuonInt-Mint)/TotalMuonInt)
      Rval = Rint/RTotalMuonInt
      RIntHist.SetBinContent(j+1,Rval)
      RMuonPoints.append(Rval)
      if(RIntHist.GetXaxis().GetBinUpEdge(j+1)==0.1):
         probs.append(RIntHist.GetBinContent(j+1))
   roccurve = ROOT.TGraph(len(MuonPoints),RMuonPoints,MuonPoints)
   plotone(roccurve,"False Positive Rate","True Positive Rate","Depth {}".format(i+1),False,0,1,"Depth{}Roc".format(i+1),0.0,1.0)
   plotcomp(IntHist,RIntHist,"Energy (GeV)", "Fraction of Events Failing Cut","Integrated Depth {}".format(i+1),"Depth{}Int.png".format(i+1),False,False,0,1,0.0,3.0)
   MuonHist.Scale(1/MuonHist.GetEntries(),"width")
   RMuonHist.Scale(1/RMuonHist.GetEntries(),"width")
   plotcomp(MuonHist, RMuonHist, "Energy (GeV)", "Normalized Rate (dN/dE)", "HE hit energy distribution, Depth {}".format(i+1),"Depth_{}.png".format(i+1),False,True,0.01,40,0,2)

op1agmiss = []
op1agmisserr = []
for j in range(6):
   op1agmiss.append(op1miss[j])
   op1agmisserr.append((op1misserr[j]/op1miss[j])**2)
   for k in range(6):
      if k<j:
         op1agmiss[k] = op1agmiss[k]*op1miss[j]
	 op1agmisserr[k] = op1agmisserr[k]+(op1misserr[j]/op1miss[j])**2
for j in range(6):
   op1agmisserr[j] = op1agmiss[j]*math.sqrt(op1agmisserr[j])
op1aghist = ROOT.TH1F("op1aghist","Missing Muon Prob",6,0.5,6.5)
for i in range(6):
   op1aghist.SetBinContent(i+1,op1agmiss[i])
   op1aghist.SetBinError(i+1,op1agmisserr[i])
op2agmiss = []
op2agmisserr = []
for j in range(6):
   op2agmiss.append(op2miss[j])
   op2agmisserr.append((op2misserr[j]/op2miss[j])**2)
   for k in range(6):
      if k<j:
         op2agmiss[k] = op2agmiss[k]*op2miss[j]
	 op2agmisserr[k] = op2agmisserr[k]+(op2misserr[j]/op2miss[j])**2
for j in range(6):
   op2agmisserr[j] = op2agmiss[j]*math.sqrt(op2agmisserr[j])   
op2aghist = ROOT.TH1F("op2aghist","Missing Muon Prob",6,0.5,6.5)
for i in range(6):
   op2aghist.SetBinContent(i+1,op2agmiss[i])
   op2aghist.SetBinError(i+1,op2agmisserr[i])

plotcomp(op1aghist, op2aghist, "Depth of First Missed Hit","Missed Muons/Event","HE Missed Muon Background","MissedBkg.png", False, True,0.000000000001,1)

rfile.Close()
