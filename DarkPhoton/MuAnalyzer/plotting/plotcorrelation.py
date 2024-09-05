 #!/usr/bin/env python

import argparse
from time import strftime
import os
import math
import ROOT
from array import array   

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

rfile = ROOT.TFile(arg.rootfile,"UPDATE")
failfail = ROOT.TGraphErrors()
passfail = ROOT.TGraphErrors()
PairCorrelation = ROOT.TMultiGraph()
for i in range(6):
   PairHist = rfile.Get("demo/DepthPairSpectra_{}{}".format(i+1,i+2))
   nbinsx = PairHist.GetNbinsX()
   nbinsy = PairHist.GetNbinsY()
   xfirst = 1
   yfirst = 1
   xupyup = PairHist.Integral(xfirst,-1,yfirst,-1)
   xlowyup = PairHist.Integral(1,xfirst,yfirst,-1)
   xlowylow = PairHist.Integral(1,xfirst,1,yfirst)
   xupylow = PairHist.Integral(xfirst,-1,1,yfirst)
   if xlowylow+xlowyup is not 0:
      a = xlowylow
      b = xlowyup
      failfailerr = math.sqrt(a*pow(1./(a+b)+a/pow(a+b,2.),2.)+b*pow(a/pow(a+b,2.),2.))
      failfail.SetPoint(i,i+1,xlowylow/(xlowylow+xlowyup))
      failfail.SetPointError(i,0,failfailerr)
   if (xupylow+xupyup)>0:
      a = xupylow
      b = xupyup
      passfailerr = math.sqrt(a*pow(1./(a+b)+a/pow(a+b,2.),2.)+b*pow(a/pow(a+b,2.),2.))
      passfail.SetPoint(i,i+1,xupylow/(xupylow+xupyup))
      passfail.SetPointError(i,0,passfailerr)

PairCorrelation.Add(failfail)
failfail.SetTitle("If fail n, Fail n+1")
failfail.SetMarkerColor(4)
failfail.SetMarkerStyle(21)
failfail.SetLineWidth(2)
PairCorrelation.Add(passfail)
PairCorrelation.SetTitle("Depth Pair Correlation")
PairCorrelation.GetXaxis().SetTitle("Depth Pair")
PairCorrelation.GetYaxis().SetTitle("Probability")
passfail.SetTitle("If pass n, Fail n+1")
passfail.SetMarkerColor(3)
passfail.SetMarkerStyle(21)
passfail.SetLineWidth(2)
rfile.cd("demo/")
ROOT.gDirectory.Delete("PairCorrelation;1")
PairCorrelation.Write("PairCorrelation")

rfile.Close()
