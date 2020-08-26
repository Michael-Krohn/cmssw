 #!/usr/bin/env python

import argparse
from time import strftime
import os
import math
import ROOT
from array import array   
from muonplotcomp import plotcomp

parser = argparse.ArgumentParser(description="Parse Root Files to create angle and energy distributions")
parser.add_argument("rootfile", help="name of files in directory")
parser.add_argument("rootfile2", help="name of second root file in directory")
parser.add_argument("rootfile3", help="name of third root file in directory")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd())
arg = parser.parse_args()

outDir = arg.outDir
nevts = 0
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
infileType = "."+arg.rootfile.split(".")[-1]

if not os.path.isdir(outDir):
     print "Output directory does not exist!"
     quit()

if not os.path.isfile(arg.rootfile):
     print "Input file does not exist!"
     quit()

if not os.path.isfile(arg.rootfile2):
     print "Input file does not exist!"
     quit()

if infileType != ".root":
     print "Input file is of incorrect type \"%s\"!"%(infileType)
     quit()

infileType = "."+arg.rootfile2.split(".")[-1]

if infileType != ".root":
     print "Input file is of incorrect type \"%s\"!"%(infileType)
     quit()

rfile = ROOT.TFile(arg.rootfile,"READ")
tree = rfile.Get("Events")
evec = ROOT.TLorentzVector()
avec = ROOT.TLorentzVector()
nvec = ROOT.TLorentzVector()

tree.SetBranchAddress("Electron",evec)
tree.SetBranchAddress("APrime",avec)
tree.SetBranchAddress("Nucleus",nvec)

grfile = ROOT.TFile(arg.rootfile2,"READ")
gtree = grfile.Get("Events")
gevec = ROOT.TLorentzVector()
gavec = ROOT.TLorentzVector()
gtree.SetBranchAddress("Electron",gevec)
gtree.SetBranchAddress("APrime",gavec)
entries = tree.GetEntries()
gentries = gtree.GetEntries()
rfile3 = ROOT.TFile(arg.rootfile3,"READ")
tree3 = rfile3.Get("Events")
evec3 = ROOT.TLorentzVector()
avec3 = ROOT.TLorentzVector()
tree3.SetBranchAddress("Electron",evec3)
tree3.SetBranchAddress("APrime",avec3)



tree.GetEntry(0)
gtree.GetEntry(0)
tree3.GetEntry(0)
ebeam1 = 45.0
ebeam2 = 45.0
ebeam3 = 45.0
mass = round(avec.M(),2)

x1 = ROOT.TH1F("M_{A'}=0.01 GeV","mA' = 0.01 GeV",100,0,1.)
x2 = ROOT.TH1F("M_{A'}=0.1 GeV","mA' = 0.1 GeV",100,0,1.)
x3 = ROOT.TH1F("M_{A'}=1.0 GeV","mA' = 1.0 GeV",100,0,1.)

for i in range(entries):
   tree.GetEntry(i)
   m_KE = evec.T()-evec.M()
   x1.Fill(m_KE/(ebeam1-evec.M()))

for i in range(gentries):
   gtree.GetEntry(i)
   tree3.GetEntry(i)
   m_KE = gevec.T()-gevec.M()
   t_KE = evec3.T()-evec3.M()
   Ef = m_KE/(ebeam2-gevec.M())
   x2.Fill(Ef)
   tEf = t_KE/(ebeam3-evec3.M())
   x3.Fill(tEf)

plotcomp(x1,x2,x3,"Outgoing Muon Energy Fraction", "Normalized Rate (1/#sigma d#sigma/dx)", "x", mass, 45.0)

rfile.Close()
