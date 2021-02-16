 #!/usr/bin/env python

import argparse
import os
import ROOT
import stackplot_tool as stack

parser = argparse.ArgumentParser(description="Create stack plots for the mu+X skim")
parser.add_argument("resultsDir", help="directory where results histograms are contained")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd()+"/plots/")
parser.add_argument("--l", "--lumi", dest="lumi", help="luminosity", default=0)
parser.add_argument("--t", "--title", dest="label", help="Data label", default="Mu+X Data Skim")
arg = parser.parse_args()

resultsDir = arg.resultsDir
#Check for trailing slash on results dir and delete
if arg.resultsDir.split("/")[-1] =="": resultsDir = arg.resultsDir[:-1]
if not os.path.isdir(resultsDir):
     print "Results directory " + resultsDir + " does not exist!"

outDir = arg.outDir
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
if not os.path.isdir(outDir):
     print "Output directory " + outDir + " does not exist!"
     quit()

#Cross Sections in pb
WJetsCx=61527.6
TTToSemiLeptonicCx=365.34
TTTo2L2NuCx=88.29
DYJetsCx=6077.22
DYJets10to50Cx=18610
WWCx=118.7

mcFiles = []
#mcFIle object hoes Skim File,Results File, Label, Cross Section, Fill Color
#Histograms are stacked in the order they are added here
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/WW/hists/WW_monitor.root","READ"),ROOT.TFile(resultsDir+"/WW/hists/WW_analyzedsplit.root","READ"), "WW",WWCx,7))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_monitor.root","READ"),ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_analyzedsplit.root","READ"), "TTBar2L2Nu",TTTo2L2NuCx,2))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_monitor.root","READ"),ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_analyzedsplit.root","READ"), "TTBarSemiLeptonic",TTToSemiLeptonicCx,5))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/WJets/hists/WJets_monitor.root","READ"), ROOT.TFile(resultsDir+"/WJets/hists/WJets_analyzedsplit.root","READ"), "WJets", WJetsCx,3))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/DYJets/hists/DYJets_monitor.root","READ"),ROOT.TFile(resultsDir+"/DYJets/hists/DYJets_analyzedsplit.root","READ"), "DYJets",DYJetsCx,4))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/DYJets_10to50/hists/DYJets10to50_monitor.root","READ"),ROOT.TFile(resultsDir+"/DYJets_10to50/hists/DYJets10to50_analyzedsplit.root","READ"), "DYJets_10to50",DYJets10to50Cx,6))

dataFile = ROOT.TFile(resultsDir+"/Data/hists/runD_analyzedsplit.root","READ")
info = stack.stackInfo(mcFiles,dataFile, arg.lumi)
info.plotAll("allEvents","Non-Isolated Events",outDir)
info.plotAll("muProbe","Non-Isolated Events, Probe Matched to Muon",outDir)
info.plotAll("nonMuonProbe","Non-Isolated Events, Probe Not Matched to Muon",outDir)

dataFile.Close()

