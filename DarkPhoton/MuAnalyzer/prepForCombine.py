 #!/usr/bin/env python
import argparse
import os
import ROOT
import math
import numpy as np

def calcScaleFactor(name,skimFile,resultFile)
   eff = resultFile.Get("demo/"+name).Integral()/skimFile.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
   eff = eff* skimFile.Get("DiMuonFilter/passingEvents/NumberOfMuonsPassingTag").Integral()/resultFile.Get("demo/allEvents/cutProgress").GetBinContent(2)
   return eff

def main():
   
   parser = argparse.ArgumentParser(description="Fit the background fraction as a function of isolation")
   parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd()+"/histograms/")
   parser.add_argument("--i", dest="inDir", help="Directory containing analyzer results", default=os.getcwd()+"/histograms/")
   arg = parser.parse_args()
   outDir = arg.outDir
   #Check for trailing slash on ouput dir and delete
   if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
   if not os.path.isdir(outDir):
        print "Output directory " + outDir + " does not exist!"
        quit()
   inDir = arg.inDir
   #Check for trailing slash on input dir and delete
   if arg.inDir.split("/")[-1] == "": inDir = arg.inDir[:-1]
   if not os.path.isdir(inDir):
        print "Input directory " + inDir + " does not exist!"
        quit()
 
   dyFile = ROOT.TFile(str(inDir)+"/ZmmSim.root","READ")
   sigFiles = []
   for filename in os.listdir(inDir):
      f = os.path.join(inDir, filename)
      # checking if it is a file
      if os.path.isfile(f):
          if "DBrem" in filename:
             sigFiles.append((ROOT.TFile(str(inDir)+"/"+filename,"READ"),filename.split("_")[-1].split(".")[0]))

   writeFile = ROOT.TFile(str(outDir)+"/darkPhotonCombineInput.root","RECREATE")

   total = writeFile.mkdir("completeDisappearance")
   partial = writeFile.mkdir("partialDisappearance")

   total.cd()
   dyHcalEnergy = dyFile.Get("demo/totalDisappear/ConeEnergy").Clone()
   dyHcalEnergy.Write("hcalEnergy_DY")
   dyNTotalEvents = dyFile.Get("demo/totalDisappear/eventCount").Clone()
   dyNTotalEvents.Write("totalEvents_DY")
   dyNEvents = dyFile.Get("demo/totalDisappear/eventCount").Clone()
   dyNEvents.Write("nEvents_DY")
   
   for sigfile, mass in sigFiles:
      dbremHcalEnergy = sigfile.Get("demo/totalDisappear/ConeEnergy").Clone()
      dbremHcalEnergy.Write("hcalEnergy_DBrem_"+mass)
      dbremNTotalEvents = sigfile.Get("demo/totalDisappear/eventCount").Clone()
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
      dbremNEvents = sigfile.Get("demo/totalDisappear/eventCount").Clone()
      dbremNEvents.Write("nEvents_DBrem_"+mass)
 
   partial.cd()
   dyHcalEnergy = dyFile.Get("demo/partialDisappear/ConeEnergy").Clone()
   dyHcalEnergy.Write("hcalEnergy_DY")
   dyNTotalEvents = dyFile.Get("demo/partialDisappear/eventCount").Clone()
   dyNTotalEvents.Write("totalEvents_DY")
   dyNEvents = dyFile.Get("demo/partialDisappear/eventCount").Clone()
   dyNEvents.Write("nEvents_DY")
   dyStaDEoverE = dyFile.Get("demo/partialDisappear/StandaloneMuondEoverE").Clone()
   dyStaDEoverE.Write("staDEoverE_DY")
   
   for sigfile, mass in sigFiles:
      dbremHcalEnergy = sigfile.Get("demo/partialDisappear/ConeEnergy").Clone()
      dbremHcalEnergy.Write("hcalEnergy_DBrem_"+mass)
      dbremNTotalEvents = sigfile.Get("demo/partialDisappear/eventCount").Clone()
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
      dbremNEvents = sigfile.Get("demo/partialDisappear/eventCount").Clone()
      dbremNEvents.Write("nEvents_DBrem_"+mass)
      dbremStaDEoverE = sigfile.Get("demo/partialDisappear/StandaloneMuondEoverE").Clone()    
      dbremStaDEoverE.Write("staDEoverE_DBrem_"+mass)

   resultsDir = "../MuFilter/results"
   muPxFiles = []
   muPxFiles.append((ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_Signal_monitor.root","READ"), ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_SignalRegion_analyzedsplit.root","READ"), "TTBar2L2Nu"))
   muPxFiles.append((ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_Signal","READ"), ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_SignalRegion_analyzedsplit.root","READ"), "TTBarSemiLeptonic"))
   muPxFiles.append((ROOT.TFile(resultsDir+"/WJets/hists/WJets_Signal_monitor.root","READ"),ROOT.TFile(resultsDir+"/WJets/hists/WJets_SignalRegion_analyzedsplit.root","READ"),"WJets"))

   for monitorFile, resultFile, name in muPxFiles:
      eff = calcScaleFactor("partialDisappear/StandaloneMuondEoverE",monitorFile,resultFile)
      ExpectedNevents = cx*1000*float(lumi)*eff
      NEvents = resultFile.Get("demo/partialDisappear/eventCount").Clone()
      if(NEvents.Integral()>0):
         NEvents.Scale(ExpectedNevents/NEvents.Integral())
      NEvents.Write("nEvents_"+name)
      StaDEoverE = resultFile.Get("demo/partialDisappear/StandaloneMuondEoverE").Clone()    
      StaDEoverE.Write("staDEoverE_"+name)

   writeFile.Close()
   dyFile.Close()
   for sigfile, mass in sigFiles:
      sigfile.Close()

if __name__ == "__main__":
   main() 
   
