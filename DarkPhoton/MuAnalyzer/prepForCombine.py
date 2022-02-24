 #!/usr/bin/env python
import argparse
import os
import ROOT
import math
import numpy as np

def calcScaleFactor(name,skimFile,resultFile):
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
   variants = ['','_PUup','_PUdown','_IDup','_IDdown','_Trigup','_Trigdown']
   sigvariants = ['','_PUup','_PUdown','_IDup','_IDdown','_Trigup','_Trigdown','_EnBin']

   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
   dyNTotalEvents.Write("totalEvents_DY")
   for variant in variants:
      dyHcalEnergy = dyFile.Get("demo/totalDisappear/ConeEnergy"+variant).Clone()
      dyHcalEnergy.Write("hcalEnergy_DY"+variant)
      dyNEvents = dyFile.Get("demo/totalDisappear/eventCount"+variant).Clone()
      dyNEvents.Write("nEvents_DY"+variant)
   
   for sigfile, mass in sigFiles:
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_EnBin")
      for variant in sigvariants:
         dbremHcalEnergy = sigfile.Get("demo/totalDisappear/ConeEnergy"+variant).Clone()
         dbremHcalEnergy.Write("hcalEnergy_DBrem_"+mass+variant)
         dbremNEvents = sigfile.Get("demo/totalDisappear/eventCount"+variant).Clone()
         dbremNEvents.Write("nEvents_DBrem_"+mass+variant)
 
   partial.cd()
   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
   dyNTotalEvents.Write("totalEvents_DY")
   for variant in variants:
      dyHcalEnergy = dyFile.Get("demo/partialDisappear/ConeEnergy"+variant).Clone()
      dyHcalEnergy.Write("hcalEnergy_DY"+variant)
      dyNEvents = dyFile.Get("demo/partialDisappear/eventCount"+variant).Clone()
      dyNEvents.Write("nEvents_DY"+variant)
      dyStaDEoverE = dyFile.Get("demo/partialDisappear/StandaloneMuondEoverE"+variant).Clone()
      dyStaDEoverE.Write("staDEoverE_DY"+variant)
   
   for sigfile, mass in sigFiles:
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_EnBin")
      for variant in sigvariants:
         dbremHcalEnergy = sigfile.Get("demo/partialDisappear/ConeEnergy"+variant).Clone()
         dbremHcalEnergy.Write("hcalEnergy_DBrem_"+mass+variant)
         dbremNEvents = sigfile.Get("demo/partialDisappear/eventCount"+variant).Clone()
         dbremNEvents.Write("nEvents_DBrem_"+mass+variant)
         dbremStaDEoverE = sigfile.Get("demo/partialDisappear/StandaloneMuondEoverE"+variant).Clone()    
         dbremStaDEoverE.Write("staDEoverE_DBrem_"+mass+variant)

   resultsDir = "../MuFilter/results"
   muPxFiles = []
   muPxFiles.append((ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_Signal_monitor.root","READ"), ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_SignalRegion_analyzedsplit.root","READ"), "TTBar2L2Nu",88.29))
   muPxFiles.append((ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_Signal_monitor.root","READ"), ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_SignalRegion_analyzedsplit.root","READ"), "TTBarSemiLeptonic",365.34))
   muPxFiles.append((ROOT.TFile(resultsDir+"/WJets/hists/WJets_Signal_monitor.root","READ"),ROOT.TFile(resultsDir+"/WJets/hists/WJets_SignalRegion_analyzedsplit.root","READ"),"WJets",61527.6))

   for monitorFile, resultFile, name, cx in muPxFiles:
      for variant in variants:
         total.cd()
         eff = calcScaleFactor("passSig/eventCount"+variant,monitorFile,resultFile)
         ExpectedNevents = cx*1000*float(59.8)*eff
         NEvents = resultFile.Get("demo/passSig/eventCount"+variant).Clone()
         if(NEvents.Integral()>0):
            NEvents.Scale(ExpectedNevents/NEvents.Integral())
         NEvents.Write("nEvents_"+name+variant)
         partial.cd()
         eff = calcScaleFactor("partialDisappear/StandaloneMuondEoverE"+variant,monitorFile,resultFile)
         ExpectedNevents = cx*1000*float(59.8)*eff
         NEvents = resultFile.Get("demo/partialDisappear/eventCount"+variant).Clone()
         if(NEvents.Integral()>0):
            NEvents.Scale(ExpectedNevents/NEvents.Integral())
         NEvents.Write("nEvents_"+name+variant)
         StaDEoverE = resultFile.Get("demo/partialDisappear/StandaloneMuondEoverE"+variant).Clone()    
         StaDEoverE.Write("staDEoverE_"+name+variant)

   writeFile.Close()
   dyFile.Close()
   for sigfile, mass in sigFiles:
      sigfile.Close()

if __name__ == "__main__":
   main() 
   
