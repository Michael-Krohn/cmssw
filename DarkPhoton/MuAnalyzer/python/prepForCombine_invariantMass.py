 #!/usr/bin/env python
 #!/usr/bin/env python
import argparse
import os
import ROOT
import pickle
import math
import pandas
import numpy as np

def calcScaleFactor(name,skimFile,resultFile):
   eff = resultFile.Get("demo/"+name).Integral()/skimFile.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
   eff = eff* skimFile.Get("DiMuonFilter/passingEvents/NumberOfMuonsPassingTag").Integral()/resultFile.Get("demo/allEvents/cutProgress").GetBinContent(2)
   return eff

def calcSigContamination():
    basedir = '/local/cms/user/revering/datafiles/combineInputs/histograms_23_6_21/'
    masses = ['map0p2','map0p4','map0p6','map0p8','map1p0']
    rates = []
    fitScale = 0.296
    offPeakLumi = 14.73
    offPeakEvents = 16.45
    for mass in masses:
        offPeak = ROOT.TFile(basedir+'ZmmDBrem_offPeakCR_'+mass+'.root')
        sig = ROOT.TFile(basedir+'ZmmDBrem_'+mass+'_fixBiasing.root')
        normrate = sig.Get('demo/totalDisappear/eventCount').Integral()/sig.Get('demo/allEvents/cutProgress').GetBinContent(2)
        offPeakRate = offPeak.Get('demo/totalDisappear/eventCount').Integral()/offPeak.Get('demo/allEvents/cutProgress').GetBinContent(2)
        rates.append(offPeakRate)
        print(offPeakRate/normrate)
    sigContRates = zip(masses,rates)

def main():
   #calcSigContamination() 
   parser = argparse.ArgumentParser(description="Fit the background fraction as a function of isolation")
   parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd()+"/histograms/")
   parser.add_argument("--i", dest="inDir", help="Directory containing analyzer results", default="/local/cms/user/revering/datafiles/combineInputs/histograms_23_6_21/")
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
 
   #dyFile = ROOT.TFile(str(inDir)+"/ZmmSim_EndcapMuFilter.root","READ")
   dyFile = ROOT.TFile(str(inDir)+"/ZmmSim_BaseEventFiltered.root","READ")
   dyShapeFile = ROOT.TFile(str(inDir)+"/DoubleMuonSkim_runA.root","READ")
   muPxFile = ROOT.TFile(str(inDir)+"/offPeakFit.root","READ")
   sigFiles = []
   for filename in os.listdir(inDir):
      f = os.path.join(inDir, filename)
      # checking if it is a file
      if os.path.isfile(f):
          #if "DBrem_map" in filename and "ecalBrem" in filename:
          #if "DBrem_map" in filename and "fixBiasing" in filename:
          if "DBrem_map" in filename:
             print "filename: ", filename
             sigFiles.append((ROOT.TFile(str(inDir)+"/"+filename,"READ"),filename.split("_")[1][-3:]))
             #sigFiles.append((ROOT.TFile(str(inDir)+"/"+filename,"READ"),filename.split("_")[-1].split(".")[0]))

   writeFile = ROOT.TFile(str(outDir)+"/darkPhotonCombineInput.root","RECREATE")

   total = writeFile.mkdir("completeDisappearance")
   partial = writeFile.mkdir("partialDisappearance")
   partialStaDphi = writeFile.mkdir("partialDisappearanceLargeStaDPhi")
   partialBDT = writeFile.mkdir("partialDisappearanceBDT")

   total.cd()
   variants = ['','_PUup','_PUdown','_IDup','_IDdown','_ISOup','_ISOdown','_Trigup','_Trigdown']
   #variants = ['']
   sigvariants = ['','_PUup','_PUdown','_IDup','_IDdown','_Trigup','_Trigdown','_ISOup','_ISOdown','_EnBin']
   #sigvariants = ['']
   sigvariantsBDT = ['','_PUup','_PUdown','_IDup','_IDdown','_Trigup','_Trigdown','_ISOup','_ISOdown','_EnBin','_cscReweight']
   biasFactors = {'0p2':800,'0p4':1200,'0p6':1600,'0p8':2500,'1p0':3200}
   #biasFactors = {'0p1':2000,'0p2':4000,'0p4':8000,'1p0':40000,'2p0':80000,'4p0':160000}
  
   dyMassShape = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
   dyMassShapeUncUp = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
   dyMassShapeUncDown = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
   dyNTotalEvents.SetBinContent(1,1451847./0.0669) #Estimation from CRs
   dyNTotalEvents.SetBinError(1,math.sqrt(1451847)/0.0669) #Estimation from CRs
   dyNTotalEvents.Write("totalEvents_DY")
   testDataShape = ROOT.TH1F("testDataShape","testDataShape",20,80,100)
   for variant in variants:
      dyTree = dyFile.Get("demo/totalDisappear/sigVariables")
      test = ROOT.TH1F("test","test",20,80,100)
      testDYShape = ROOT.TH1F("testDYShape","testDYShape",20,80,100)
      testDYShapeUp = ROOT.TH1F("testDYShapeUp","testDYShapeUp",20,80,100)
      testDYShapeDown = ROOT.TH1F("testDYShapeDown","testDYShapeDown",20,80,100)
      for entry in dyTree:
         if entry.hcalIso < 15:
             if variant is '':
               test.Fill(entry.invMass,entry.EventWeight)
             elif variant is '_PUup':
               test.Fill(entry.invMass,entry.PUupWeight)
             elif variant is '_PUdown':
               test.Fill(entry.invMass,entry.PUdownWeight)
             elif variant is '_IDup':
               test.Fill(entry.invMass,entry.IDupWeight)
             elif variant is '_IDdown':
               test.Fill(entry.invMass,entry.IDdownWeight)
             elif variant is '_ISOup':
               test.Fill(entry.invMass,entry.ISOupWeight)
             elif variant is '_ISOdown':
               test.Fill(entry.invMass,entry.ISOdownWeight)
             elif variant is '_Trigup':
               test.Fill(entry.invMass,entry.TrigUpWeight)
             elif variant is '_Trigdown':
               test.Fill(entry.invMass,entry.TrigDownWeight)
             elif variant is '_EnBin':
               test.Fill(entry.invMass,entry.EnBinWeight)
      test.Scale(59.8/37.4)
      binError = 0
      for i in range(1,test.GetXaxis().GetNbins() + 1):
          binError = binError + test.GetBinError(i)**2
      binError = math.sqrt(binError)
      dyMassShape.Scale(test.Integral()/dyMassShape.Integral())
      dyMassShapeUncUp.Scale((test.Integral()+binError)/dyMassShapeUncUp.Integral())
      dyMassShapeUncDown.Scale((test.Integral()-binError)/dyMassShapeUncDown.Integral())
      j = 1
      setBinError = 0.
      for i in range(1,dyMassShape.GetXaxis().GetNbins() + 1):
          if dyMassShape.GetXaxis().GetBinCenter(i) > 80 and dyMassShape.GetXaxis().GetBinCenter(i) < 100:
            testDYShape.SetBinContent(j,dyMassShape.GetBinContent(i))
            testDYShape.SetBinError(j,dyMassShape.GetBinError(i))
            testDYShapeUp.SetBinContent(j,dyMassShapeUncUp.GetBinContent(i))
            testDYShapeUp.SetBinError(j,dyMassShapeUncUp.GetBinError(i))
            testDYShapeDown.SetBinContent(j,dyMassShapeUncDown.GetBinContent(i))
            testDYShapeDown.SetBinError(j,dyMassShapeUncDown.GetBinError(i))
            #if binError/math.sqrt(20.) < dyMassShape.GetBinContent(i):
            #  testDYShape.SetBinError(j,binError/math.sqrt(20.))
            #  setBinError = setBinError + (binError/math.sqrt(20.))**2
            #else:
            #  testDYShape.SetBinError(j,dyMassShape.GetBinContent(i))
            #  setBinError = setBinError + dyMassShape.GetBinContent(i)**2
            if variant is '':
              print "FILLING testData"
              testDataShape.SetBinContent(j,dyMassShape.GetBinContent(i))
              #testDataShape.SetBinError(j,binError/math.sqrt(20.))
            j = j + 1
      print "binError: ", binError
      print "setBinError: ", math.sqrt(setBinError)
      testDYShape.Write("nEvents_DY"+variant)
      if variant is '':
        testDYShapeUp.Write("nEvents_DY_StatUncup")
        testDYShapeDown.Write("nEvents_DY_StatUncdown")
   testDataShape.Write("nEvents_Data")
      #dyMassShape.Write("nEvents_DY"+variant)
#      test.Write("nEvents_DY"+variant)
   muPXNevents = muPxFile.Get("offPeakShape").Clone()
   muPXNeventsUncUp = muPxFile.Get("offPeakShape").Clone()
   muPXNeventsUncDown = muPxFile.Get("offPeakShape").Clone()
   muPXScale = 17.701145585223223
   muPXuncertainty = 0.6079596221266597
   muPXNevents.Scale(muPXScale/muPXNevents.Integral())
   muPXNeventsUncUp.Scale((muPXScale+muPXuncertainty)/muPXNeventsUncUp.Integral())
   muPXNeventsUncDown.Scale((muPXScale-muPXuncertainty)/muPXNeventsUncDown.Integral())
   for i in range(1,muPXNevents.GetXaxis().GetNbins() + 1):
       muPXNevents.SetBinError(i,0)
   #    muPXNevents.SetBinError(i,muPXuncertainty/math.sqrt(20.))
   muPXNevents_fitVariation = ROOT.TH1F("MuPXNevents_fitVariation","",1,0,1)
   muPXNevents_fitVariation.SetBinContent(1,8.2)
   muPXNevents_fitVariation.SetBinError(1,4)
#   muPXNevents_fitVariation.Write("nEvents_muPX_fitVariationUp")
#   muPXNevents_fitVariation.Write("nEvents_muPX_fitVariationDown")
   muPXNeventsUncUp.Write("nEvents_muPX_Fitup")
   muPXNeventsUncDown.Write("nEvents_muPX_Fitdown")
   for variant in variants:
      muPXNevents.Write("nEvents_muPX"+variant)
   for sigfile, mass in sigFiles:
      sigPeakShiftUp = ROOT.TH1F("sigPeakShiftUp","sigPeakShiftUp",20,80,100)
      sigPeakShiftDown = ROOT.TH1F("sigPeakShiftDown","sigPeakShiftDown",20,80,100)
      if 'fixBiasing' in sigfile.GetName():
        print "mass: ", mass
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
        print "BinContent: ", sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1)
        dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1))
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_EnBin")
        sigShapeUncUp = ROOT.TH1F("sigShapeUncUp","sigShapeUncUp",20,80,100)
        sigShapeUncDown = ROOT.TH1F("sigShapeUncDown","sigShapeUncDown",20,80,100)
        for variant in sigvariants:
           sigShape = ROOT.TH1F("sigShape","sigShape",20,80,100)
           print "variant: ", variant
           dyMassShape = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
           dyMassShapeUncUp = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
           dyMassShapeUncDown = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
           dbremTree = sigfile.Get("demo/totalDisappear/sigVariables")
           test = ROOT.TH1F("test","test",20,80,100)
           for entry in dbremTree:
              if entry.hcalIso < 15:
                if variant is '':
                  test.Fill(entry.invMass,entry.EventWeight)
                elif variant is '_PUup':
                  test.Fill(entry.invMass,entry.PUupWeight)
                elif variant is '_PUdown':
                  test.Fill(entry.invMass,entry.PUdownWeight)
                elif variant is '_IDup':
                  test.Fill(entry.invMass,entry.IDupWeight)
                elif variant is '_IDdown':
                  test.Fill(entry.invMass,entry.IDdownWeight)
                elif variant is '_ISOup':
                  test.Fill(entry.invMass,entry.ISOupWeight)
                elif variant is '_ISOdown':
                  test.Fill(entry.invMass,entry.ISOdownWeight)
                elif variant is '_Trigup':
                  test.Fill(entry.invMass,entry.TrigUpWeight)
                elif variant is '_Trigdown':
                  test.Fill(entry.invMass,entry.TrigDownWeight)
                elif variant is '_EnBin':
                  test.Fill(entry.invMass,entry.EnBinWeight)
           dyMassShape.Scale(test.Integral()/dyMassShape.Integral())
           binError = 0
           for i in range(1,test.GetXaxis().GetNbins() + 1):
              binError = binError + test.GetBinError(i)**2
           binError = math.sqrt(binError)
           dyMassShapeUncUp.Scale((test.Integral()+binError)/dyMassShapeUncUp.Integral())
           dyMassShapeUncDown.Scale((test.Integral()-binError)/dyMassShapeUncDown.Integral())
           j = 1
           setBinError = 0.
           for i in range(1,dyMassShape.GetXaxis().GetNbins() + 1):
              if dyMassShape.GetXaxis().GetBinCenter(i) > 80 and dyMassShape.GetXaxis().GetBinCenter(i) < 100:
                sigShape.SetBinContent(j,dyMassShape.GetBinContent(i))
                sigShape.SetBinError(j,dyMassShape.GetBinError(i))
                sigShapeUncUp.SetBinContent(j,dyMassShapeUncUp.GetBinContent(i))
                sigShapeUncUp.SetBinError(j,dyMassShapeUncUp.GetBinError(i))
                sigShapeUncDown.SetBinContent(j,dyMassShapeUncDown.GetBinContent(i))
                sigShapeUncDown.SetBinError(j,dyMassShapeUncDown.GetBinError(i))
                #if binError/math.sqrt(20.) < dyMassShape.GetBinContent(i):
                #  sigShape.SetBinError(j,binError/math.sqrt(20.))
                #  setBinError = setBinError + (binError/math.sqrt(20.))**2
                #else:
                #  sigShape.SetBinError(j,dyMassShape.GetBinContent(i))
                #  setBinError = setBinError + dyMassShape.GetBinContent(i)**2
                if variant is '':
                  sigPeakShiftUp.SetBinContent(j,dyMassShape.GetBinContent(i-1))
                  #sigPeakShiftUp.SetBinError(j,binError/math.sqrt(20.))
                  sigPeakShiftDown.SetBinContent(j,dyMassShape.GetBinContent(i+1))
                  #sigPeakShiftDown.SetBinError(j,binError/math.sqrt(20.))
                j = j + 1
           sigShape.Write("nEvents_DBrem_"+mass+variant)
           if variant is '':
             sigShapeUncUp.Write("nEvents_DBrem_"+mass+"_StatUncup")
             sigShapeUncDown.Write("nEvents_DBrem_"+mass+"_StatUncdown")
        sigPeakShiftUp.Write("nEvents_DBrem_"+mass+"_Peakup")
        sigPeakShiftDown.Write("nEvents_DBrem_"+mass+"_Peakdown")
           #dbremNEvents.Scale(1/2000.)
           #dbremNEvents.Write("nEvents_DBrem_"+mass+variant)
      else:
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
        dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1))
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_ecalBrem")
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_ecalBrem_EnBin")
        sigShapeUncUp = ROOT.TH1F("sigShapeUncUp","sigShapeUncUp",20,80,100)
        sigShapeUncDown = ROOT.TH1F("sigShapeUncDown","sigShapeUncDown",20,80,100)
        for variant in sigvariants:
           sigShape = ROOT.TH1F("sigShape","sigShape",20,80,100)
           dbremTree = sigfile.Get("demo/totalDisappear/sigVariables")
           test = ROOT.TH1F("test","test",20,80,100)
           dyMassShape = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
           dyMassShapeUncUp = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
           dyMassShapeUncDown = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
           for entry in dbremTree:
              if entry.hcalIso < 15:
                if variant is '':
                  test.Fill(entry.invMass,entry.EventWeight)
                elif variant is '_PUup':
                  test.Fill(entry.invMass,entry.PUupWeight)
                elif variant is '_PUdown':
                  test.Fill(entry.invMass,entry.PUdownWeight)
                elif variant is '_IDup':
                  test.Fill(entry.invMass,entry.IDupWeight)
                elif variant is '_IDdown':
                  test.Fill(entry.invMass,entry.IDdownWeight)
                elif variant is '_ISOup':
                  test.Fill(entry.invMass,entry.ISOupWeight)
                elif variant is '_ISOdown':
                  test.Fill(entry.invMass,entry.ISOdownWeight)
                elif variant is '_Trigup':
                  test.Fill(entry.invMass,entry.TrigUpWeight)
                elif variant is '_Trigdown':
                  test.Fill(entry.invMass,entry.TrigDownWeight)
                elif variant is '_EnBin':
                  test.Fill(entry.invMass,entry.EnBinWeight)
           dyMassShape.Scale(test.Integral()/dyMassShape.Integral())
           binError = 0
           for i in range(1,test.GetXaxis().GetNbins() + 1):
              binError = binError + test.GetBinError(i)**2
           binError = math.sqrt(binError)
           dyMassShapeUncUp.Scale((test.Integral()+binError)/dyMassShapeUncUp.Integral())
           dyMassShapeUncDown.Scale((test.Integral()-binError)/dyMassShapeUncDown.Integral())
           print "mass: ", mass
           print "sigfile.GetName(): ", sigfile.GetName()
           print "dyMassShape.Integral(): ", dyMassShape.Integral()
           j = 1
           setBinError = 0.
           for i in range(1,dyMassShape.GetXaxis().GetNbins() + 1):
              if dyMassShape.GetXaxis().GetBinCenter(i) > 80 and dyMassShape.GetXaxis().GetBinCenter(i) < 100:
                sigShape.SetBinContent(j,dyMassShape.GetBinContent(i))
                sigShape.SetBinError(j,dyMassShape.GetBinError(i))
                sigShapeUncUp.SetBinContent(j,dyMassShapeUncUp.GetBinContent(i))
                sigShapeUncUp.SetBinError(j,dyMassShapeUncUp.GetBinError(i))
                sigShapeUncDown.SetBinContent(j,dyMassShapeUncDown.GetBinContent(i))
                sigShapeUncDown.SetBinError(j,dyMassShapeUncDown.GetBinError(i))
                #if binError/math.sqrt(20.) < dyMassShape.GetBinContent(i):
                #  sigShape.SetBinError(j,binError/math.sqrt(20.))
                #  setBinError = setBinError + (binError/math.sqrt(20.))**2
                #else:
                #  sigShape.SetBinError(j,dyMassShape.GetBinContent(i))
               #   setBinError = setBinError + dyMassShape.GetBinContent(i)**2
                if variant is '':
                  sigPeakShiftUp.SetBinContent(j,dyMassShape.GetBinContent(i-1))
               #   sigPeakShiftUp.SetBinError(j,binError/math.sqrt(20.))
                  sigPeakShiftDown.SetBinContent(j,dyMassShape.GetBinContent(i+1))
               #   sigPeakShiftDown.SetBinError(j,binError/math.sqrt(20.))
                j = j + 1
           sigShape.Write("nEvents_DBrem_"+mass+"_ecalBrem"+variant)
           if variant is '':
             sigShapeUncUp.Write("nEvents_DBrem_"+mass+"_ecalBrem"+"_StatUncup")
             sigShapeUncDown.Write("nEvents_DBrem_"+mass+"_ecalBrem"+"_StatUncdown")
        sigPeakShiftUp.Write("nEvents_DBrem_"+mass+"_ecalBrem"+"_Peakup")
        sigPeakShiftDown.Write("nEvents_DBrem_"+mass+"_ecalBrem"+"_Peakdown")

   #Fill partial region signal files 
   partial.cd()
   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
   dyNTotalEvents.SetBinContent(1,9819550)
   #dyNTotalEvents.SetBinContent(1,dyFile.Get("demo/allEvents/cutProgress").GetBinContent(1))
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
        dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1))
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_EnBin")
        for variant in sigvariants:
           dbremHcalEnergy = sigfile.Get("demo/partialDisappear/ConeEnergy"+variant).Clone()
           dbremHcalEnergy.Write("hcalEnergy_DBrem_"+mass+variant)
           dbremNEvents = sigfile.Get("demo/partialDisappear/eventCount"+variant).Clone()
           #dbremNEvents.Scale(1/2000.)
           dbremNEvents.Write("nEvents_DBrem_"+mass+variant)
           dbremStaDEoverE = sigfile.Get("demo/partialDisappear/StandaloneMuondEoverE"+variant).Clone()    
           dbremStaDEoverE.Write("staDEoverE_DBrem_"+mass+variant)

   #Fill partial region large dPhi signal
   partialStaDphi.cd()
   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
   dyNTotalEvents.SetBinContent(1,9819550)
   #dyNTotalEvents.SetBinContent(1,dyFile.Get("demo/allEvents/cutProgress").GetBinContent(1))
   dyNTotalEvents.Write("totalEvents_DY")
   for variant in variants:
      dyHcalEnergy = dyFile.Get("demo/partialDisappear_largeCSCDefl/ConeEnergy"+variant).Clone()
      dyHcalEnergy.Write("hcalEnergy_DY"+variant)
      dyNEvents = dyFile.Get("demo/partialDisappear_largeCSCDefl/eventCount"+variant).Clone()
      dyNEvents.Write("nEvents_DY"+variant)
      dyStaDEoverE = dyFile.Get("demo/partialDisappear_largeCSCDefl/StandaloneMuondEoverE"+variant).Clone()
      dyStaDEoverE.Write("staDEoverE_DY"+variant)
   
   for sigfile, mass in sigFiles:
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
      dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1))
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
      dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_EnBin")
      for variant in sigvariants:
         dbremHcalEnergy = sigfile.Get("demo/partialDisappear_largeCSCDefl/ConeEnergy"+variant).Clone()
         dbremHcalEnergy.Write("hcalEnergy_DBrem_"+mass+variant)
         dbremNEvents = sigfile.Get("demo/partialDisappear_largeCSCDefl/eventCount"+variant).Clone()
         #dbremNEvents.Scale(1/2000.)
         dbremNEvents.Write("nEvents_DBrem_"+mass+variant)
         dbremStaDEoverE = sigfile.Get("demo/partialDisappear_largeCSCDefl/StandaloneMuondEoverE"+variant).Clone()    
         dbremStaDEoverE.Write("staDEoverE_DBrem_"+mass+variant)

   resultsDir = "../MuFilter/results"
   muPxFiles = []
#   muPxFiles.append((ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_Signal_monitor.root","READ"), ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_SignalRegion_analyzedsplit.root","READ"), "TTBar2L2Nu",88.29))
#   muPxFiles.append((ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_Signal_monitor.root","READ"), ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_SignalRegion_analyzedsplit.root","READ"), "TTBarSemiLeptonic",365.34))
#   muPxFiles.append((ROOT.TFile(resultsDir+"/WJets/hists/WJets_Signal_monitor.root","READ"),ROOT.TFile(resultsDir+"/WJets/hists/WJets_SignalRegion_analyzedsplit.root","READ"),"WJets",61527.6))

   for monitorFile, resultFile, name, cx in muPxFiles:
      for variant in variants:
         total.cd()
         eff = calcScaleFactor("passSig/eventCount"+variant,monitorFile,resultFile)
         ExpectedNevents = cx*1000*float(59.8)*eff
         NEvents = resultFile.Get("demo/passSig/eventCount"+variant).Clone()
         if(NEvents.Integral()<0):
            print("Broken file! Name is: "+name+".\n")
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
         partialStaDphi.cd()
         eff = calcScaleFactor("partialDisappearLargeStaDphi/StandaloneMuondEoverE"+variant,monitorFile,resultFile)
         ExpectedNevents = cx*1000*float(59.8)*eff
         NEvents = resultFile.Get("demo/partialDisappearLargeStaDphi/eventCount"+variant).Clone()
         if(NEvents.Integral()>0):
            NEvents.Scale(ExpectedNevents/NEvents.Integral())
         NEvents.Write("nEvents_"+name+variant)
         StaDEoverE = resultFile.Get("demo/partialDisappearLargeStaDphi/StandaloneMuondEoverE"+variant).Clone()    
         StaDEoverE.Write("staDEoverE_"+name+variant)
         hcalEnergy = dyFile.Get("demo/partialDisappear_largeCSCDefl/ConeEnergy"+variant).Clone()
         hcalEnergy.Write("hcalEnergy_"+name+variant)

   #Fill partial region bdt fitting region
   partialBDT.cd()
   #dictList = pickle.load(open('darkPhotonBDT/bdtPredictions.pkl','rb'))
   dictList = pickle.load(open('bdtPredictions1p0.pkl','rb'))
   #dictList = pickle.load(open('bdtPredictions1p0_fixBiasing.pkl','rb'))
   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
#   dyNTotalEvents.SetBinContent(1,9819550)
   dyNTotalEvents.SetBinContent(1,dyFile.Get("demo/allEvents/cutProgress").GetBinContent(1)*14.5) #14.5 is the inverse of the filter efficiency
   #dyNTotalEvents.SetBinContent(1,dyFile.Get("demo/allEvents/cutProgress").GetBinContent(1))
   dyNTotalEvents.Write("totalEvents_DY")
   for sigfile, mass in sigFiles:
      if 'fixBiasing' in sigfile.GetName():
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
        dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1))
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_EnBin")
      else:
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
        dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1))
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_ecalBrem")
        dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_ecalBrem"+"_EnBin")

      print "mass: ", mass
      for prediction in dictList:
         if 'name' not in prediction.keys():
            continue
         if prediction['name']==mass:
            for variant in sigvariantsBDT:
               sigBDTScore = ROOT.TH1F("bdtScore"+mass+variant,";BDT Score;Event Weight",1000,0,1)
               dyBDTScore = ROOT.TH1F("dyBDTScore"+mass+variant,";BDT Score;Event Weight",1000,0,1)
               muPXBDTScore = ROOT.TH1F("muPXbdtScore",";BDT Score;Event Weight",1000,0,1)
               ecalBremString = ''
               if 'fixBiasing' in sigfile.GetName():
                  ecalBremString = ''
               else:
                  print "ECAL BREM"
                  ecalBremString = '_ecalBrem'
               if variant=='':
                  weights = prediction['EventWeight'+ecalBremString]
                  dyWeights = dictList[-1]['EventWeight']
               elif 'Trigup' in variant:
                  weights = prediction['TrigUpWeight'+ecalBremString]
                  dyWeights = dictList[-1]['TrigUpWeight']
               elif 'Trigdown' in variant:
                  weights = prediction['TrigDownWeight'+ecalBremString]
                  dyWeights = dictList[-1]['TrigDownWeight']
               elif 'cscReweight' in variant:
                  weights = prediction['EventWeight_cscReweight'+ecalBremString]
                  dyWeights = dictList[-1]['EventWeight_cscReweight']
               else:
                  weights = prediction[variant.split("_")[-1]+"Weight"+ecalBremString]
                  dyWeights = dictList[-1][variant.split("_")[-1]+"Weight"]
               summedSigEventWeight = 0
               print len(prediction['sigPreds'+ecalBremString])
               print len(weights)
               for idx,predValue in enumerate(prediction['sigPreds'+ecalBremString]):
                  #if predValue > 0.99: 
                      #summedSigEventWeight += weights[idx]
                  sigBDTScore.Fill(predValue,weights[idx])
               adyWeights = []
               dyBDT = []
               for idx,dyPredValue in enumerate(prediction['mcPreds']):
                  dyBDTScore.Fill(dyPredValue,dyWeights[idx])
                  #if dyWeights[idx] > 2: print dyWeights[idx]
                  #if dyPredValue > 0.6:
                  #      print "dyWeights[idx]: ", dyWeights[idx]
#                  if dyPredValue < 0.1: continue
#                  if dyPredValue < 0.05: continue
                  dyBDT.append(dyPredValue)
                  adyWeights.append(dyWeights[idx])
               print "variant: ", variant
               print "summedSigEventWeight: ", summedSigEventWeight
               for idx,muPXPredValue in enumerate(prediction['muPXPreds']):
                  muPXBDTScore.Fill(muPXPredValue)

               #print "mass: ", mass
               #print "variant: ", variant
               #print "len(weights): ", len(dyWeights)
               #print "dyWeights: ", dyWeights
               #print "Number of DY entries: ", dyBDTScore.GetEntries()
               #print "DY Integral: ", dyBDTScore.Integral()
               dyBDTScoreKDE = dyBDTScore.Clone("dyBDTScoreKDE"+mass+variant)
               dyBDTnp = np.array(dyBDT)
               dyWeightsnp = np.array(adyWeights)
#               dykde = ROOT.TKDE(len(dyBDT),dyBDTnp,dyWeightsnp,0.1,1)
#              dykde = ROOT.TKDE(len(dyBDT),dyBDTnp,dyWeightsnp,0.05,1)
               dykde = ROOT.TKDE(len(dyBDT),dyBDTnp,dyWeightsnp,0,1)
               dykde.SetTuneFactor(15)
               #print "dyBDTScore.GetBinContent(1): ", dyBDTScore.GetBinContent(1)
               #print "dyBDTScore.GetBinContent(2): ", dyBDTScore.GetBinContent(2)
               for ibin in range(1,dyBDTScoreKDE.GetNbinsX() + 1):
#                  if ibin == 1:
#                  if ibin < 11:
#                     dyBDTScoreKDE.SetBinContent(ibin,dyBDTScore.GetBinContent(ibin))
#                     dyBDTScoreKDE.SetBinError(ibin,dyBDTScore.GetBinError(ibin))
#                  else:
                     dyBDTScoreKDE.SetBinContent(ibin,dykde.GetValue(dyBDTScoreKDE.GetBinCenter(ibin)))
#                     if ibin > 990: 
#                        print "dykde.GetValuer(dyBDTScoreKDE.GetBinCenter(ibin)): ", dykde.GetValue(dyBDTScoreKDE.GetBinCenter(ibin))
#                        print "ibin: ", ibin
                     dyBDTScoreKDE.SetBinError(ibin,dykde.GetError(dyBDTScoreKDE.GetBinCenter(ibin)))
                     
               dyBDTScoreKDE.Scale(dyBDTScore.Integral()/dyBDTScoreKDE.Integral())                
#               c1 = ROOT.TCanvas("c1","c1",600,600)
#               graph = dykde.GetGraphWithErrors(100,0,1)
               #print "dykde.Integral(): ", graph.GetErrorY(99)
#               fl = dykde.GetLowerFunction(0.684, 100,0,1)
#               fu = dykde.GetUpperFunction(0.684, 100,0,1)
#               fl_2sigma = dykde.GetLowerFunction(0.654)
#               fu_2sigma = dykde.GetUpperFunction(0.954)
#               dyBDTScoreKDE.SetLineColor(ROOT.kGreen)
#               dyBDTScoreKDE.Draw("hist e same")
#               dyBDTScoreKDE.GetYaxis().SetRangeUser(0.0001,8000)
#               f_dykde = dykde.GetFunction(100,0,1) 
#               f_dykde.SetLineColor(ROOT.kBlack)
#               f_dykde.Draw("same") 
#               fl.Draw("same")
#               fu.Draw("same")
#               fl_2sigma.SetLineColor(ROOT.kViolet)
#               fu_2sigma.SetLineColor(ROOT.kViolet)
#               fl_2sigma.Draw("same")
#               fu_2sigma.Draw("same")
#               dyBDTScore.Draw("hist e same")
#               c1.SetLogy()
#               c1.SaveAs("bdtScoreRawKDE_DY_"+mass+variant+".pdf")

               dyBDTScoreKDE.SetEntries(dyBDTScore.GetEntries())
               if 'fixBiasing' in sigfile.GetName():
                 sigBDTScore.Write("bdtScore_DBrem_"+mass+variant)
                 if 'EnBin' not in variant:
                   dyBDTScore.Write("bdtScore_DY_"+mass+variant)
                   dyBDTScoreKDE.Write("bdtScoreKDE_DY_"+mass+variant)
                 if 'cscReweight' not in variant:
                   dyNEvents = dyFile.Get("demo/partialDisappear/eventCount"+variant).Clone()
                   sigNEvents = sigfile.Get("demo/partialDisappear/eventCount"+variant).Clone()
                   dyNEvents.Scale(prediction['dyEff'])
                   sigNEvents.Scale(prediction['sigEff'])
                   dyNEvents.Write("nEvents_DY_"+mass+variant)
                   sigNEvents.Write("nEvents_DBrem_"+mass+variant)
                 else:
                   print "variant!!!!!!!", variant
                   dyNEvents = dyFile.Get("demo/partialDisappear/eventCount").Clone()
                   sigNEvents = sigfile.Get("demo/partialDisappear/eventCount").Clone()
                   dyNEvents.Scale(prediction['dyEff'])
                   sigNEvents.Scale(prediction['sigEff'])
                   dyNEvents.Write("nEvents_DY_"+mass+variant)
                   sigNEvents.Write("nEvents_DBrem_"+mass+variant)
               else:
                 sigBDTScore.Write("bdtScore_DBrem_"+mass+"_ecalBrem"+variant)
                 if 'cscReweight' not in variant:
                   sigNEvents = sigfile.Get("demo/partialDisappear/eventCount"+variant).Clone()
                   sigNEvents.Scale(prediction['sigEff'])
                   sigNEvents.Write("nEvents_DBrem_"+mass+"_ecalBrem"+variant)
                 else:
                   print "variant!!!!!!!", variant
                   sigNEvents = sigfile.Get("demo/partialDisappear/eventCount").Clone()
                   sigNEvents.Scale(prediction['sigEff'])
                   sigNEvents.Write("nEvents_DBrem_"+mass+"_ecalBrem"+variant)                    
   muPXBDTScore.Write("bdtScore_muPX")
  
            
   resultsDir = "../MuFilter/results"
   muPxFiles = []
#   muPxFiles.append((ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_Signal_monitor.root","READ"), ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_SignalRegion_analyzedsplit.root","READ"), "TTBar2L2Nu",88.29))
##   muPxFiles.append((ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_Signal_monitor.root","READ"), ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_SignalRegion_analyzedsplit.root","READ"), "TTBarSemiLeptonic",365.34))
##   muPxFiles.append((ROOT.TFile(resultsDir+"/WJets/hists/WJets_Signal_monitor.root","READ"),ROOT.TFile(resultsDir+"/WJets/hists/WJets_SignalRegion_analyzedsplit.root","READ"),"WJets",61527.6))

   for monitorFile, resultFile, name, cx in muPxFiles:
      for variant in variants:
         total.cd()
         eff = calcScaleFactor("passSig/eventCount"+variant,monitorFile,resultFile)
         ExpectedNevents = cx*1000*float(59.8)*eff
         NEvents = resultFile.Get("demo/passSig/eventCount"+variant).Clone()
         if(NEvents.Integral()<0):
            print("Broken file! Name is: "+name+".\n")
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
         partialStaDphi.cd()
         eff = calcScaleFactor("partialDisappearLargeStaDphi/StandaloneMuondEoverE"+variant,monitorFile,resultFile)
         ExpectedNevents = cx*1000*float(59.8)*eff
         NEvents = resultFile.Get("demo/partialDisappearLargeStaDphi/eventCount"+variant).Clone()
         if(NEvents.Integral()>0):
            NEvents.Scale(ExpectedNevents/NEvents.Integral())
         NEvents.Write("nEvents_"+name+variant)
         StaDEoverE = resultFile.Get("demo/partialDisappearLargeStaDphi/StandaloneMuondEoverE"+variant).Clone()    
         StaDEoverE.Write("staDEoverE_"+name+variant)
         hcalEnergy = dyFile.Get("demo/partialDisappear_largeCSCDefl/ConeEnergy"+variant).Clone()
         hcalEnergy.Write("hcalEnergy_"+name+variant)

   writeFile.Close()
   dyFile.Close()
   for sigfile, mass in sigFiles:
      sigfile.Close()

if __name__ == "__main__":
   main() 
   
