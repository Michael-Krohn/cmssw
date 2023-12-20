#########################################################################
#Script that puts the histograms that will eventually be used for combine
#into a single file:
#   -Handles the partial and complete disappearance regions
#   -Handles the DY, muPX, and signal predictions
#########################################################################

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
    basedir = '/local/cms/user/revering/datafiles/combineInputs/histograms_23_12_18/'
    masses = ['map0p01','map0p025','map0p05','map0p1','map0p2','map0p4','map0p6','map0p8','map1p0']
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
   parser.add_argument("--i", dest="inDir", help="Directory containing analyzer results", default="/local/cms/user/revering/datafiles/combineInputs/histograms_23_12_18/")
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
 
   dyFile_cutProgress = ROOT.TFile(str(inDir)+"/ZmmSim_EndcapMuFilter.root","READ")
   dyFile = ROOT.TFile(str(inDir)+"/ZmmSim_BaseEventFiltered.root","READ")
   dyShapeFile = ROOT.TFile(str(inDir)+"/DoubleMuonSkim_runA.root","READ")
   muPxFile = ROOT.TFile(str(inDir)+"/offPeakFit.root","READ")
   sigFiles = []
   for filename in os.listdir(inDir):
      f = os.path.join(inDir, filename)
      # checking if it is a file
      if os.path.isfile(f):
          if "DBrem_map" in filename:
             print "filename: ", filename
             print "mass: ", filename.split("_")[1][3:]
             sigFiles.append((ROOT.TFile(str(inDir)+"/"+filename,"READ"),filename.split("_")[1][3:]))

   writeFile = ROOT.TFile(str(outDir)+"/darkPhotonCombineInput.root","RECREATE")

   total = writeFile.mkdir("completeDisappearance")
   partialBDT = writeFile.mkdir("partialDisappearanceBDT")

   total.cd()

   variants = ['','_PUup','_PUdown','_IDup','_IDdown','_ISOup','_ISOdown','_Trigup','_Trigdown']
   sigvariants = ['','_PUup','_PUdown','_IDup','_IDdown','_Trigup','_Trigdown','_ISOup','_ISOdown','_EnBin']
   sigvariantsBDT = ['','_PUup','_PUdown','_IDup','_IDdown','_Trigup','_Trigdown','_ISOup','_ISOdown','_EnBin','_cscReweight']
   biasFactors = {'0p01':5,'0p025':50,'0p05':200,'0p1':300,'0p2':800,'0p4':1200,'0p6':1600,'0p8':2500,'1p0':3200}

   #DY histograms in the complete disappearance region
   dyMassShape = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
   dyMassShapeUncUp = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
   dyMassShapeUncDown = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
   dyNTotalEvents.SetBinContent(1,1451847./0.0669) #Estimation from CRs
   dyNTotalEvents.SetBinError(1,math.sqrt(1451847)/0.0669) #Estimation from CRs
   dyNTotalEvents.Write("totalEvents_DY")

   testDataShape = ROOT.TH1F("testDataShape","testDataShape",20,80,100)
   test = ROOT.TH1F("test","test",20,80,100)
   testDYShape = ROOT.TH1F("testDYShape","testDYShape",20,80,100)
   testDYShapeUp = ROOT.TH1F("testDYShapeUp","testDYShapeUp",20,80,100)
   testDYShapeDown = ROOT.TH1F("testDYShapeDown","testDYShapeDown",20,80,100)
   for variant in variants:
      print "variant: ", variant
      test.Reset()
      testDYShape.Reset()
      testDYShapeUp.Reset()
      testDYShapeDown.Reset()
      dyTree = dyFile.Get("demo/totalDisappear/sigVariables")
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
      test.Scale(59.8/37.4) ###Lumi scale the DY sample
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
            if variant is '':
              testDataShape.SetBinContent(j,dyMassShape.GetBinContent(i))
            j = j + 1
      testDYShape.Write("muonTrackMass_DY"+variant)
      if variant is '':
        testDYShapeUp.Write("muonTrackMass_DY_DYStatUncCompleteup")
        testDYShapeDown.Write("muonTrackMass_DY_DYStatUncCompletedown")
   testDataShape.Write("muonTrackMass_Data")

   #Handling the muPX histograms in the complete disappearance region
   muPXNevents = muPxFile.Get("offPeakShape").Clone()
   muPXNeventsUncUp = muPxFile.Get("offPeakShape").Clone()
   muPXNeventsUncDown = muPxFile.Get("offPeakShape").Clone()
   muPXNeventsSlopeUncUp = muPxFile.Get("offPeakUpSlopeShape").Clone()
   muPXNeventsSlopeUncDown = muPxFile.Get("offPeakDownSlopeShape").Clone()
   muPXScale = 17.662631808415185
   muPXuncertainty = 0.5972814854461866
   muPXNevents.Scale(muPXScale/muPXNevents.Integral())
   muPXNeventsUncUp.Scale((muPXScale+muPXuncertainty)/muPXNeventsUncUp.Integral())
   muPXNeventsUncDown.Scale((muPXScale-muPXuncertainty)/muPXNeventsUncDown.Integral())
   muPXNeventsSlopeUncUp.Scale(muPXScale/muPXNeventsSlopeUncUp.Integral())
   muPXNeventsSlopeUncDown.Scale(muPXScale/muPXNeventsSlopeUncDown.Integral())
   for i in range(1,muPXNevents.GetXaxis().GetNbins() + 1):
       muPXNevents.SetBinError(i,0.001)
       muPXNeventsUncUp.SetBinError(i,0.001)
       muPXNeventsUncDown.SetBinError(i,0.001)
   muPXNevents_fitVariation = ROOT.TH1F("MuPXNevents_fitVariation","",1,0,1)
   muPXNevents_fitVariation.SetBinContent(1,8.2)
   muPXNevents_fitVariation.SetBinError(1,4)
   muPXNeventsUncUp.Write("muonTrackMass_muPX_FitNormUncup")
   muPXNeventsUncDown.Write("muonTrackMass_muPX_FitNormUncdown")
   muPXNeventsSlopeUncUp.Write("muonTrackMass_muPX_Slopeup")
   muPXNeventsSlopeUncDown.Write("muonTrackMass_muPX_Slopedown")
   muPXNevents.Write("muonTrackMass_muPX")

   #Handling the signal histograms in the complete disappearance region
   sigShape = ROOT.TH1F("sigShape","sigShape",20,80,100)
   sigShapeUncUp = ROOT.TH1F("sigShapeUncUp","sigShapeUncUp",20,80,100)
   sigShapeUncDown = ROOT.TH1F("sigShapeUncDown","sigShapeUncDown",20,80,100)
   sigPeakShiftUp = ROOT.TH1F("sigPeakShiftUp","sigPeakShiftUp",20,80,100)
   sigPeakShiftDown = ROOT.TH1F("sigPeakShiftDown","sigPeakShiftDown",20,80,100)
   for sigfile, mass in sigFiles:
      sigPeakShiftUp.Reset()
      sigPeakShiftDown.Reset()
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
      dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/allEvents/cutProgress").GetBinContent(1))
      dbremNTotalEvents_En = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
      sigShapeUncUp.Reset()
      sigShapeUncDown.Reset()
      for variant in sigvariants:
          sigShape.Reset()
          dyMassShape = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
          dyMassShapeUncUp = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
          dyMassShapeUncDown = dyShapeFile.Get("DiMuonFilter/allEvents/MuonsTrackMass").Clone()
          dbremTree = sigfile.Get("demo/totalDisappear/sigVariables")
          test.Reset()
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
          if "0p01" in mass and 'fixBiasing' in sigfile.GetName():
              dyMassShapeUncDown.Scale(10e6*(test.Integral()-binError)/dyMassShapeUncDown.Integral())
          else:
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
                if variant is '':
                  sigPeakShiftUp.SetBinContent(j,dyMassShape.GetBinContent(i-1))
                  sigPeakShiftDown.SetBinContent(j,dyMassShape.GetBinContent(i+1))
                j = j + 1
          if 'fixBiasing' in sigfile.GetName():
              sigShape.Write("muonTrackMass_DBrem_"+mass+variant)
              if variant is '':
                  sigShapeUncUp.Write("muonTrackMass_DBrem_"+mass+"_SigStatUncCompleteup")
                  sigShapeUncDown.Write("muonTrackMass_DBrem_"+mass+"_SigStatUncCompletedown")
          else:
              sigShape.Write("muonTrackMass_DBrem_"+mass+"_ecalBrem"+variant)
              if variant is '':
                  sigShapeUncUp.Write("muonTrackMass_DBrem_"+mass+"_ecalBrem"+"_SigStatUncCompleteup")
                  sigShapeUncDown.Write("muonTrackMass_DBrem_"+mass+"_ecalBrem"+"_SigStatUncCompletedown")
      if 'fixBiasing' in sigfile.GetName():
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
        dbremNTotalEvents_En.Write("totalEvents_DBrem_"+mass+"_EnBin")
        sigPeakShiftUp.Write("muonTrackMass_DBrem_"+mass+"_Peakup")
        sigPeakShiftDown.Write("muonTrackMass_DBrem_"+mass+"_Peakdown")
      else:
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_ecalBrem")
        dbremNTotalEvents_En.Write("totalEvents_DBrem_"+mass+"_ecalBrem_EnBin")
        if variant is '':
            sigShapeUncUp.Write("muonTrackMass_DBrem_"+mass+"_ecalBrem"+"_SigStatUncCompleteup")
            sigShapeUncDown.Write("muonTrackMass_DBrem_"+mass+"_ecalBrem"+"_SigStatUncCompletedown")
        sigPeakShiftUp.Write("muonTrackMass_DBrem_"+mass+"_ecalBrem"+"_Peakup")
        sigPeakShiftDown.Write("muonTrackMass_DBrem_"+mass+"_ecalBrem"+"_Peakdown")

   #Fill partial region bdt fitting region
   partialBDT.cd()
   print "STARTING partialBDT"
   dictList = pickle.load(open('bdtPredictions_addedCut.pkl','rb'))

   dyNTotalEvents = dyFile.Get("demo/allEvents/genEventCount").Clone()
   dyNTotalEvents.SetBinContent(1,(58.6/59.)*(11082440./dyFile_cutProgress.Get("demo/allEvents/cutProgress").GetBinContent(19))*(1.08/4.106)) #SigLumi/SkimLumi * NSelectedProbes / dyFile.Get("demo/allEvents/cutProgress").GetBinContent(19) * (totalDY / peakDY) / (EndcapMuFilter / BAseEventFiler)
   dyNTotalEvents.Write("totalEvents_DY")

   muPXBDTScore = ROOT.TH1F("muPXbdtScore",";BDT Score;Event Weight",1000,0,1)

   for sigfile, mass in sigFiles:
      dbremNTotalEvents = sigfile.Get("demo/allEvents/genEventCount").Clone()
      dbremNTotalEvents.SetBinContent(1,sigfile.Get("demo/noBremEvents/cutProgress").GetBinContent(19))
      dbremNTotalEvents_En = sigfile.Get("demo/allEvents/genEventCount_EnBin").Clone()
      if 'fixBiasing' in sigfile.GetName():
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass)
        dbremNTotalEvents_En.Write("totalEvents_DBrem_"+mass+"_EnBin")
      else:
        dbremNTotalEvents.Write("totalEvents_DBrem_"+mass+"_ecalBrem")
        dbremNTotalEvents_En.Write("totalEvents_DBrem_"+mass+"_ecalBrem"+"_EnBin")

      print "mass: ", mass
      for prediction in dictList:
         if 'name' not in prediction.keys():
            continue
         if prediction['name']==mass:
            for variant in sigvariantsBDT:
               sigBDTScore = ROOT.TH1F("bdtScore"+mass+variant,";BDT Score;Event Weight",1000,0,1)
               dyBDTScore = ROOT.TH1F("dyBDTScore"+mass+variant,";BDT Score;Event Weight",1000,0,1)
               muPXBDTScore.Reset()

               ecalBremString = ''
               if 'fixBiasing' in sigfile.GetName():
                  ecalBremString = ''
               else:
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

               for idx,predValue in enumerate(prediction['sigPreds'+ecalBremString]):
                  sigBDTScore.Fill(predValue,weights[idx])
               for idx,dyPredValue in enumerate(prediction['mcPreds']):
                  dyBDTScore.Fill(dyPredValue,dyWeights[idx])
               for idx,muPXPredValue in enumerate(prediction['muPXPreds']):
                  muPXBDTScore.Fill(muPXPredValue)

               if 'fixBiasing' in sigfile.GetName():
                 sigBDTScore.Write("bdtScore_DBrem_"+mass+variant)
                 if 'EnBin' not in variant:
                   dyBDTScore.Write("bdtScore_DY_"+mass+variant)
               else:
                 sigBDTScore.Write("bdtScore_DBrem_"+mass+"_ecalBrem"+variant)
   muPXBDTScore.Write("bdtScore_muPX")
  
            
   resultsDir = "../MuFilter/results"
   muPxFiles = []

   writeFile.Close()
   dyFile.Close()
   for sigfile, mass in sigFiles:
      sigfile.Close()

if __name__ == "__main__":
   main() 
   
