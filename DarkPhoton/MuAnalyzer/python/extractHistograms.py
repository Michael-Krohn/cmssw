####################################################################
#Creates the Combine workspaces. Scales all of the histograms 
#properly and gives them the proper names for Combine
####################################################################

import os,ROOT
import argparse
import math
from array import array

parser = argparse.ArgumentParser(description='option')
parser.add_argument('-i', dest='inputFile')
args = parser.parse_args()

iFile = args.inputFile

masses = [
"DBrem_2p0",
"DBrem_1p0",
"DBrem_0p8",
"DBrem_0p6",
"DBrem_0p4",
"DBrem_0p2",
"DBrem_0p1",
"DBrem_0p05",
"DBrem_0p025",
"DBrem_0p01",
] 


f_origin = ROOT.TFile(iFile)
channels = f_origin.GetListOfKeys()

DYxsec = 398 #### For the private RECO sample
#DYxsec = 1583 #### For the private RECO sample, for the old version

lumi = 59800.

biasFactors = {'0p01':5,'0p025':50,'0p05':200,'0p1':300,'0p2':800,'0p4':1200,'0p6':1600,'0p8':2500,'1p0':3200}

invariantMassBinning = [80,82,84,86,88,90,92,94,96,98,100]
ainvariantMassBinning = array("d",invariantMassBinning)

for channel in channels:

    print channel.GetName()

    directory = f_origin.Get(channel.GetName())
    histoKeys = directory.GetListOfKeys()


    for mass in masses:

      #hcalEnergyBinning = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.84,0.88,0.9,0.92,0.94,0.96,0.98,1.0]
      hcalEnergyBinning = [0.98,1.0]
      ahcalEnergyBinning = array("d",hcalEnergyBinning)

      out = ROOT.TFile('CombineInput/'+channel.GetName()+'_'+mass+'.root','RECREATE')
      DYnominal = ROOT.TH1F()

      for key in histoKeys:
	if (mass in key.GetName()) or ('DY' in key.GetName()) or ('muPX' in key.GetName()) or ('Data' in key.GetName()):

	  #### Partial Disappearnace Region ####
	  if ("partialDisappearanceBDT" in channel.GetName() and "bdtScore" in key.GetName()):
	    if mass[6:] not in key.GetName() and 'muPX' not in key.GetName():continue

	    totalEventsName = key.GetName().replace("bdtScore","totalEvents")
	    totalEvents = ROOT.TH1F()

	    #### Grab Up/Down histos ####
	    if 'DY' in totalEventsName:
	      totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName[0:14])
	    elif 'DBrem' in totalEventsName:
	      totalEventsName = totalEventsName.split("_")[0]+'_'+totalEventsName.split("_")[1]+'_'+totalEventsName.split("_")[2]
              totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
	      
            h = f_origin.Get(channel.GetName()+'/'+key.GetName())
	    if 'DY' in key.GetName():
	      h.Scale(totalEvents.Integral())
	    elif mass in key.GetName():
              h.Scale((58.6/59.)*11082440.*1.08/biasFactors[mass[6:]])#SigLumi/SkimLumi * NSelectedProbes *totalDY / peakDY / biasFactor ###With epsilon^2 = 1

	    h = h.Rebin(len(ahcalEnergyBinning) - 1, h.GetName(), ahcalEnergyBinning)

            if 'DY' in key.GetName():
                h.SetName(key.GetName()[9:].replace("up","Up").replace("down","Down").replace('_'+key.GetName().split('_')[2],""))
            else:
                h.SetName(key.GetName()[9:].replace("up","Up").replace("down","Down"))
	    if 'DBrem' in key.GetName() and 'EnBin' in key.GetName():
	      h.SetName(key.GetName()[9:]+'Up')
            if 'DY' in key.GetName() and 'cscReweight' in key.GetName():
              h.SetName(h.GetName()+'Down')

            out.cd()
	    print "Saving: ", h.GetName()
            if h.GetName() == 'DY':
                DYnominal = h.Clone()
            h.Write()

	    ###Temporary until we have data
            if 'muPX' in key.GetName() and 'up' not in key.GetName() and 'down' not in key.GetName():
              fakeData = h.Clone("data_obs")
              fakeData.SetDirectory(0)
              #fakeData.Write()
            elif 'DY' in key.GetName() and 'up' not in key.GetName() and 'down' not in key.GetName() and 'KDE' not in key.GetName() and 'csc' not in key.GetName():
              fakeData.Add(h)
              fakeData.Write()

            ###Creating UP cscReweight histo
            if 'DY' in key.GetName() and 'cscReweight' in key.GetName():
                cscReweightUp = h.Clone("DY_cscReweightUp")
                cscReweightUp.SetBinContent(1,abs(h.Integral() - DYnominal.Integral()) + DYnominal.Integral())
                cscReweightUp.Write()

	    ###Setting nominal DBrem histo to EnBinDown
	    if 'DBrem' in key.GetName() and 'up' not in key.GetName() and 'down' not in key.GetName() and 'EnBin' not in key.GetName():
              h.Scale(1.001)
	      h.SetName(key.GetName()[9:]+'_EnBinDown')
	      h.Write()


	  #### Complete disappearance region ####
	  elif ("complete" in channel.GetName() and "muonTrackMass" in key.GetName()):
	    totalEventsName = key.GetName().replace("muonTrackMass","totalEvents")
	    totalEvents = ROOT.TH1F()

	    #### Grab Up/Down histos ####
	    if ('up' in totalEventsName or 'down' in totalEventsName) and 'DY' in totalEventsName:
	      totalEventsName = totalEventsName[:14]
	      totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
	    elif 'DY' in totalEventsName:
	      totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
	    elif ('up' in totalEventsName or 'down' in totalEventsName or 'EnBin' in totalEventsName) and 'DBrem' in totalEventsName:
	      totalEventsName = totalEventsName[:21]
	      totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
	    elif 'DBrem' in totalEventsName:
              totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
            elif 'Data' in totalEventsName:
              totalEvents = f_origin.Get(channel.GetName()+'/totalEvents_DY')  
            h = f_origin.Get(channel.GetName()+'/'+key.GetName())
	    if 'DY' in key.GetName():
	      h.Scale(DYxsec*lumi/totalEvents.GetBinContent(1))
	    elif mass in key.GetName():
              h.Scale(DYxsec*lumi/biasFactors[mass[6:]]) ### With epsilon^2 = 1

            h.SetName(key.GetName()[14:].replace("up","Up").replace("down","Down"))
	    if 'DBrem' in key.GetName() and 'EnBin' in key.GetName():
	      h.SetName(key.GetName()[14:]+'Up')
            out.cd()
	    h = h.Rebin(len(ainvariantMassBinning) - 1, h.GetName(), ainvariantMassBinning)
            h.Write()

            ###Temporary until we have data
            if 'Data' in key.GetName() and 'up' not in key.GetName() and 'down' not in key.GetName():
              h.Scale(DYxsec*lumi/totalEvents.GetBinContent(1))
              fakeData = h.Clone("data_obs")
            elif 'muPX' in key.GetName() and 'up' not in key.GetName() and 'down' not in key.GetName():
              fakeData.Add(h)
              fakeData.Write()

	    ###Setting nominal DBrem histo to EnBinDown
	    if 'DBrem' in key.GetName() and 'up' not in key.GetName() and 'down' not in key.GetName() and 'EnBin' not in key.GetName():
	      h.SetName(key.GetName()[14:]+'_EnBinDown')
              h.Scale(1.001)
	      h = h.Rebin(len(ainvariantMassBinning) - 1, h.GetName(), ainvariantMassBinning)
	      h.Write()

      out.Close()
