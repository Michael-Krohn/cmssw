##########################################################
#Script that merges the HCAL and ECAL Brem signal sample
#histograms and saves all other histograms in the output
#file
##########################################################

import os,ROOT
import argparse

parser = argparse.ArgumentParser(description='option')
parser.add_argument('-i', dest='inputFile')
args = parser.parse_args()

iFile = args.inputFile

f_origin = ROOT.TFile(iFile)
channels = f_origin.GetListOfKeys()

out = ROOT.TFile('Modified/'+args.inputFile,'RECREATE')

for channel in channels:
    print "channel.GetName(): ", channel.GetName()

    directory = f_origin.Get(channel.GetName())
    histoKeys = directory.GetListOfKeys()

    i = 0
    DBremHistos = {}
    print "channel: ", channel.GetName()
    if channel.GetName() == 'partialDisappearance': continue
    if channel.GetName() == 'partialDisappearanceLargeStaDPhi': continue
    out.mkdir(channel.GetName())
    print "len(histoKeys): ", len(histoKeys)
    for key in histoKeys:
        print "key: ", key.GetName()
        ####Handling the complete disappearance region  signal samples
        if 'complete' in channel.GetName() and 'DBrem' in key.GetName() and "muonTrackMass" in key.GetName():
                h = ROOT.TH1F()
                firstSample = True
                if key.GetName().replace('_ecalBrem','') not in DBremHistos:
                    DBremHistos[key.GetName().replace('_ecalBrem','')] = f_origin.Get(channel.GetName()+'/'+key.GetName())
                else:
                    h = f_origin.Get(channel.GetName()+'/'+key.GetName())
                    firstSample = False
                    
                totalEventsName = key.GetName().replace("muonTrackMass","totalEvents")
                totalEvents = ROOT.TH1F()
                if ('up' in totalEventsName or 'down' in totalEventsName or 'EnBin' in totalEventsName) and 'DBrem' in totalEventsName:
                    if 'ecalBrem' in key.GetName():
                        totalEventsName = totalEventsName.split("_")[0]+'_'+totalEventsName.split("_")[1]+'_'+totalEventsName.split("_")[2]+'_'+totalEventsName.split("_")[3]
                    else:
                        totalEventsName = totalEventsName.split("_")[0]+'_'+totalEventsName.split("_")[1]+'_'+totalEventsName.split("_")[2]
                    totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
                elif 'DBrem' in totalEventsName:
                    totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
                
                if firstSample:
                    DBremHistos[key.GetName().replace('_ecalBrem','')].Scale(1./totalEvents.GetBinContent(1))
                else:
                    h.Scale(1./totalEvents.GetBinContent(1))
                    DBremHistos[key.GetName().replace('_ecalBrem','')].Add(h)
                    DBremHistos[key.GetName().replace('_ecalBrem','')].SetName(key.GetName().replace('_ecalBrem',''))
        ####Handling the partial disappearance region  signal samples
        elif 'partialDisappearanceBDT' in channel.GetName() and 'DBrem' in key.GetName() and "bdtScore" in key.GetName():
                h = ROOT.TH1F()                
                firstSample = True
                if key.GetName().replace('_ecalBrem','') not in DBremHistos:
                    DBremHistos[key.GetName().replace('_ecalBrem','')] = f_origin.Get(channel.GetName()+'/'+key.GetName())
                else:
                    h = f_origin.Get(channel.GetName()+'/'+key.GetName())
                    firstSample = False

                totalEventsName = key.GetName().replace("bdtScore","totalEvents")
                totalEvents = ROOT.TH1F()
                if ('up' in totalEventsName or 'down' in totalEventsName or 'EnBin' in totalEventsName or 'cscReweight' in totalEventsName) and 'DBrem' in totalEventsName:
                    if 'ecalBrem' not in totalEventsName:
                        totalEventsName = totalEventsName.split("_")[0]+'_'+totalEventsName.split("_")[1]+'_'+totalEventsName.split("_")[2]
                    else:
                        totalEventsName = totalEventsName.split("_")[0]+'_'+totalEventsName.split("_")[1]+'_'+totalEventsName.split("_")[2]+'_'+totalEventsName.split("_")[3]
                    totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)
                elif 'DBrem' in totalEventsName:
                    totalEvents = f_origin.Get(channel.GetName()+'/'+totalEventsName)

                if firstSample:
                    DBremHistos[key.GetName().replace('_ecalBrem','')].Scale(1./totalEvents.GetBinContent(1))
                else:
                    h.Scale(1./totalEvents.GetBinContent(1))
                    DBremHistos[key.GetName().replace('_ecalBrem','')].Add(h)
                    DBremHistos[key.GetName().replace('_ecalBrem','')].SetName(key.GetName().replace('_ecalBrem',''))                    
        ####Handling all other histograms 
	else:
          h = f_origin.Get(channel.GetName()+'/'+key.GetName())
	  h.SetName(key.GetName())

          out.cd()
	  out.cd(channel.GetName())
          print "else statement Writing: ", h.GetName()
          h.Write()
   

    for key in DBremHistos:
        print "Writing key: ", key
        print "Writing: ", DBremHistos[key].GetName()
        out.cd(channel.GetName())
        DBremHistos[key].Write()
