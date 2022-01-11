 #!/usr/bin/env python

import ROOT

class mcFile:
   def __init__(self, skimfile, resultsfile, name, cross, fillcolor):
      self.skim=skimfile
      self.result=resultsfile
      self.label=name
      self.cx=cross
      self.color=fillcolor

class stackInfo:
   def __init__(self, mc, lumi,data=None):
      self.mcFiles=mc
      self.dataFile=data
      self.lumi=lumi
   def plotAll(self,basedir,tag,outDir):
      #Histogram name, X-axis label, Log Scale, Title, Rebin
      self.stack(basedir+"/cutProgress","CutProgress",True,tag,outDir,False)
      self.stack(basedir+"/MuonsTrackMass","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
      self.stack(basedir+"/TaggingMuonEta","Tagging Muon #eta",True,tag,outDir)
      self.stack(basedir+"/TaggingMuonPt","Tagging Muon Pt",False,tag,outDir,True)
      self.stack(basedir+"/TaggingMuonPhi","Tagging Muon #phi",False,tag,outDir)
      self.stack(basedir+"/ProbeTrackEta","Probe Track #eta",True,tag,outDir)
      self.stack(basedir+"/ProbeTrackPt","Probe Track Pt",True,tag,outDir)
      self.stack(basedir+"/ProbeTrackPhi","Probe Track #phi",True,tag,outDir)
      self.stack(basedir+"/ProbeHcalIsolation","Probe Hcal Energy in #Delta R cone of 0.3 (GeV)",True,tag,outDir)
      self.stack(basedir+"/ProbeTrackIsolation","Probe Track Isolation",True,tag,outDir,True)
      self.stack(basedir+"/ProbeEcalIsolation","Ecal Energy within cone (GeV)",False,tag,outDir)
      self.stack(basedir+"/NumberOfMuonsPassingTag","Number of possible tag muons",True,tag,outDir)
      self.stack(basedir+"/TagTrackIsolation","Tag Muon Track Isolation",True,tag,outDir,True)
      self.stack(basedir+"/TagEcalIsolation","Tag Muon ECAL Isolation",True,tag,outDir)
      self.stack(basedir+"/Njets","Number of Jets in Event",False,tag,outDir)
      self.stack(basedir+"/ProbeJetDr","#Delta R between Probe and nearest jet",True,tag,outDir)
      self.stack(basedir+"/ProbeCaloJetE","Energy of Largest Calo Jet within #Delta R of 0.2 to probe", True, tag,outDir)
      self.stack(basedir+"/TagCaloGetE","Energy of Largest Calo Jet within #Delta R of 0.2 to tag", True, tag, outDir, True)
      self.stack(basedir+"/NPV","Number of Primary Vertices",False,tag,outDir)  
      self.stack(basedir+"/TagVtxIndex","Index of Tag Muon Primary Vertex", True, tag, outDir)
      self.stack(basedir+"/AdjustedMuonTrackMass","Adjusted Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
      self.stack(basedir+"/Non Muon Pdg Id","PDG ID of Closest Non-muon to probe",True,tag,outDir,False)
      self.stack(basedir+"/ProbeDCA","Probe Transverse IP",True,tag,outDir,False)
      self.stack(basedir+"/ProbeLongDCA","Probe Longitudinal IP",True,tag,outDir,False)
      self.stack(basedir+"/nExtraElectrons","Extra Electrons",True,tag,outDir,False)
      self.stack(basedir+"/nExtraMuons","Extra Muons",True,tag,outDir,False)
      self.cutFlowText(basedir+"/cutProgress",outDir)

   def calcScaleFactor(self,name,skimFile,resultFile):
      #Use # of muons passing tag for total weight because I know that there was no over/underflow
      eff = resultFile.Get("demo/"+name).Integral()/skimFile.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
      #Adjust efficiency to account for failed jobs
      eff = eff* skimFile.Get("DiMuonFilter/passingEvents/NumberOfMuonsPassingTag").Integral()/resultFile.Get("demo/allEvents/cutProgress").GetBinContent(2)
      return eff

   def stack(self, name, xtitle, log, tag, outDir,rebin=False):
      stackplot = ROOT.THStack("stack","")     
      leg = ROOT.TLegend(0.8,0.6,0.95,0.85)
      leg.SetBorderSize(0)
      hists = []
      NtotalEvents=0
      #nextname = name[:-1]+str(int(name[-1])+1)
      #print nextname
      for mcFile in self.mcFiles:
         eff = self.calcScaleFactor(name,mcFile.skim,mcFile.result)
         #neff = self.calcScaleFactor(nextname,mcFile.skim,mcFile.result)
         Nevents = mcFile.cx*1000*float(self.lumi)*eff
         #Nevents = Nevents+mcFile.cx*1000*float(self.lumi)*neff
         #NtotalEvents=NtotalEvents+Nevents
         hist=(mcFile.result.Get("demo/"+name).Clone())
         #hist.Add(mcFile.result.Get("demo/"+nextname))
         if(hist.Integral()>0):
            hist.Scale(Nevents/hist.Integral())
         hist.SetFillColor(mcFile.color)
         hist.SetLineWidth(0)
         if(rebin):
            hist.Rebin(2)
         hist.SetDirectory(0)
         hists.append(hist)
         leg.AddEntry(hist, mcFile.label,"f")

      for hist in hists:
         stackplot.Add(hist)

      plot= ROOT.TCanvas("stackplot","stackplot",2)
      ROOT.gStyle.SetTextFont(43)
      pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
      pad1.SetTopMargin(0.05)
      pad1.SetBottomMargin(0.1)
      pad1.SetRightMargin(0.03)
      pad1.SetLeftMargin(0.13)
      pad1.Draw()
      pad1.cd()
      if self.dataFile is not None:
         dataHist = self.dataFile.Get("demo/"+name).Clone()
         print dataHist.Integral()
         #dataHist.Add(self.dataFile.Get("demo/"+nextname))
         #print dataHist.Integral()
         #print("Expected " + str(NtotalEvents)+" found " + str(dataHist.Integral())+".")
         dataHist.SetLineWidth(0)
         dataHist.SetLineColor(1)
         dataHist.SetMarkerStyle(21)
         dataHist.SetMarkerSize(0.5)
         if(rebin):
            dataHist.Rebin(2)
          
         dataHist.SetStats(0)
         dataHist.Draw("psame") 
      stackplot.Draw("hist")
      if(log):
         pad1.SetLogy()
         ROOT.gPad.SetLogy()
      if self.dataFile is not None:
         maximum=max(stackplot.GetMaximum(),dataHist.GetMaximum())
         minimum=min(stackplot.GetMinimum(),dataHist.GetMinimum(0))
      else:
         maximum=stackplot.GetMaximum()
         minimum=stackplot.GetMinimum()
      if(minimum==0 and self.dataFile is not None):
         minimum=dataHist.GetMinimum(0)
      if not log:
         minimum=0      
      if log and minimum==0:
         minimum=0.5     
 
      stackplot.SetMaximum(maximum*1.2)
      stackplot.SetMinimum(minimum*0.2)
      stackplot.Draw("hist")
      if self.dataFile is not None:
         dataHist.Draw("psame") 
  
         leg.AddEntry(dataHist, "Data", "p") 
      leg.Draw()
      
      label = ROOT.TPaveText(0.4,0.95,0.99,0.99,"brNDC")
      label.SetTextSize(12)
      label.AddText(tag)
      label.SetFillColor(0)
      label.SetBorderSize(0)
      label.Draw()
      
      stackplot.GetYaxis().SetTitle("Events at "+self.lumi+" fb^{-1}")
      stackplot.GetYaxis().CenterTitle()
      stackplot.GetYaxis().SetTitleSize(20)
      stackplot.GetYaxis().SetTitleFont(43)
      stackplot.GetYaxis().SetTitleOffset(1.35)
      stackplot.GetYaxis().SetLabelFont(43)
      stackplot.GetYaxis().SetLabelSize(15)
      stackplot.GetXaxis().SetLabelFont(43)
      stackplot.GetXaxis().SetLabelSize(15)
      stackplot.GetXaxis().SetTickSize(0.07)
      xlabel = ROOT.TPaveText(0.2,0.0,0.99,0.06,"brNDC")
      xlabel.SetTextSize(20)
      xlabel.AddText(xtitle)
      xlabel.SetFillColor(0)
      xlabel.SetBorderSize(0)
      xlabel.Draw()
      
      lab_x0 = 0.20
      lab_y0 = 0.90
      tag1 = ROOT.TLatex(lab_x0,lab_y0,"CMS")
      tag1.SetNDC()
      tag1.SetTextFont(62)
      tag2 = ROOT.TLatex(lab_x0+0.085, lab_y0, "Internal")
      tag2.SetNDC()
      tag2.SetTextFont(52)
      tag1.SetTextSize(0.04)
      tag2.SetTextSize(0.03)
      tag1.Draw()
      tag2.Draw()
      
      img = ROOT.TImage.Create()
      img.FromPad(plot)
      img.WriteImage(outDir+"/"+name.split("/")[0]+name.split("/")[1]+"_stack.png")
      plot.Print(outDir+"/"+name.split("/")[0]+name.split("/")[1]+"_stack.pdf")
      #plot.Print(outDir+"/"+name.split("/")[0]+name.split("/")[1]+"_stack.C")

   def cutFlowText(self, name, outDir):
      outfile = open(outDir+"/CutProgress.txt","w")
      cutList = [ 'All Events',
                  'Tight ID muon',
                  'Tight PF Iso muon',
                  'Muon pt > 26 GeV and |eta|<2.',
                  'Muon is global',
                  'Track With >20 GeV pt >0.2 dR from muon',
                  'Track |eta|>1.45 and <2.4',
                  'Track is High Purity',
                  'Track Belongs to PV',
                  'Track Longitudinal Impact Parameter <0.05',
                  'Track Transverse Impact Parameter <0.005',
                  'Tag and Probe Opposite Charge',
                  'Tag and Probe valid fitted vertex',
                  'Fitted vtx reduced chi^2<3.0',
                  'Tag Probe Invariant Mass >50  <150 GeV',
                  'Probe track ecal iso > 30 GeV or tk iso>0.1',
                  'No other probes with dR>0.3 pair with track']
     
      hists = [] 
      for mcFile in self.mcFiles:
         eff = self.calcScaleFactor(name,mcFile.skim,mcFile.result)
         #neff = self.calcScaleFactor(nextname,mcFile.skim,mcFile.result)
         Nevents = mcFile.cx*1000*float(self.lumi)*eff
         #Nevents = Nevents+mcFile.cx*1000*float(self.lumi)*neff
         #NtotalEvents=NtotalEvents+Nevents
         hist=(mcFile.result.Get("demo/"+name).Clone())
         #hist.Add(mcFile.result.Get("demo/"+nextname))
         if(hist.Integral()>0):
            hist.Scale(Nevents/hist.Integral())
         hists.append(hist)

      stackplot = hists[0].Clone()    
      for hist in hists:
         if(hist!=hists[0]):
            stackplot.Add(hist)

      outfile.write("Cut Name,WW,TTTo2L2Nu,TTToSemiLeptonic,WJets,DYJets,DYJets_10to50\n")
      for n in range(2,len(cutList)+2):
         outfile.write(cutList[n-2])
         for hist in hists:
            if(hist.GetBinContent(n)>0):
               outfile.write(","+str(hist.GetBinContent(n)))
            else:
               outfile.write(",0")
         outfile.write("\n")

#      outfile.write("Data:\n")
#      dataHist=self.dataFile.Get("demo/allEvents/cutProgress").Clone()
#      for n in range(1,len(cutList)+2):
#         if(dataHist.GetBinContent(n)>0):
#            outfile.write(cutList[n-2]+","+str(dataHist.GetBinContent(n))+"\n")
    
      outfile.close()
