 #!/usr/bin/env python

import ROOT
import copy
import math

class mcFile:
   def __init__(self, skimfile, resultsfile, name, cross, fillcolor):
      self.skim=skimfile
      self.result=resultsfile
      self.label=name
      self.cx=cross
      self.color=fillcolor

class stackInfo:
   def __init__(self, mc, data, lumi):
      self.mcFiles=mc
      self.dataFile=data
      self.lumi=lumi
   def plotAll(self,basedir,tag,outDir):
      #Histogram name, X-axis label, Log Scale, Title, Output Directory, Rebin
      self.stack(basedir+"/MuonsTrackMass","Invariant Mass of Tag/Probe Pair (GeV)",True,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass0","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass1","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass2","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass3","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass4","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass5","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass6","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass7","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass8","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalMuonsTrackMass9","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass0","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass1","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass2","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass3","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass4","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass5","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass6","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass7","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass8","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/HcalFractionMuonsTrackMass9","Invariant Mass of Tag/Probe Pair (GeV)",False,tag,outDir,True)

#      self.stack(basedir+"/ProbeCaloJetDr","#Delta R (probe, Calo Jet)",False,tag,outDir,True)
#      self.stack(basedir+"/TrackVsGenMuonE","",False,tag,outDir,True)
#      self.stack(basedir+"/GlobalVsGenMuonE","",False,tag,outDir,True)

      self.stack(basedir+"/TaggingMuonEta","Tagging Muon #eta",True,tag,outDir,True)
      self.stack(basedir+"/TaggingMuonPt","Tagging Muon Pt (GeV)",False,tag,outDir,True)
      self.stack(basedir+"/TaggingMuonPhi","Tagging Muon #phi",False,tag,outDir,True)
      self.stack(basedir+"/ProbeTrackEta","Probe Track #eta",True,tag,outDir, True)
      self.stack(basedir+"/ProbeTrackPt","Probe Track Pt (GeV)",True,tag,outDir,True)
      self.stack(basedir+"/ProbeTrackPhi","Probe Track #phi",True,tag,outDir, True)
      self.stack(basedir+"/ProbeHcalIsolation","Probe Hcal Energy in #Delta R cone of 0.2 (GeV)",True,tag,outDir, True)
      self.stack(basedir+"/ProbeTrackIsolation","Probe Track Isolation (GeV)",True,tag,outDir,True)
      self.stack(basedir+"/ProbeEcalIsolation","Probe Ecal Energy within cone 0.4 (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/ProbeEcalIsolation_dR0p35","Probe Ecal Energy within cone 0.35 (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/ProbeEcalIsolation_dR0p3","Probe Ecal Energy within cone 0.3 (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/ProbeEcalIsolation_dR0p2","Probe Ecal Energy within cone 0.2 (GeV)",False,tag,outDir,True)
#      self.stack(basedir+"/TagTrackIsolation","Tag Muon Track Isolation",True,tag,outDir,True)
#      self.stack(basedir+"/TagEcalIsolation","Tag Muon ECAL Isolation",False,tag,outDir)
#      self.stack(basedir+"/ProbeCaloJetDr","#Delta R between Probe and nearest Calo jet",True,tag,outDir)
#      self.stack(basedir+"/GlobalMuonDR","dR between Global muon and track",True,tag,outDir)
#      self.stack(basedir+"/GlobalMuonE","Energy of Global muon matched to track",True,tag,outDir,True)
#      self.stack(basedir+"/GlobalMuonEta","#eta of Global muon matched to track",True,tag,outDir,True)
#      self.stack(basedir+"/GlobalMuonPhi","#phi of Global muon matched to track",True,tag,outDir,True)
#      self.stack(basedir+"/GlobalMuonPt","p_{T} of Global muon matched to track",True,tag,outDir,True)
#      self.stack(basedir+"/CaloJetHad","HcalE of Calo Jet",True,tag,outDir,True)
#      self.stack(basedir+"/ProbePlusHcalE","Energy of probe + calo jet HcalE",True,tag,outDir)
  

   def getRatio(self, hist, reference):
        ratio = hist.Clone("%s_ratio"%hist.GetName())
        ratio.SetDirectory(0)
        ratio.SetLineColor(hist.GetLineColor())
        for xbin in xrange(1,reference.GetNbinsX()+1):
                ref = reference.GetBinContent(xbin)
                val = hist.GetBinContent(xbin)

                refE = reference.GetBinError(xbin)
                valE = hist.GetBinError(xbin)

                try:
                        ratio.SetBinContent(xbin, val/ref)
                        ratio.SetBinError(xbin, math.sqrt( (val*refE/(ref**2))**2 + (valE/ref)**2 ))
                except ZeroDivisionError:
                        #ratio.SetBinContent(xbin, 1.0)
                        ratio.SetBinError(xbin, 0.0)

        return ratio

   def TOTerror(self, hmc, ratio ):
      hmc.Sumw2()
      den1 = hmc.Clone ("den1");
      den2 = hmc.Clone ("den2");

      nvar = hmc.GetNbinsX();

      x1 = []
      y1 = []
      exl1 = []
      eyl1= []
      exh1= []
      eyh1= []

      for km in range(1,nvar+1):
        delta = hmc.GetBinError(km)
        den1.SetBinError(km,0)

      # ratio from variation and nominal
      ratiop = hmc.Clone("ratiop");
      ratiom = hmc.Clone("ratiom");

      ratiop.Divide(den1);
      ratiom.Divide(den1);
      #den1.Divide(ratiop)
      #den2.Divide(ratiom)

      return ratiop;


   def calcScaleFactor(self,name,skimFile,resultFile):
      #Use # of muons passing tag for total weight because I know that there was no over/underflow
      #eff = resultFile.Get("demo/"+name).Integral()/skimFile.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
      eff = resultFile.Get("demo/"+name).Integral()/skimFile.Get("demo/allEventsNoCuts/eventCount").Integral()
      return eff

   def stack(self, name, xtitle, log, tag, outDir,rebin=False):
      stackplot = ROOT.THStack("stack","")     
      leg = ROOT.TLegend(0.6,0.75,0.85,0.88)
      leg.SetBorderSize(0)
      leg.SetTextSize(0.035)
      leg.SetTextFont(42)
      leg2 = ROOT.TLegend(0.15,0.75,0.4,0.88)
      leg2.SetBorderSize(0)
      leg2.SetTextSize(0.035)
      leg2.SetTextFont(42)
      hists = []
      count = 0
      for mcFile in self.mcFiles:
         eff = self.calcScaleFactor(name,mcFile.skim,mcFile.result)
	 print "nevents: ", eff
	 print "mcFile xsec: ", mcFile.cx
         Nevents = mcFile.cx*1000*float(self.lumi)*eff
         hist=(mcFile.result.Get("demo/"+name).Clone())
         if(hist.Integral()>0):
            hist.Scale(Nevents/hist.Integral())
         hist.SetFillColor(mcFile.color)
         hist.SetLineWidth(0)
         if(rebin):
            hist.Rebin(2)
            #hist.Rebin(4)
         hist.SetDirectory(0)
         hists.append(hist)
	 if count < 2:
           leg2.AddEntry(hist, mcFile.label,"f")
	 else:
           leg.AddEntry(hist, mcFile.label,"f")
	 count += 1

      dataHist = self.dataFile.Get("demo/"+name).Clone()
      for hist in hists:
	 print "MC hist integral(): ", hist.Integral()
	 hist.Sumw2()
         #hist.Scale(dataHist.Integral()/hist.Integral())
         stackplot.Add(hist)

      plot= ROOT.TCanvas("stackplot","stackplot",2)
      ROOT.gStyle.SetTextFont(43)
      pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1)
      pad2 = ROOT.TPad('pad2','pad2',0,0.0,1.0,0.3)
      pad1.SetBottomMargin(0)
      pad2.SetTopMargin(0.)
      pad2.SetBottomMargin(0.35)

      pad1.SetFillStyle(4000)
      pad1.SetFrameFillStyle(1000)
      pad1.SetFrameFillColor(0)
      pad2.SetFillStyle(4000)
      pad2.SetFrameFillStyle(1000)
      pad2.SetFrameFillColor(0)
      pad1.Draw()
      pad2.Draw()
      pad1.cd()
     
      dataHist.Sumw2()
      dataHist.SetLineWidth(0)
      dataHist.SetLineColor(1)
      dataHist.SetMarkerStyle(21)
      dataHist.SetMarkerSize(0.5)
      if(rebin):
         dataHist.Rebin(2)
         #dataHist.Rebin(4)
       
      dataHist.SetStats(0)
      stackplot.Draw("hist")
      print "dataHist.Integral(): ", dataHist.Integral()
      dataHist.Draw("psame") 
      if(log):
         pad1.SetLogy()
         ROOT.gPad.SetLogy()
      maximum=max(stackplot.GetMaximum(),dataHist.GetMaximum())
      minimum=min(stackplot.GetMinimum(),dataHist.GetMinimum(0))
      if(minimum==0):
         minimum=dataHist.GetMinimum(0)
      if not log:
         minimum=0      
      
      if not log:
         stackplot.SetMaximum(maximum*1.2)
      else:
         stackplot.SetMaximum(maximum*6)
      stackplot.SetMinimum(minimum*0.2)
      stackplot.GetXaxis().SetTitleOffset(2)
      stackplot.GetXaxis().SetTitle("")
      stackplot.GetXaxis().SetTickSize(0)
      stackplot.Draw("hist")
      dataHist.Draw("psame") 
  
      newDataHist = copy.deepcopy(dataHist)
      allMC = stackplot.GetStack().Last().Clone()

      herr = allMC.Clone('herr')
      theErrorGraph = ROOT.TGraphErrors(herr)
      theErrorGraph.SetFillColor(ROOT.kGray+2)
      theErrorGraph.SetFillStyle(3002)
      herr.SetFillColor(ROOT.kGray+2)
      herr.SetFillStyle(3002)
      herr.SetMarkerColor(1111)

      theErrorGraph.Draw('SAME2')

      leg.AddEntry(dataHist, "Data", "p") 
      leg.Draw()
      leg2.Draw()
      
      label = ROOT.TPaveText(0.4,0.905,0.99,0.96,"brNDC")
      label.SetTextSize(14)
      label.AddText(tag)
      label.SetFillColor(0)
      label.SetBorderSize(0)
#      label.Draw()
      
      stackplot.GetYaxis().SetTitle("Events")
      stackplot.GetYaxis().CenterTitle()
      stackplot.GetYaxis().SetTitleSize(20)
      stackplot.GetYaxis().SetTitleFont(43)
      stackplot.GetYaxis().SetTitleOffset(1.35)
      stackplot.GetYaxis().SetLabelFont(43)
      stackplot.GetYaxis().SetLabelSize(15)
#      stackplot.GetXaxis().SetLabelFont(43)
#      stackplot.GetXaxis().SetLabelSize(15)
#      stackplot.GetXaxis().SetTickSize(0.07)
      stackplot.GetXaxis().SetTitleOffset(2)
      stackplot.GetXaxis().SetTitle("")
      stackplot.GetXaxis().SetTickSize(0)
#      xlabel = ROOT.TPaveText(0.2,0.0,0.99,0.06,"brNDC")
#      xlabel.SetTextSize(20)
#      xlabel.AddText(xtitle)
#      xlabel.SetFillColor(0)
#      xlabel.SetBorderSize(0)
#      xlabel.Draw()
      
      lab_x0 = 0.10
      lab_y0 = 0.905
      tag1 = ROOT.TLatex(lab_x0,lab_y0,"CMS")
      tag1.SetNDC()
      tag1.SetTextFont(62)
      tag2 = ROOT.TLatex(lab_x0+0.06, lab_y0, "Internal")
      tag2.SetNDC()
      tag2.SetTextFont(52)
      tag1.SetTextSize(0.04)
      tag2.SetTextSize(0.03)
      tag1.Draw()
      tag2.Draw()

      plot.Modified()

      pad2.cd()

      ratio = self.getRatio(newDataHist,allMC)
      herr3= self.TOTerror(allMC,ratio);
      toterree = ROOT.TGraphErrors(herr3)

      ratio.SetStats(0)
      ratio.GetYaxis().SetRangeUser(0,3)
      ratio.GetYaxis().SetNdivisions(504)
      ratio.GetYaxis().SetTitle("Data/Simulation")
      ratio.GetXaxis().SetTitle(xtitle)

      ratio.GetXaxis().SetTitleSize(0.14)
      ratio.GetXaxis().SetTitleOffset(1.0)
      ratio.GetXaxis().SetLabelOffset(0.005)
      ratio.GetYaxis().SetTitleOffset(0.5)
      ratio.GetYaxis().SetLabelSize(0.12)
      ratio.GetYaxis().SetTitleSize(0.11)
      ratio.GetXaxis().SetLabelSize(0.11)

      ratio.SetLineStyle(1)
      ratio.SetLineWidth(2)

      line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0,
                      ratio.GetXaxis().GetXmax(), 1.0)
      line.SetLineColor(ROOT.kGray)
      line.SetLineStyle(2)
      line.Draw()
      tKsChi = ROOT.TLatex()
      tKsChi.SetNDC()
      tKsChi.SetTextFont(42)
      tKsChi.SetTextSize(0.09)

      ratio.Draw("PE")

      toterree.SetFillColor(ROOT.kGray+2)
      toterree.SetLineColor(ROOT.kGray+2)
      toterree.SetFillStyle(3002)
#        toterree.Draw("2 same")
      line.Draw("same")

      leg4 = ROOT.TLegend(0.7,0.89,0.5,0.8)
      leg4.SetFillStyle(0)
      leg4.SetBorderSize(0)
      leg4.SetTextSize(0.05)
      leg4.SetTextFont(42)
      leg4.AddEntry(toterree,"MC uncert. (stat.)","fl")
#      leg4.Draw()

      
      img = ROOT.TImage.Create()
      img.FromPad(plot)
      img.WriteImage(outDir+"/"+name.split("/")[0]+name.split("/")[1]+"_stack.png")
      plot.Print(outDir+"/"+name.split("/")[0]+name.split("/")[1]+"_stack.pdf")
      plot.Print(outDir+"/"+name.split("/")[0]+name.split("/")[1]+"_stack.C")
