 #!/usr/bin/env python

import ROOT

class histfile:
   def __init__(self, mcfile, name, fillcolor):
      self.mcfile=mcfile
      self.label=name
      self.color=fillcolor

class stackInfo:
   def __init__(self, mc):
      self.mcFiles=mc
   def plotAll(self,basedir,tag,outDir):
      #Histogram name, X-axis label, Log Scale, Title, Rebin
      self.stack(basedir+"/ConeHits","Hits along probe trajectory in HE",False,tag,outDir)
      self.stack(basedir+"/ConeEnergy","Energy of cells along probe trajectory in HE",True,tag,outDir)
      self.stack(basedir+"/dRtoNearestCSC","#Delta R from probe to nearest CSC hit",True,tag,outDir)
      self.stack(basedir+"/StandaloneMuondEoverE","(Standalone P - Probe P)/Probe P",True,tag,outDir)
      self.stack(basedir+"/DrtoNearestStandalone","#Delta R to Nearest Standalone Mion",True,tag,outDir)
      self.stack(basedir+"/signalCategory","Signal Category",True,tag,outDir)
      self.stack(basedir+"/CaloConeEnergy","Energy in HE Cone (GeV)",True,tag,outDir)
      self.stack(basedir+"/Dark Brem R","R of Dark Brem Location",True,tag,outDir)
      self.twoDshow(basedir+"/HEEStandDEOverE","Energy in HE Hits GeV","(Standalone P - Probe P)/Probe P)", False,tag,outDir,1)
      self.twoDshow(basedir+"/FailLocation","Probe i#eta","Probe i#phi", False, tag, outDir,1)
      self.stack("allEvents/cutProgress","Cut Progress",True,tag,outDir,False,True)
      self.stack("allEvents/ProbeEcalIsolation","Ecal Energy in 0.4 #Delta R  Cone", True, tag, outDir,False,False,35)
      self.stack("allEvents/TrackIsolation","Relative Track Isolation (Track pt / #Sigma pt in 0.3 #Delta R)", True, tag, outDir,False,False,0.3)
      self.stack("allEvents/MinCaloJetDr","#Delta R from Probe to Nearest Calo Jet", False, tag, outDir,False,False,3)
      self.stack("allEvents/ProbePFJetDr","#Delta R from Probe to Nearest PF Jet", True, tag, outDir)
      self.cutFlowText(outDir)

   def stack(self, name, xtitle, log, tag, outDir,rebin=False, scale=False, xmax=-1):
      leg = ROOT.TLegend(0.7,0.35,0.95,0.45)
      leg.SetBorderSize(0)
      hists = []
      NtotalEvents=0
      #nextname = name[:-1]+str(int(name[-1])+1)
      #print nextname
      for mcFile in self.mcFiles:
         hist=(mcFile.mcfile.Get("demo/"+name).Clone())
         if scale:
            hist.Scale(1./hist.GetBinContent(2))
         else:
            if(hist.Integral()>0):
               hist.Scale(1./hist.Integral())
         hist.SetLineWidth(2)
         hist.SetLineColor(mcFile.color)
         hist.SetDirectory(0)
         hist.SetStats(0)
         hists.append(hist)
         leg.AddEntry(hist, mcFile.label,"l")

      plot= ROOT.TCanvas("plot","plot",2)
      ROOT.gStyle.SetTextFont(43)
      pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
      pad1.SetTopMargin(0.05)
      pad1.SetBottomMargin(0.1)
      pad1.SetRightMargin(0.03)
      pad1.SetLeftMargin(0.13)
      pad1.Draw()
      pad1.cd()
  
      for hist in hists:
         hist.Draw("histsame")
      if(log):
         pad1.SetLogy()
         ROOT.gPad.SetLogy()
 
      leg.Draw()
      
      label = ROOT.TPaveText(0.4,0.96,0.99,0.99,"brNDC")
      label.SetTextSize(20)
      label.AddText(tag)
      label.SetFillColor(0)
      label.SetBorderSize(0)
      label.Draw()
      
      hists[0].GetYaxis().SetTitle("Fraction of events")
      hists[0].GetYaxis().CenterTitle()
      hists[0].GetYaxis().SetTitleSize(20)
      hists[0].GetYaxis().SetTitleFont(43)
      hists[0].GetYaxis().SetTitleOffset(1.35)
      hists[0].GetYaxis().SetLabelFont(43)
      hists[0].GetYaxis().SetLabelSize(15)
      hists[0].GetXaxis().SetLabelFont(43)
      hists[0].GetXaxis().SetLabelSize(15)
      hists[0].GetXaxis().SetTickSize(0.07)
      hists[0].SetMinimum(0.00001)
      if xmax>0:
         hists[0].GetXaxis().SetRangeUser(0,xmax)
      xlabel = ROOT.TPaveText(0.2,0.0,0.99,0.06,"brNDC")
      xlabel.SetTextSize(20)
      xlabel.AddText(xtitle)
      xlabel.SetFillColor(0)
      xlabel.SetBorderSize(0)
      xlabel.Draw()
      
      lab_x0 = 0.20
      lab_y0 = 0.96
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

   def twoDshow(self, name, xtitle, ytitle, log, tag, outDir, index):
      mcFile=self.mcFiles[index]
      hist=(mcFile.mcfile.Get("demo/"+name).Clone())
      if(hist.Integral()>0):
         hist.Scale(1./hist.Integral())
 
      hist.SetStats(0)
      plot= ROOT.TCanvas("plot","plot",2)
      ROOT.gStyle.SetTextFont(43)
      pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
      pad1.SetTopMargin(0.05)
      pad1.SetBottomMargin(0.1)
      pad1.SetRightMargin(0.03)
      pad1.SetLeftMargin(0.13)
      pad1.Draw()
      pad1.cd()
  
      hist.Draw("colz")
      if(log):
         pad1.SetLogy()
         ROOT.gPad.SetLogy()
 
      label = ROOT.TPaveText(0.4,0.96,0.99,0.99,"brNDC")
      label.SetTextSize(20)
      label.AddText(tag)
      label.SetFillColor(0)
      label.SetBorderSize(0)
      label.Draw()
      
      hist.GetYaxis().SetTitle(ytitle)
      hist.GetYaxis().CenterTitle()
      hist.GetYaxis().SetTitleSize(20)
      hist.GetYaxis().SetTitleFont(43)
      hist.GetYaxis().SetTitleOffset(1.35)
      hist.GetYaxis().SetLabelFont(43)
      hist.GetYaxis().SetLabelSize(15)
      hist.GetXaxis().SetLabelFont(43)
      hist.GetXaxis().SetLabelSize(15)
      hist.GetXaxis().SetTickSize(0.07)
      xlabel = ROOT.TPaveText(0.2,0.0,0.99,0.06,"brNDC")
      xlabel.SetTextSize(20)
      xlabel.AddText(xtitle)
      xlabel.SetFillColor(0)
      xlabel.SetBorderSize(0)
      xlabel.Draw()
      
      lab_x0 = 0.20
      lab_y0 = 0.96
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
 
   def cutFlowText(self, outDir):
      outfile = open(outDir+"/CutProgress.txt","w")
      cutList = [ 'AllEvents',
                  'Probe Track with >20 GeV pt',
                  'Probe Track is in endcap',
                  'Probe Track is high purity',
                  'Probe Track dxy<0.05',
                  'Probe Track Belongs to Primary Vertex',
                  'Tag Muon passes tight ID',
                  'Tag Muon passes tight PF Iso',
                  'Tag Muon pt >26 GeV and |eta|<2.4',
                  'Tag Muon is global muon',
                  'DR between tag and probe >0.2',
                  'Tag and Probe have opposite charge',
                  'Fitted tag/probe vertex is valid',
                  'Fitted tag/probe vertex has chi^2/ndof<3',
                  'Invariant mass of pair between 50 and 150 GeV',
                  'Probe track isolation <0.05',
                  'Probe ecal isolation <10 GeV',
                  'No other probes more than 0.3 DR from selected probe pair with tag',
                  'Closest HE cell to probe trajectory between 18 and 28',
                  '>2 hits along trajectory in HE >0.1 GeV',
                  '<10 GeV in HE hits along trajectory',
                  'No calo jet within DR of 1.0',
                  'DR to nearest CSC hit >0.05',
                  'DR to nearest Standalone muon >0.5']


      for mcFile in self.mcFiles:
         hist=(mcFile.mcfile.Get("demo/allEvents/cutProgress").Clone())
         outfile.write("\n")
         for n in range(1,len(cutList)+1):
             if(hist.GetBinContent(n)>0):
                 outfile.write(cutList[n-1]+","+str(hist.GetBinContent(n))+"\n")
      hist=self.mcFiles[2].mcfile.Get("demo/genAccepted/cutProgress").Clone()
      outfile.write("\n")
      for n in range(1,len(cutList)+1):
          if(hist.GetBinContent(n)>0):
              outfile.write(cutList[n-1]+","+str(hist.GetBinContent(n))+"\n")

      hist=self.mcFiles[0].mcfile.Get("demo/genAccepted/cutProgress").Clone()
      outfile.write("\n")
      for n in range(1,len(cutList)+1):
          if(hist.GetBinContent(n)>0):
              outfile.write(cutList[n-1]+","+str(hist.GetBinContent(n))+"\n")

      outfile.close()          

