import ROOT
import tdrstyle

def plotcomp(TH1, TH2, xtitle, ytitle, title, imagetitle, ratio,log, y_min, y_max,x_min=None,x_max=None):
   plot = ROOT.TCanvas(title,title,2) 
   x_low = 0
   y_low = 0
   x_up = 1.0
   y_up = 1.0
   if ratio:
      div_line = 0.45*(y_up-y_low)
      pad1 = ROOT.TPad("pad1","pad1",x_low,div_line,x_up,y_up)
   else:
      pad1 = ROOT.TPad("pad1","pad1",x_low,y_low,x_up,y_up)
   if ratio:
      pad1.SetBottomMargin(0.)
   else:
      pad1.SetBottomMargin(0.1)
      pad1.SetLeftMargin(0.15)
   pad1.SetRightMargin(0.03)
   pad1.Draw()
   pad1.cd()
   #tdrstyle.setTDRStyle()
   TH1.SetTitle("")
   TH1.SetStats(0)
   TH1.SetMaximum(y_max)
   TH1.SetMinimum(y_min)
   TH1.SetLineColor(2)
   #TH1.Scale(1/TH1.GetEntries())
   TH1.Draw("HistE")
   #TH2.Scale(1/TH2.GetEntries())
   TH2.Draw("sameHistE")
   TH1.GetYaxis().SetTitle(ytitle)
   TH1.GetYaxis().SetTitleSize(23)
   TH1.GetYaxis().SetTitleFont(43)
   TH1.GetYaxis().SetTitleOffset(1.5)
   TH1.GetYaxis().SetLabelFont(43)
   TH1.GetYaxis().SetLabelSize(23)
   TH1.GetXaxis().SetTitleFont(43)
   TH1.GetXaxis().SetTitle(xtitle)
   TH1.GetXaxis().SetTitleSize(23)
   TH1.GetXaxis().SetTitleOffset(1.0)
   TH1.GetXaxis().SetLabelFont(43)
   TH1.GetXaxis().SetLabelSize(23)
   if x_min is not None:
      TH1.GetXaxis().SetRangeUser(x_min,x_max)
   TH1.SetLineWidth(2)
   TH2.SetLineWidth(2)
   if log:
      pad1.SetLogy()
   leg_x0 = 0.2*x_up
   leg_y0 = 0.6*y_up
   leg = ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.45,leg_y0+0.25)
   leg.SetBorderSize(0)
   leg.AddEntry(TH1, "0.1 GeV Threshold","l")
   leg.AddEntry(TH2, "0.2 GeV Threshold","l")
   leg.Draw()
   ROOT.gStyle.SetTextFont(43)
   label = ROOT.TPaveText(0.25,0.9,0.99,0.96,"brNDC")
   label.SetTextSize(25)
   label.AddText(title)
   label.SetFillColor(0)
   label.SetBorderSize(0)
   #label.Draw() 
   tag2 = ROOT.TLatex(0.16,0.91,"CMS") 
   tag2.SetNDC()   
   tag2.SetTextFont(62)  
   tag3 = ROOT.TLatex(0.28,0.91,"Preliminary") 
   tag3.SetNDC()  
   tag3.SetTextFont(52)  
   tag2.SetTextSize(0.055) 
   tag3.SetTextSize(0.045)    
   tag2.Draw()
   tag3.Draw()
   plot.cd()
   if ratio:
      pad2 = ROOT.TPad("pad2","pad2",x_low, y_low+0.06, x_up, div_line)
      pad2.SetTopMargin(0.0)
      pad2.SetBottomMargin(0.1)
      pad2.SetRightMargin(0.03)
      pad2.Draw()
      pad2.cd()
      ratio = TH2.Clone("h3")
      nbins = ratio.GetNbinsX()
      for i in range(nbins):
#      Nentries1 = 0
#      Nentries2 = 0
#      for j in range(nbins-i):
#         Nentries1 = Nentries1 + TH1.GetBinContent(i+j)
	 #Nentries2 = Nentries2 + TH2.GetBinContent(i+j+1)
         if (TH1.Integral(i+1,nbins)) > 0:
            ratio.SetBinContent(i+1,TH2.Integral(i+1,nbins)/(TH1.Integral(i+1,nbins)))
      ratio.SetMinimum(0)
      ratio.SetMaximum(0.43)
      ratio.SetStats(0)
#     pad2.SetLogy()
#     ratio.SetMarkerStyle(2)
      ratio.SetLineColor(1)
      ratio.Draw()
   

      ratio.SetLineWidth(2)
      ratio.SetTitle("")
      ratio.GetYaxis().SetTitle("Random/Muon")
      ratio.GetYaxis().CenterTitle()
      ratio.GetYaxis().SetTitleSize(14)
      ratio.GetYaxis().SetTitleFont(43)
      ratio.GetYaxis().SetTitleOffset(1.55)
      ratio.GetYaxis().SetLabelFont(43)
      ratio.GetYaxis().SetLabelSize(15)

      #ratio.GetXaxis().SetTitle(xtitle)
      #ratio.GetXaxis().SetTitleSize(20)
      #ratio.GetXaxis().SetTitleFont(43)
      #ratio.GetXaxis().SetTitleOffset(3.)
      ratio.GetXaxis().SetLabelFont(43)
      ratio.GetXaxis().SetLabelSize(15)
      ratio.GetXaxis().SetTickSize(0.07)
  
   plot.cd()
   xlabel = ROOT.TPaveText(0.4,0.0,0.99,0.07,"brNDC")
   xlabel.SetTextSize(35)
   xlabel.AddText(xtitle)
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   img.WriteImage(imagetitle)
   plot.Print("/local/cms/user/revering/dphoton/DBremAdd/CMSSW_10_2_6/src/DarkPhoton/MuAnalyzer/python/"+imagetitle+".C")
   plot.Print("/local/cms/user/revering/dphoton/DBremAdd/CMSSW_10_2_6/src/DarkPhoton/MuAnalyzer/python/"+imagetitle+".pdf")
   return
