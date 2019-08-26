import ROOT

def peakfind(TH1):
   return peak

def plotcomp(TH1, TH2, xtitle, ytitle, layer):
   plot = ROOT.TCanvas(layer,"Plot comparison",2) 
   x_low = 0
   y_low = 0
   x_up = 1.0
   y_up = 1.0
   div_line = 0
   pad1 = ROOT.TPad("pad1","pad1",x_low,div_line,x_up,y_up)
   pad1.SetBottomMargin(0.1)
   pad1.SetRightMargin(0.03)
   pad1.Draw()
   pad1.cd()
   TH1.SetTitle("")
   TH1.SetStats(0)
   TH1.SetMaximum(500000)
   TH1.SetMinimum(0.5)
   TH1.SetLineColor(2)
   TH1.Scale(1/TH1.GetEntries())
   TH1.Draw("Hist")
   TH2.Scale(1/TH2.GetEntries())
   TH2.Draw("histsame")
   TH1.GetYaxis().SetTitle(ytitle)
   TH1.GetYaxis().SetTitleSize(20)
   TH1.GetYaxis().SetTitleFont(43)
   TH1.GetYaxis().SetTitleOffset(1.05)
   TH1.GetYaxis().SetLabelFont(43)
   TH1.GetYaxis().SetLabelSize(15)
   pad1.SetLogy()
   leg_x0 = 0.6*x_up
   leg_y0 = 0.5*y_up
   leg = ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.35,leg_y0+0.15)
   leg.SetBorderSize(0)
   leg.AddEntry(TH1, "Depth Before Missing Hit","l")
   leg.AddEntry(TH2, "Depth Before Found Hit","l")
   leg.Draw()
   ROOT.gStyle.SetTextFont(43)
   label = ROOT.TPaveText(0.7,0.9,0.99,0.96,"brNDC")
   label.SetTextSize(20)
   label.AddText(layer)
   label.SetFillColor(0)
   label.SetBorderSize(0)
   label.Draw() 
   plot.cd()
   TH1.SetLineWidth(2)
   TH2.SetLineWidth(2)
  # ratio.SetTitle("")
  # ratio.GetYaxis().SetTitle("Random/Muon")
  # ratio.GetYaxis().CenterTitle()
  # ratio.GetYaxis().SetTitleSize(14)
  # ratio.GetYaxis().SetTitleFont(43)
  # ratio.GetYaxis().SetTitleOffset(1.55)
  # ratio.GetYaxis().SetLabelFont(43)
  # ratio.GetYaxis().SetLabelSize(15)

   #ratio.GetXaxis().SetLabelFont(43)
   #ratio.GetXaxis().SetLabelSize(15)
   #ratio.GetXaxis().SetTickSize(0.07)
  
   plot.cd()
   xlabel = ROOT.TPaveText(0.5,0.0,0.99,0.06,"brNDC")
   xlabel.SetTextSize(20)
   xlabel.AddText(xtitle)
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   img.WriteImage(layer.replace(" ","_")+".png")
   return
