import ROOT

def plotcomp(TH1, TH2, TH3, xtitle, ytitle, th_x, mA, ebeam):
   plot = ROOT.TCanvas("Plot comparison"+th_x,"Plot comparison",2) 
   div_line = 0.
   pad1 = ROOT.TPad("pad1","pad1",0,div_line,1.0,1.0)
   pad1.SetBottomMargin(0.1)
   pad1.SetRightMargin(0.03)
   pad1.Draw()
   pad1.cd()
   TH1.SetTitle("")
   TH1.SetStats(0)
   TH1.SetLineColor(2)
   TH1.Scale(1/TH1.GetEntries(),"width")
   TH1.Draw("Hist")
   TH2.Scale(1/TH2.GetEntries(),"width")
   TH2.Draw("sameHist")
   TH3.Scale(1/TH3.GetEntries(),"width")
   TH3.Draw("sameHist")
   TH3.SetLineColor(3)
   TH1.GetYaxis().SetTitle(ytitle)
   TH1.GetYaxis().SetTitleSize(17)
   TH1.GetYaxis().SetTitleFont(43)
   TH1.GetYaxis().SetTitleOffset(1.15)
   TH1.GetYaxis().SetLabelFont(43)
   TH1.GetYaxis().SetLabelSize(15)
   TH1.GetXaxis().SetTitle("#theta")
   TH1.GetXaxis().SetTitleSize(17)
   TH1.GetXaxis().SetTitleFont(43)

   TH1.SetMinimum(0.00001)
   if th_x == "x":
      TH1.SetMaximum(20.0)
   if th_x == "theta":
      TH1.SetMaximum(200.0)
   pad1.SetLogy()
   if th_x == "theta":
      leg_x0 = 0.5
      leg_y0 = 0.6
   if th_x == "x":
      leg_x0 = 0.3
      leg_y0 = 0.3
   leg = ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.25,leg_y0+0.15)
   leg.SetBorderSize(0)
   leg.AddEntry(TH1, TH1.GetName(),"l")
   leg.AddEntry(TH2, TH2.GetName(),"l")
   leg.AddEntry(TH3, TH3.GetName(),"l")
   leg.Draw()
   ROOT.gStyle.SetTextFont(43)
   label = ROOT.TPaveText(0.1,0.9,0.99,0.96,"brNDC")
   label.SetTextSize(17)
   label.AddText("Muon Recoil Angular Distribution, Beam energy = "+str(ebeam)+" GeV")
   label.SetFillColor(0)
   label.SetBorderSize(0)
   label.Draw() 
   plot.cd()
   TH1.SetLineWidth(2)
   TH2.SetLineWidth(2)
   TH3.SetLineWidth(2)
   #ratio.GetXaxis().SetTitle(xtitle)
   #ratio.GetXaxis().SetTitleSize(20)
   #ratio.GetXaxis().SetTitleFont(43)
   #ratio.GetXaxis().SetTitleOffset(3.)
   TH1.GetXaxis().SetLabelFont(43)
   TH1.GetXaxis().SetLabelSize(15)
   TH1.GetXaxis().SetTickSize(0.07)
  
   plot.cd()
   xlabel = ROOT.TPaveText(0.5,0.0,0.99,0.06,"brNDC")
   xlabel.SetTextSize(20)
   xlabel.AddText(xtitle)
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   if th_x == "theta":
      img.WriteImage("newersimtheta_map_"+str(mA)+"_ebeam_"+str(ebeam)+".png")
   if th_x == "x":
      img.WriteImage("newersimx_map_"+str(mA)+"_ebeam_"+str(ebeam)+".png")
   return plot
