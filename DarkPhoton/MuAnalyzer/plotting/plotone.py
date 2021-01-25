import ROOT

def plotone(hist, xtitle, ytitle, title, log, ylow, yup,imgname,outDir="",xlow=None,xup=None):
   plot = ROOT.TCanvas(imgname+"_plot",imgname+"_plot",2)
   ROOT.gStyle.SetTextFont(43)
   x_low = 0
   y_low = 0
   x_up = 1.0
   y_up = 1.0
   pad1 = ROOT.TPad("pad1","pad1",x_low,y_low,x_up,y_up)
   pad1.SetTopMargin(0.05)
   pad1.SetBottomMargin(0.1)
   pad1.SetRightMargin(0.03)
   pad1.SetLeftMargin(0.13)
   pad1.Draw()
   pad1.cd()
   hist.SetTitle("")
   hist.SetStats(0)
   hist.Scale(1/hist.GetEntries())
   hist.SetMaximum(yup)
   hist.SetMinimum(ylow)
   if(log):
      pad1.SetLogy()
   hist.SetLineColor(1)
   hist.SetLineWidth(2)
   hist.Draw("Hist")
   label = ROOT.TPaveText(0.7,0.96,0.99,0.99,"brNDC")
   label.SetTextSize(20)
   label.AddText(title)
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
   if xlow is not None:
      hist.GetXaxis().SetRangeUser(xlow,xup)
   xlabel = ROOT.TPaveText(0.4,0.0,0.99,0.06,"brNDC")
   xlabel.SetTextSize(20)
   xlabel.AddText(xtitle)
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
 
   lab_x0 = 0.13
   lab_y0 = 0.96 
   tag1 = ROOT.TLatex(lab_x0,lab_y0,"CMS")
   tag1.SetNDC()
   tag1.SetTextFont(62)
   tag2 = ROOT.TLatex(lab_x0+0.085, lab_y0, "Preliminary")
   tag2.SetNDC()
   tag2.SetTextFont(52)
   tag1.SetTextSize(0.04)
   tag2.SetTextSize(0.03)
   tag1.Draw()
   tag2.Draw()
 
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   img.WriteImage(outDir+"/"+imgname+".png")
   plot.Print(outDir+"/"+imgname+".pdf")
   plot.Print(outDir+"/"+imgname+".C")
   
