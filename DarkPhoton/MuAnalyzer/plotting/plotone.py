import ROOT

def plotone(hist, xtitle, ytitle, title, log, ylow, yup,imgname,xlow=None,xup=None):
   plot = ROOT.TCanvas(imgname+"_plot",imgname+"_plot",2)
   x_low = 0
   y_low = 0
   x_up = 1.0
   y_up = 1.0
   pad1 = ROOT.TPad("pad1","pad1",x_low,y_low,x_up,y_up)
   pad1.SetBottomMargin(0.1)
   pad1.SetRightMargin(0.03)
   pad1.Draw()
   pad1.cd()
   hist.SetTitle("")
   #hist.SetStats(0)
   hist.SetMaximum(yup)
   hist.SetMinimum(ylow)
   if(log):
      pad1.SetLogy()
   hist.SetLineColor(1)
   hist.SetLineWidth(2)
   hist.Draw()
   label = ROOT.TPaveText(0.5,0.9,0.99,0.96,"brNDC")
   label.SetTextSize(20)
   label.AddText(title)
   label.SetFillColor(0)
   label.SetBorderSize(0)
   label.Draw()
   hist.GetYaxis().SetTitle(ytitle)
   hist.GetYaxis().CenterTitle()
   hist.GetYaxis().SetTitleSize(14)
   hist.GetYaxis().SetTitleFont(43)
   hist.GetYaxis().SetTitleOffset(1.55)
   hist.GetYaxis().SetLabelFont(43)
   hist.GetYaxis().SetLabelSize(15)
   hist.GetXaxis().SetLabelFont(43)
   hist.GetXaxis().SetLabelSize(15)
   hist.GetXaxis().SetTickSize(0.07)
   if xlow is not None:
      hist.GetXaxis().SetRangeUser(xlow,xup)
   xlabel = ROOT.TPaveText(0.5,0.0,0.99,0.06,"brNDC")
   xlabel.SetTextSize(20)
   xlabel.AddText(xtitle)
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
   
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   img.WriteImage(imgname+".png")
   
