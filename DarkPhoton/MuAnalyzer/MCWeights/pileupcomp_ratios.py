import ROOT
import tdrstyle 

dyf = ROOT.TFile("DYJets_Private_MCWeights.root","READ")
nomData = ROOT.TFile("/hdfs/cms/user/revering/dphoton/PileupHists/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root","READ")
downData = ROOT.TFile("/hdfs/cms/user/revering/dphoton/PileupHists/PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root","READ")
upData = ROOT.TFile("/hdfs/cms/user/revering/dphoton/PileupHists/PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root","READ")

nomHist = nomData.Get("pileup")
downHist = downData.Get("pileup")
upHist = upData.Get("pileup")
dyHist = dyf.Get("MCWeights/TrueNumInteractions")

nomHist.Scale(1/nomHist.Integral())
dyHist.Scale(1/dyHist.Integral())
downHist.Scale(1/downHist.Integral())
upHist.Scale(1/upHist.Integral())

nomHist.Divide(dyHist)
upHist.Divide(dyHist)
downHist.Divide(dyHist)

#for i in range(0,nomHist.GenNbinsX()):
#   if dyHist.GetBinContent(i)!=0:
#      dyHist.SetBinContent(i,1/dyHist.GetBinContent(i))
#   if wjhist.GetBinContent(i)!=0:
#      wjhist.SetBinContent(i,1/wjhist.GetBinContent(i))
#   if tthist.GetBinContent(i)!=0:
#      tthist.SetBinContent(i,1/tthist.GetBinContent(i))


plot = ROOT.TCanvas("puDists","puDists",2) 
x_low = 0
y_low = 0
x_up = 1.0
y_up = 1.0
pad1 = ROOT.TPad("pad1","pad1",x_low,y_low,x_up,y_up)
pad1.SetBottomMargin(0.1)
pad1.SetLeftMargin(0.15)
pad1.SetRightMargin(0.03)
pad1.Draw()
pad1.cd()
tdrstyle.setTDRStyle()
nomHist.SetTitle("")
nomHist.SetStats(0)
nomHist.SetMaximum(5)
nomHist.SetMinimum(0)
upHist.SetLineColor(1)
nomHist.SetLineColor(1)
downHist.SetLineColor(1)

upHist.SetLineStyle(2)
downHist.SetLineStyle(3)
nomHist.Draw("Hist")
upHist.Draw("sameHist")
downHist.Draw("sameHist")
nomHist.GetYaxis().SetTitle("Pileup Reweighting Factor")
nomHist.GetYaxis().SetTitleSize(23)
nomHist.GetYaxis().SetTitleFont(43)
nomHist.GetYaxis().SetTitleOffset(1.5)
nomHist.GetYaxis().SetLabelFont(43)
nomHist.GetYaxis().SetLabelSize(23)
nomHist.GetXaxis().SetTitleFont(43)
nomHist.GetXaxis().SetTitle("Number of Pileup")
nomHist.GetXaxis().SetTitleSize(23)
nomHist.GetXaxis().SetTitleOffset(1.0)
nomHist.GetXaxis().SetLabelFont(43)
nomHist.GetXaxis().SetLabelSize(23)
nomHist.SetLineWidth(2)
nomHist.SetLineWidth(2)
upHist.SetLineWidth(2)
downHist.SetLineWidth(2)

leg_x0 = 0.2*x_up
leg_y0 = 0.6*y_up
leg = ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.35,leg_y0+0.25)
leg.SetBorderSize(0)
leg.AddEntry(nomHist, "Data","l")
leg.AddEntry(upHist, "Data +1 s.d.","l")
leg.AddEntry(downHist, "Data -1 s.d.","l")
leg.Draw()
ROOT.gStyle.SetTextFont(43)
label = ROOT.TPaveText(0.25,0.9,0.99,0.96,"brNDC")
label.SetTextSize(25)
label.AddText("Pileup Distributions")
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
tag4 = ROOT.TLatex(0.7,0.91,"59.8 fb^{-1} (13 TeV)") 
tag4.SetNDC()  
tag4.SetTextFont(62)  
tag4.SetTextSize(0.036)    
tag4.Draw()

plot.cd()

plot.cd()
xlabel = ROOT.TPaveText(0.4,0.0,0.99,0.05,"brNDC")
xlabel.SetTextSize(26)
xlabel.AddText("Number of Pileup")
xlabel.SetFillColor(0)
xlabel.SetBorderSize(0)
xlabel.Draw()
img = ROOT.TImage.Create()
img.FromPad(plot)
img.WriteImage("PileupRatioDists.png")
plot.Print("PileupRatioDists.C")
plot.Print("PileupRatioDists.pdf")
