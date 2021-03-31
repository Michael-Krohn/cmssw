 #!/usr/bin/env python

import argparse
import os
import ROOT
import stackplot_tool as stack
import math

parser = argparse.ArgumentParser(description="Create stack plots for the mu+X skim")
parser.add_argument("resultsDir", help="directory where results histograms are contained")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd()+"/plots/")
parser.add_argument("--l", "--lumi", dest="lumi", help="luminosity", default=0)
parser.add_argument("--t", "--title", dest="label", help="Data label", default="Mu+X Data Skim")
parser.add_argument("--mc", "--mc", dest="mc", help="Fit MC", default=False)
parser.add_argument("--func", "--function", dest="function", help="Function to fit", default="landau")
arg = parser.parse_args()

resultsDir = arg.resultsDir
#Check for trailing slash on results dir and delete
if arg.resultsDir.split("/")[-1] =="": resultsDir = arg.resultsDir[:-1]
if not os.path.isdir(resultsDir):
     print "Results directory " + resultsDir + " does not exist!"

outDir = arg.outDir
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
if not os.path.isdir(outDir):
     print "Output directory " + outDir + " does not exist!"
     quit()

#Cross Sections in pb
#WJetsCx=61527.6*0.67
#TTToSemiLeptonicCx=365.34*3.915
#TTTo2L2NuCx=88.29*3.915
#DYJetsCx=6077.22*0.858
WJetsCx=61527.6
TTToSemiLeptonicCx=365.34
TTTo2L2NuCx=88.29
DYJetsCx=6077.22
DYJets10to50Cx=18610
WWCx=118.7

mcFiles = []
#mcFIle object goes Skim File,Results File, Label, Cross Section, Fill Color
#Histograms are stacked in the order they are added here
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/WW/hists/WW_monitor.root","READ"),ROOT.TFile(resultsDir+"/WW/hists/WW_analyzedsplit.root","READ"), "WW",WWCx,7))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_monitor.root","READ"),ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_analyzedsplit.root","READ"), "TTBar2L2Nu",TTTo2L2NuCx,2))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_monitor.root","READ"),ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_analyzedsplit.root","READ"), "TTBarSemiLeptonic",TTToSemiLeptonicCx,5))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/WJets/hists/WJets_monitor.root","READ"), ROOT.TFile(resultsDir+"/WJets/hists/WJets_analyzedsplit.root","READ"), "WJets", WJetsCx,3))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/DYJets/hists/DYJets_monitor.root","READ"),ROOT.TFile(resultsDir+"/DYJets/hists/DYJets_analyzedsplit.root","READ"), "DYJets",DYJetsCx,4))
mcFiles.append(stack.mcFile(ROOT.TFile(resultsDir+"/DYJets_10to50/hists/DYJets10to50_monitor.root","READ"),ROOT.TFile(resultsDir+"/DYJets_10to50/hists/DYJets10to50_analyzedsplit.root","READ"), "DYJets_10to50",DYJets10to50Cx,6))

dataFile = ROOT.TFile(resultsDir+"/Data/hists/2018_analyzedsplit.root","READ")

if not arg.mc:
   data = dataFile.Get("demo/nonMuonProbe/AdjustedMuonTrackMass")
mc = ROOT.TObjArray(3)
wJets = mcFiles[3].result.Get("demo/allEvents/AdjustedMuonTrackMass")
wJetsNevents=WJetsCx*1000*float(arg.lumi)*mcFiles[3].result.Get("demo/allEvents/AdjustedMuonTrackMass").Integral()/mcFiles[3].skim.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
wJets.Scale(wJetsNevents/wJets.Integral())
dyJets = mcFiles[4].result.Get("demo/allEvents/AdjustedMuonTrackMass")
dyJetsNevents=DYJetsCx*1000*float(arg.lumi)*mcFiles[4].result.Get("demo/allEvents/AdjustedMuonTrackMass").Integral()/mcFiles[4].skim.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
dyJets.Scale(dyJetsNevents/dyJets.Integral())
ttbar = mcFiles[1].result.Get("demo/allEvents/AdjustedMuonTrackMass")
ttbarNevents=TTTo2L2NuCx*1000*float(arg.lumi)*mcFiles[1].result.Get("demo/allEvents/AdjustedMuonTrackMass").Integral()/mcFiles[1].skim.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
ttbar.Scale(ttbarNevents/ttbar.Integral())
ttbarsemi = mcFiles[2].result.Get("demo/allEvents/AdjustedMuonTrackMass")
ttbarseminEvents=TTToSemiLeptonicCx*1000*float(arg.lumi)*mcFiles[2].result.Get("demo/allEvents/AdjustedMuonTrackMass").Integral()/mcFiles[2].skim.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
ttbarsemi.Scale(ttbarseminEvents/ttbarsemi.Integral())
ttbarNevents=ttbarNevents+ttbarseminEvents
ttbar.Add(ttbarsemi)
mc.Add(ttbar)
mc.Add(dyJets)
mc.Add(wJets)

mcStack = ROOT.THStack("stack","")
hists = []
for mcfile in mcFiles:
   eff = mcfile.result.Get("demo/nonMuonProbe/AdjustedMuonTrackMass").Integral()/mcfile.skim.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
   Nevents =  mcfile.cx*1000*float(arg.lumi)*eff
   hist=(mcfile.result.Get("demo/nonMuonProbe/AdjustedMuonTrackMass").Clone())
   if(hist.Integral()>0):
      hist.Scale(Nevents/hist.Integral())
   hist.SetFillColor(mcfile.color)
   hist.SetLineWidth(0)
   hist.SetDirectory(0)
   hists.append(hist)
   
for hist in hists:
   mcStack.Add(hist)
   if arg.mc:
      if hist == hists[0]:
         data = hist.Clone()
      else:
         data.Add(hist)
function=arg.function
#expo,poly,erf,landau, or gauss

#func = ROOT.TF1('func', ROOT.TF1Convolution("expo","TMath::Erf(x)",50,150,True),50,150,f_conv.GetNpar())
#if function=="expo":
#   opt_offset = 7000
#   opt_expo = -0.044
#   for i in range(0,50):
#      func = ROOT.TF1('func', str(opt_offset)+"+[0]*(x-[1])*e**([2]*(x))",50,150)
#      r1 = data.Fit('func','RS')
#      opt_expo = r1.Parameter(2)
#      func2 = ROOT.TF1('func2', "[2]+[0]*(x-[1])*e**("+str(opt_expo)+"*x)",50,150)
#      r2 = data.Fit('func2','RS')
#      opt_offset=r2.Parameter(2)
#      print("Iteration one, opt expo is "+str(opt_expo)+", opt offset is "+str(opt_offset)+".\n")
if function=="erf":
   opt_offset = 7000
   opt_expo = -0.044
   for i in range(0,50):
      func = ROOT.TF1('func', str(opt_offset)+"+[0]*(erf([2]*x))*e**([1]*(x))",50,150)
      r1 = data.Fit('func','S')
      opt_expo = r1.Parameter(1)
      func2 = ROOT.TF1('func2', "[1]+[0]*(erf([2]*x))*e**("+str(opt_expo)+"*x)",50,150)
      func2.SetParameter(2,40)
      r2 = data.Fit('func2','S')
      opt_offset=r2.Parameter(1)
      print("Iteration "+str(i)+", opt expo is "+str(opt_expo)+", opt offset is "+str(opt_offset)+".\n")

func = ROOT.TF1()

label = ""
if function=="expo":
   func = ROOT.TF1('func', "[3]+[0]*(x-[1])*e**([2]*(x))",50,150)
   func.SetParameters(12000,37,-0.03,5000)
   label = "f(x) = A + B(x-C)e^{-Dx}"
if function=="erf":
   func = ROOT.TF1('func', str(opt_offset)+"+[0]*(erf([1]*x))*e**([2]*(x))",50,150)
   label = "f(x) = A + B(erf(x))e^{-Cx}"
if function=="poly":
   func = ROOT.TF1('func', "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4",50,150)
   label = "f(x) = A + Bx + Cx^{2} + Dx^{3}" 
if function=="landau":
   func = ROOT.TF1('func', "landau(0)+[3]",50,150)
   func.SetParameters(27000,74,2,0)
   label = "Landau Distribution"
#   for i in range(87,94):
#      print("Bin number is: "+str(data.FindBin(i)))
#      data.SetBinError(data.FindBin(i),0)

if function=="gauss":
   func = ROOT.TF1('func', "gaus",50,150)
   label = "Gaussian Distribution"
   #for i in range(87,94):
      #data.SetBinError(i,0)

r1=data.Fit('func','RS')
finalfunc = data.GetFunction('func').Clone()
chi2  = finalfunc.GetChisquare()
canv = ROOT.TCanvas("fitplot","fitplot",2)
ROOT.gStyle.SetTextFont(43)
div_line = 0.35
pad1 = ROOT.TPad("pad1","pad1",0,div_line,1,1)
pad1.SetTopMargin(0.05)
pad1.SetBottomMargin(0.)
pad1.SetRightMargin(0.03)
#pad1.SetLeftMargin(0.13)
pad1.Draw()
pad1.cd()
data.SetStats(0)
data.SetLineWidth(0)
data.SetMarkerStyle(21)
data.SetMarkerSize(0.5)

ROOT.TGaxis().SetMaxDigits(3)
if not arg.mc:
   data.Draw("p")
   data.GetXaxis().SetTitle("Invariant Mass of Tag and Probe (GeV)")
   data.GetXaxis().SetTitleSize(0)
   data.GetXaxis().SetLabelFont(43)
   data.GetXaxis().SetRangeUser(50,150)
   data.GetXaxis().SetLabelSize(15)
   data.GetXaxis().SetTickSize(0.05)
   data.GetYaxis().SetTitle("Events at "+arg.lumi+" fb^{-1}")
   data.GetYaxis().SetTitleSize(0.05)
   data.GetYaxis().SetLabelFont(43)
   data.GetYaxis().SetLabelSize(15)
   data.SetMinimum(10.01)
else:
   mcStack.Draw("hist")
   finalfunc.Draw("same")
   mcStack.GetXaxis().SetTitle("Invariant Mass of Tag and Probe (GeV)")
   mcStack.GetXaxis().SetTitleSize(0)
   mcStack.GetXaxis().SetLabelFont(43)
   mcStack.GetXaxis().SetLabelSize(15)
   mcStack.GetXaxis().SetTickSize(0.05)
   mcStack.GetXaxis().SetRangeUser(50,150)
   mcStack.GetYaxis().SetTitle("Events at "+arg.lumi+" fb^{-1}")
   mcStack.GetYaxis().SetTitleSize(0.05)
   mcStack.GetYaxis().SetLabelFont(43)
   mcStack.GetYaxis().SetLabelSize(15)
   mcStack.SetMinimum(10.01)

plottitle = ROOT.TPaveText(0.4,0.955,0.99,0.99,"brNDC")
plottitle.SetTextSize(14)
if arg.mc:
   plottitle.AddText("MC, Probe Not Matched to Muon")
else:
   plottitle.AddText("Run A+D Data, Probe Not Matched to Muon")
plottitle.SetFillColor(0)
plottitle.SetBorderSize(0)
plottitle.Draw()

lab_x0 = 0.2
lab_y0 = 0.96
tag1 = ROOT.TLatex(lab_x0,lab_y0,"CMS")
tag1.SetTextFont(62)
tag1.SetNDC()
tag2 = ROOT.TLatex(lab_x0+0.075,lab_y0,"Internal")
tag2.SetTextFont(52)
tag2.SetNDC()
tag1.SetTextSize(0.05)
tag2.SetTextSize(0.04)
tag1.Draw()
tag2.Draw()

tag3 = ROOT.TLatex(0.6,0.6,label)
tag3.SetTextFont(52)
tag3.SetTextSize(0.04)
tag3.SetNDC()
tag3.Draw()
chitag = ROOT.TLatex(0.6,0.55,"Reduced #chi^{2}= "+str(round(chi2/finalfunc.GetNDF(),2)))
chitag.SetTextFont(52)
chitag.SetTextSize(0.04)
chitag.SetNDC()
chitag.Draw()
canv.cd()
pad2 = ROOT.TPad("pad2","pad2",0.0,0.06,1.0,div_line)
pad2.SetTopMargin(0.0)
pad2.SetBottomMargin(0.1)
pad2.SetRightMargin(0.03)
pad2.Draw()
pad2.cd()
ratio = data.Clone("h3")
ratio.SetStats(0)
if not arg.mc:
   ratio.Divide(finalfunc)
if arg.mc:
   for i in range(0,ratio.GetNbinsX()):
      if(ratio.GetBinContent(i)>0):
         ratio.SetBinContent(i,(ratio.GetBinContent(i)-finalfunc.Eval(ratio.GetBinLowEdge(i)+0.5))/math.sqrt(finalfunc.Eval(ratio.GetBinLowEdge(i)+0.5)))
   totalpull=0
   for i in range(0,ratio.GetNbinsX()):
      totalpull=totalpull+ratio.GetBinContent(i)
      ratio.SetBinContent(i,totalpull)
ratio.SetLineColor(1)
ratio.GetFunction("func").Delete()
ratio.Draw("p")
ratio.GetXaxis().SetRangeUser(50,150)
midline = ROOT.TLine(50,0,150,0)
midline.SetLineStyle(7)
midline.SetLineWidth(2)
midline.Draw()
ratio.SetTitle("")
if arg.mc:
   ratio.GetYaxis().SetTitle("Integrated Obs-Ex/sqrt(Ex)")
else:
   ratio.GetYaxis().SetTitle("Data/Fit")
ratio.SetMinimum(-99.9)
ratio.SetMaximum(99.9)
ratio.GetYaxis().CenterTitle()
ratio.GetYaxis().SetTitleSize(15)
ratio.GetYaxis().SetTitleFont(43)
ratio.GetYaxis().SetTitleOffset(1.72)
ratio.GetYaxis().SetLabelFont(43)
ratio.GetYaxis().SetLabelSize(15)
ratio.GetXaxis().SetLabelFont(43)
ratio.GetXaxis().SetTitleSize(0)
ratio.GetXaxis().SetLabelSize(15)
ratio.GetXaxis().SetTickSize(0.07)
canv.cd()
xlabel = ROOT.TPaveText(0.3,0.0,0.99,0.06,"brNDC")
xlabel.SetTextSize(20)
xlabel.AddText("Invariant Mass of Tag and Probe (GeV)")
xlabel.SetFillColor(0)
xlabel.SetBorderSize(0)
xlabel.Draw()

if arg.mc:
   canv.Print(function+"_mcfitresult.pdf")
else:
   canv.Print(function+"_fitresult.pdf")

dataFile.Close()
