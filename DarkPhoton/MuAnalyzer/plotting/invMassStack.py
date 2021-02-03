 #!/usr/bin/env python

import argparse
import os
import ROOT
import plotone
import plot2d

def calcScaleFactor(monitor, analyzed):
   #Filter efficiency first
   eff = analyzed.Get("demo/MuonsTrackMass").Integral()/monitor.Get("DiMuonFilter/allEvents/NumberOfMuonsPassingTag").Integral()
   return eff

parser = argparse.ArgumentParser(description="Create stack plots for the mu+X skim")
parser.add_argument("resultsDir", help="directory where results histograms are contained")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd()+"/plots/")
parser.add_argument("--l", "--lumi", dest="lumi", help="luminosity", default=0)
parser.add_argument("--t", "--title", dest="label", help="Data label", default="Mu+X Data Skim")
arg = parser.parse_args()


WJetsCx=61527.6
TTToSemiLeptonicCx=364.34
TTTo2L2NuCx=88.29
DYJetsCx=6077.22
DYJets10to50Cx=18610

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

dataAnalyzed = ROOT.TFile(resultsDir+"/Data/hists/MuPxSkim_analyzed.root","READ")
wJetsMonitor = ROOT.TFile(resultsDir+"/WJets/hists/WJets_monitor.root","READ")
wJetsAnalyzed = ROOT.TFile(resultsDir+"/WJets/hists/WJets_analyzed.root","READ")
TTBarMonitor = ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_monitor.root","READ")
TTBarAnalyzed = ROOT.TFile(resultsDir+"/TTToSemiLeptonic/hists/TTToSemiLeptonic_analyzed.root","READ")
TTBar2LMonitor = ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_monitor.root","READ")
TTBar2LAnalyzed = ROOT.TFile(resultsDir+"/TTTo2L2Nu/hists/TTTo2L2Nu_analyzed.root","READ")
DYJetsMonitor = ROOT.TFile(resultsDir+"/DYJets/hists/DYJets_monitor.root","READ")
DYJetsAnalyzed = ROOT.TFile(resultsDir+"/DYJets/hists/DYJets_analyzed.root","READ")
DYJets10to50Monitor = ROOT.TFile(resultsDir+"/DYJets_10to50/hists/DYJets10to50_monitor.root","READ")
DYJets10to50Analyzed = ROOT.TFile(resultsDir+"/DYJets_10to50/hists/DYJets10to50_analyzed.root","READ")

wJetsEff = calcScaleFactor(wJetsMonitor,wJetsAnalyzed)
TTBarEff = calcScaleFactor(TTBarMonitor,TTBarAnalyzed)
TTBar2LEff = calcScaleFactor(TTBar2LMonitor,TTBar2LAnalyzed)
DYJetsEff = calcScaleFactor(DYJetsMonitor,DYJetsAnalyzed)
DYJets10to50Eff = calcScaleFactor(DYJets10to50Monitor,DYJets10to50Analyzed)

wJetsNevents = WJetsCx*1000*float(arg.lumi)*wJetsEff
TTBarNevents = TTToSemiLeptonicCx*1000*float(arg.lumi)*TTBarEff
TTBar2LNevents = TTTo2L2NuCx*1000*float(arg.lumi)*TTBar2LEff
DYJetsNevents = DYJetsCx*1000*float(arg.lumi)*DYJetsEff
DYJets10to50Nevents = DYJets10to50Cx*1000*float(arg.lumi)*DYJets10to50Eff

stackplot = ROOT.THStack("stack","")

wJetsMass = wJetsAnalyzed.Get("demo/MuonsTrackMass")
wJetsMass.Scale(wJetsNevents/wJetsMass.Integral())
wJetsMass.SetFillColor(3)
wJetsMass.SetLineWidth(0)
DYJetsMass = DYJetsAnalyzed.Get("demo/MuonsTrackMass")
DYJetsMass.Scale(DYJetsNevents/DYJetsMass.Integral())
DYJetsMass.SetFillColor(4)
DYJetsMass.SetLineWidth(0)
DYJets10to50Mass = DYJets10to50Analyzed.Get("demo/MuonsTrackMass")
DYJets10to50Mass.Scale(DYJets10to50Nevents/DYJets10to50Mass.Integral())
DYJets10to50Mass.SetFillColor(6)
DYJets10to50Mass.SetLineWidth(0)
TTBarMass = TTBarAnalyzed.Get("demo/MuonsTrackMass")
TTBarMass.Scale(TTBarNevents/TTBarMass.Integral())
TTBarMass.SetFillColor(5)
TTBarMass.SetLineWidth(0)
TTBar2LMass = TTBar2LAnalyzed.Get("demo/MuonsTrackMass")
TTBar2LMass.Scale(TTBar2LNevents/TTBar2LMass.Integral())
TTBar2LMass.SetFillColor(2)
TTBar2LMass.SetLineWidth(0)

stackplot.Add(TTBar2LMass)
stackplot.Add(TTBarMass)
stackplot.Add(wJetsMass)
stackplot.Add(DYJetsMass)
stackplot.Add(DYJets10to50Mass)

plot= ROOT.TCanvas("stackplot","stackplot",2)
ROOT.gStyle.SetTextFont(43)
pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
pad1.SetTopMargin(0.05)
pad1.SetBottomMargin(0.1)
pad1.SetRightMargin(0.03)
pad1.SetLeftMargin(0.13)
pad1.Draw()
pad1.cd()

stackplot.Draw("hist")
dataMass = dataAnalyzed.Get("demo/MuonsTrackMass")
dataMass.SetLineWidth(0)
dataMass.SetLineColor(1)
dataMass.SetMarkerStyle(21)
dataMass.Draw("psame")

print("Data events: "+str(dataMass.Integral())+ "\nExpected Events: " + str(wJetsNevents+TTBarNevents+TTBar2LNevents+DYJetsNevents+DYJets10to50Nevents))

leg = ROOT.TLegend(0.6,0.6,0.85,0.85)
leg.SetBorderSize(0)
leg.AddEntry(wJetsMass, "WJets","f")
leg.AddEntry(TTBarMass, "TTBarSemiLeptonic","f")
leg.AddEntry(TTBar2LMass, "TTBar2L2Nu","f")
leg.AddEntry(DYJetsMass,"DYJets","f")
leg.AddEntry(DYJets10to50Mass,"DYJets_10to50","f")
leg.AddEntry(dataMass,"Data","p")
leg.Draw()

label = ROOT.TPaveText(0.4,0.96,0.99,0.99,"brNDC")
label.SetTextSize(20)
label.AddText("Non Isolated Events")
label.SetFillColor(0)
label.SetBorderSize(0)
label.Draw()

stackplot.GetYaxis().SetTitle("Events at "+arg.lumi+" fb^{-1}")
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
xlabel.AddText("Invariant Mass of Tag and Probe (GeV)")
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
img.WriteImage(outDir+"/stackplot.png")
plot.Print(outDir+"/stackplot.pdf")
plot.Print(outDir+"/stackplot.C")

dataAnalyzed.Close()
wJetsMonitor.Close()
wJetsAnalyzed.Close()
TTBarMonitor.Close()
TTBarAnalyzed.Close()
DYJetsMonitor.Close()
DYJetsAnalyzed.Close() 
