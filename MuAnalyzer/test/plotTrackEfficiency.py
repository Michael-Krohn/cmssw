from optparse import OptionParser
import ROOT as rt
from array import *
import os
import random
import sys
import math
import tdrstyle

def setStyle():
    tdrstyle.setTDRStyle()
    rt.gStyle.SetPadTopMargin(0.10)
    rt.gStyle.SetPadLeftMargin(0.16)
    rt.gStyle.SetPadRightMargin(0.10)
    rt.gStyle.SetPaintTextFormat("1.3f")
    rt.gStyle.SetOptFit(0000)
    rt.gROOT.SetBatch()
#    rt.gStyle.SetOptStat(0)
#    rt.gStyle.SetOptFit(0000)
#    rt.gStyle.SetOptTitle(0)
#    rt.gStyle.SetPaintTextFormat("1.2g")
    #rt.gStyle.SetPalette()
#    rt.gStyle.SetNumberContours(999)
#    rt.gROOT.SetBatch()
#    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)
    
#    rt.gStyle.SetStatY(1.9)
#    rt.gStyle.SetStatX(1.9)
#    rt.gStyle.SetStatW(0.1)
#    rt.gStyle.SetStatH(0.1)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d','--dir',dest="outDir",default="plots",type="string",
                  help="Output directory to store results")
    parser.add_option('-l','--lumi',dest="lumi", default=300.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--numerator',dest="numerator",default="demo/histoDir/recoAK8pT_DoubleBTrig_np4",type="string", 
                  help="numerator trigger")
    parser.add_option('--denominator',dest="denominator",default="demo/histoDir/recoAK8pT",type="string",
                  help="denominator trigger")
    parser.add_option('-i',dest="inputFile",default="AnalysisOutput_DoubleB.root",type="string",
                  help="denominator trigger")


    
    (options,args) = parser.parse_args()
     

    tfile = rt.TFile.Open(options.inputFile,'read')

    num = tfile.Get(options.numerator)
    denom = tfile.Get(options.denominator)

    print "num.GetNbinsX(): ", num.GetNbinsX()
    print "denom.GetNbinsX(): ", denom.GetNbinsX()

    print "num.Integral(): ", num.Integral()
    print "denom.Integral(): ", denom.Integral()

    # rebin
    denom.SetDirectory(0)
    num.SetDirectory(0)
    denom.RebinY(10)
    num.RebinY(10)
    denom.RebinX(10)
    num.RebinX(10)
#    denom.RebinY(20)
#    num.RebinY(20)
#    denom.RebinX(25)
#    num.RebinX(25)

#    num.Rebin(2)
#    denom.Rebin(2)
     
    setStyle()
    c = rt.TCanvas("c_"+options.numerator.split('/')[-1],"c_"+options.numerator.split('/')[-1],500,400)
    c.SetRightMargin(0.15)

    pEff = rt.TEfficiency(num, denom)
    pEff.SetDirectory(0)
    pEff.SetName('eff_'+options.numerator.split('/')[-1])
    xaxisTitle = "#eta"
    yaxisTitle = "#phi"
    if 'recoAK8pT' in options.numerator:
	xaxisTitle = "Lead AK8 jet p_{T} (GeV)"
    elif 'recoAK8csv' in options.numerator:
        xaxisTitle = "CSV"
    elif 'recoAK8DoubleB' in options.numerator:
        xaxisTitle = "Double-b"
    elif 'recoAK8Maxcsv' in options.numerator:
        xaxisTitle = "CSV"
    elif 'recoAK8MaxDoubleB' in options.numerator:
        xaxisTitle = "Lead AK8 Double-b"
    elif 'recoAK8eta' in options.denominator:
        xaxisTitle = "Lead AK8 jet #eta"
    elif 'recoAK8phi' in options.denominator:
        xaxisTitle = "Lead AK8 jet #phi"
    
    pEff.SetTitle("Efficiency;%s;%s"%(xaxisTitle,yaxisTitle))
    
    pEff.SetMarkerSize(0.8)
    pEff.SetMarkerStyle(20)
#    pEff.Draw("apez")
    pEff.Draw("colztext")
    
    rt.gPad.Update()        
#    pEff.GetPaintedGraph().SetMarkerStyle(8)
#    pEff.GetPaintedGraph().SetMarkerSize(20)        
#    pEff.GetPaintedGraph().SetMinimum(0)
#    pEff.GetPaintedGraph().SetMaximum(1.1)
#    pEff.GetPaintedGraph().GetXaxis().SetTitle('Wide m_{jj}')
#    pEff.GetPaintedGraph().GetXaxis().SetRangeUser(250,900)
 
    lumi = 8.672 
   
    rt.gPad.Update()      
    tag1 = rt.TLatex(0.63,0.91,"%.1f fb^{-1} (13 TeV)"%lumi)
    tag1.SetNDC(); tag1.SetTextFont(42)
    tag1.SetTextSize(0.045)
    tag2 = rt.TLatex(0.16,0.91,"CMS")
    tag2.SetNDC()
    tag2.SetTextFont(62)
    tag3 = rt.TLatex(0.25,0.91,"Preliminary")
    tag3.SetNDC()
    tag3.SetTextFont(52)
    tag2.SetTextSize(0.055)
    tag3.SetTextSize(0.045)
#    tag1.Draw()
    tag2.Draw()
    tag3.Draw()

    c.Print(options.outDir+"/"+"HltEff_" + options.numerator.split('/')[-1] + ".png")
    c.Print(options.outDir+"/"+"HltEff_" + options.numerator.split('/')[-1] + ".pdf")
    c.Print(options.outDir+"/"+"HltEff_" + options.numerator.split('/')[-1] + ".C")


    rootFile = rt.TFile.Open('HltEff.root','update')
    tdirectory = rootFile.GetDirectory('HltEff')
    if tdirectory==None:
        rootFile.mkdir('HltEff')
        tdirectory = rootFile.GetDirectory('HltEff')
    
    tdirectory.cd()
    pEff.Write()
    c.Write()
