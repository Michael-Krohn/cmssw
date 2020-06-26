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
    rt.gStyle.SetPadTopMargin(0.07)
    rt.gStyle.SetPadLeftMargin(0.14)
    rt.gStyle.SetPadRightMargin(0.10)
    rt.gStyle.SetPadBottomMargin(0.11)
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
    parser.add_option('-i',dest="inputFile",default="DoubleMuon_data.root",type="string", help="Input root file to read from")
    
    (options,args) = parser.parse_args()

    tfile = rt.TFile.Open(options.inputFile,'read')

    DBLoc = tfile.Get("demo/DarkBremLocation")

    DBLoc.SetDirectory(0)
     
    setStyle()
    c = rt.TCanvas("DarkBremLocation","DarkBremLocation",400,400)
    c.SetRightMargin(0.15)

    xaxisTitle = "Z (cm)"
    yaxisTitle = "R (cm)"
    
    DBLoc.SetTitle("Dark Brem Location;%s;%s"%(xaxisTitle,yaxisTitle))
    
    DBLoc.SetLineColor(2)
    DBLoc.SetLineWidth(2)
    DBLoc.GetXaxis().SetTitleSize(0.05)
    DBLoc.GetXaxis().SetLabelSize(0.05)
    DBLoc.GetYaxis().SetTitleSize(0.05)
    DBLoc.GetYaxis().SetLabelSize(0.05)
    DBLoc.GetZaxis().SetTitleSize(0.05)
    DBLoc.GetZaxis().SetLabelSize(0.05)
    DBLoc.GetZaxis().SetTitleOffset(0.85)
    DBLoc.GetYaxis().SetTitleOffset(1.4)
    DBLoc.Draw("colz")
    DBLoc.GetYaxis().SetRangeUser(0,300)
    DBLoc.GetXaxis().SetRangeUser(300,600)
#    leg = rt.TLegend(0.35,0.7,0.7,0.8)
#    leg.SetBorderSize(0)
#    leg.AddEntry(DBLoc, "Before Dark Bremsstrahlung","l")
#    leg.AddEntry(MuFinal, "After Dark Bremsstrahlung", "l")
#    leg.Draw()

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
    tag2 = rt.TLatex(0.13,0.94,"CMS")
    tag2.SetNDC()
    tag2.SetTextFont(62)
    tag3 = rt.TLatex(0.27,0.94,"Preliminary")
    tag3.SetNDC()
    tag3.SetTextFont(52)
    tag2.SetTextSize(0.065)
    tag3.SetTextSize(0.05)
#    tag1.Draw()
    tag2.Draw()
    tag3.Draw()

    c.Print(options.outDir+"DBlocation.png")
    c.Print(options.outDir+"DBlocation.pdf")
    c.Print(options.outDir+"DBlocation.C")


    rootFile = rt.TFile.Open('DBlocation.root','update')
    tdirectory = rootFile.GetDirectory('Data')
    if tdirectory==None:
        rootFile.mkdir('Data')
        tdirectory = rootFile.GetDirectory('Data')
    
    tdirectory.cd()
    c.Write()
