import os,ROOT
import mylib
import canvas_margin
import tdrstyle
import CMS_lumi, tdrstyle
import math
from array import array


## SampleGroup ##
class SampleGroup:
  def __init__(self, Name, Type, Samples, Year, Color=0, Style=1, TLatexAlias="", LatexAlias=""):
    self.Name = Name
    self.Type = Type
    self.Samples = Samples
    self.Year = Year
    self.Color = Color
    self.Style = Style
    self.TLatexAlias = TLatexAlias
    self.LatexAlias = LatexAlias

  def Print(self):
    print 'Sample group name = '+self.Name
    print '  Type = '+self.Type
    print '  Samples = ',
    print self.Samples
    print '  Year = '+str(self.Year)
    print '  Color = '+str(self.Color)
    print '  TLatexAlias = '+str(self.TLatexAlias)
    print '  LatexAlias = '+str(self.LatexAlias)

## DBremSignalInfo ##
class DBremSignalInfo:
  def __init__(self, mDBrem):
    self.mDBrem = mDBrem
    self.xsec = 1.
    self.xsecScale = 1.
    self.Color = ROOT.kBlack
    self.Style = 3
    self.useOfficial = False
    self.mDBrem_FileString = str(mDBrem).replace('.','p')

    print "DEFINING DBremSignalInfo"
    print "self.mDBrem: ", self.mDBrem
    self.TLatexAlias = "m_{A'} = %2.1lf GeV"%(self.mDBrem)

  def Print(self):
    print '(%d, %d, %f, %f)'%(self.mDBrem, self.xsec, self.kfactor)
  def GetTLatexAlias(self):
    if self.xsecScale!=1.:
      self.TLatexAlias = "%2.0f"%(self.xsecScale)
      print "self.TLatexAlias: ", self.TLatexAlias
      #self.TLatexAlias.replace("e-0","#times 10^{-")
#      self.TLatexAlias = self.TLatexAlias.replace("e-0","#times 10^{-")
#      print "self.TLatexAlias: ", self.TLatexAlias
#      self.TLatexAlias += "}"
      self.TLatexAlias += " #times m_{A'} = %2.1lf GeV"%(self.mDBrem)

    return self.TLatexAlias

## Variable ##
class Variable:
  def __init__(self, Name, TLatexAlias, Unit):
    self.Name = Name
    self.TLatexAlias = TLatexAlias
    self.Unit = Unit
  def Print(self):
    print '(%s, %s, %s)' % (self.Name, self.TLatexAlias, self.Unit)

## Region ##
class Region:
  def __init__(self, Name, UnblindData=True, Logy=-1, TLatexAlias=""):
    self.Name = Name
    self.UnblindData = UnblindData
    self.Logy = Logy
    self.TLatexAlias = TLatexAlias

    self.DrawData = True
    self.DrawRatio = True

  def Print(self):
    print '(%s, UnblindData=%s, Logy=%f, %s)'%(self.Name, self.UnblindData, self.Logy, self.TLatexAlias)

## Systematic ##
class Systematic:
  def __init__(self, Name, Direction, Year):
    self.Name = Name
    self.Direction = Direction
    self.Year = Year ## if <0, it's correalted
  def FullName(self):
    if self.Year>0:
      return 'Run'+str(self.Year)+'_'+self.Name
    else:
      return self.Name
  def Print(self):
    print '(%s, %s, %d)'%(self.Name, self.Direction, self.Year)

## Plotter ##
class Plotter:

  def __init__(self):

    self.DoDebug = False

    self.DataYear = 2016
    self.DataDirectory = "2016"

    self.SampleGroups = []
    self.RegionsToDraw = []
    self.VariablesToDraw = []
    self.SignalsToDraw = []
    self.Signals = ["2p0","1p0","0p4","0p2"]

    self.Systematics = []
    self.InputDirectory = ""
    self.Filename_prefix = ""
    self.Filename_suffix = ""
    self.Filename_skim = ""
    self.OutputDirectory = ""

    self.ScaleMC = False

    self.ExtraLines = ""

    self.ErrorFromShape = False
    self.AddErrorLinear = False

    self.NoErrorBand = False

  def PrintBorder(self):
    print '--------------------------------------------------------------------------'
  def PrintSamples(self):
    self.PrintBorder()
    print '[Plotter.PrintSamples()] Printing samples'
    for s in self.SampleGroups:
      s.Print()
    self.PrintBorder()

  def PrintRegions(self):
    self.PrintBorder()
    print '[Plotter.PrintRegions()] Printing regions to be drawn'
    for s in self.RegionsToDraw:
      s.Print()
    self.PrintBorder()

  def PrintVariables(self):
    self.PrintBorder()
    print '[Plotter.PrintVariables()] Printing variables to be drawn'
    for s in self.VariablesToDraw:
      s.Print()
    self.PrintBorder()

  def PrintSystematics(self):
    self.PrintBorder()
    print '[Plotter.PrintSystematics()] Printing systematics'
    for s in self.Systematics:
      s.Print()
    self.PrintBorder()

  ## Benning related

  def SetBinningFilepath(self, RebinFilepath, XaxisFilepath, YaxisFilepath):
    self.RebinFilepath = RebinFilepath
    self.XaxisFilepath = XaxisFilepath
    self.YaxisFilepath = YaxisFilepath
  def ReadBinningInfo(self, Region):
    ## Rebin
    Rebins = dict()
    for line in open(self.RebinFilepath).readlines():
      words = line.split()
      if Region!=words[0]:
        continue
      Rebins[words[1]] = int(words[2])
    ## xaxis
    XaxisRanges = dict()
    for line in open(self.XaxisFilepath).readlines():
      words = line.split()
      if Region!=words[0]:
        continue
      XaxisRanges[words[1]] = [float(words[2]), float(words[3])]

    return Rebins, XaxisRanges
  def Rebin(self, hist, region, var, nRebin):
    if var=='WRCand_Mass':
      return mylib.RebinWRMass(hist, region, self.DataYear)
    elif var=='ToBeCorrected_Jet_Pt':
      return mylib.RebinJetPt(hist, region, self.DataYear)
#    elif var=='hcalEnergy':
#      print "Rebinning hcalEnergy"
#      return mylib.RebinHcalEnergy(hist)
    else:
      if nRebin>0:
        hist.Rebin(nRebin)
        return hist
      else:
        return hist
  def ZeroDataCheckCut(self,var,xlow,xhigh):
    ## TODO
    return False

  def Draw(self):

   ROOT.gErrorIgnoreLevel = ROOT.kFatal

   tdrstyle.setTDRStyle()
   ROOT.TH1.AddDirectory(False)

   sampleOrder = []
   Signal = '0p2'

   c1 = ROOT.TCanvas('c1', '', 1300, 800)

   c1_up_left = ROOT.TPad("c1_up_left", "", 0, 0.25, 0.75, 1)
   c1_down_left = ROOT.TPad("c1_down_left", "", 0, 0, 0.75, 0.25)
   c1_up_right = ROOT.TPad("c1_up_right", "", 0.75, 0.25, 1, 1)
   c1_down_right = ROOT.TPad("c1_down_right", "", 0.75, 0, 1, 0.25)
   c1, c1_up_left, c1_down_left, c1_up_right, c1_down_right = canvas_margin.canvas_margin_4(c1, c1_up_left, c1_down_left, c1_up_right, c1_down_right)

   c1.Draw()
   c1_up_left.Draw()
   c1_down_left.Draw()
   c1_up_right.Draw()
   c1_down_right.Draw()

   for Region in self.RegionsToDraw:
    #for Signal in self.Signals:

    print '## Drawing '+Region.Name
    if('partialDisappearance' in Region.Name): 
#        sampleOrder = [2,1,0,3]
      sampleOrder = [0,1,2,3,4]
      self.SampleGroups = [self.SampleGroups[i] for i in sampleOrder]
    elif('completeDisappearance' in Region.Name): 
      sampleOrder = [3,2,1,0,4]
#        sampleOrder = [0,1,2,3]
      self.SampleGroups = [self.SampleGroups[i] for i in sampleOrder]

      ## Read binning data
      #Rebins, XaxisRanges = self.ReadBinningInfo(Region.Name)

      ## Input/Output directotry
    Indir = self.InputDirectory
    Outdir = self.OutputDirectory+'/'+Region.Name+'/'
    if self.ScaleMC:
      Outdir = self.OutputDirectory+'/ScaleMC/'+Region.Name+'/'
    print '##   Outputs => '+Outdir
    os.system('mkdir -p '+Outdir)

    ## Data file
    f_Data = ROOT.TFile(Indir+'/'+self.Filename_prefix+self.Filename_skim+self.Filename_suffix+'DY_'+Signal+'.root')

    ## Loop over variables

    for Variable in self.VariablesToDraw:
      if ('partialDisappearance' in Region.Name) and ('nEvents' in Variable.Name): continue
        #if ('completeDisappearance' in Region.Name) and ('staDEoverE' in Variable.Name): continue
      if ('completeDisappearance' in Region.Name) and ('hcalEnergy' in Variable.Name): continue

      ## BinInfo
      nRebin = 0
      #xMin = 0
      #xMax = 1
      xMin = 80
      xMax = 100
      xMin_complete = 80
      xMax_complete = 100
      if ('partialDisappearanceLargeStaDPhi' in Region.Name):
        xMin = 0
	xMax = 20
      elif ('partialDisappearanceBDT' in Region.Name):
        xMin = 0.98
	xMax = 1
      elif ('partialDisappearance' in Region.Name):
        xMin = -1
        #nRebin = Rebins[Variable.Name]
        #xMin = XaxisRanges[Variable.Name][0]
        #xMax = XaxisRanges[Variable.Name][1]
      yMax = -999

      if self.DoDebug:
          print '[DEBUG] Trying to draw variable = '+Variable.Name
          print '[DEBUG] (xMin,xMax) = (%s,%s)'%(xMin,xMax)

        ## xtitle
      xtitle_complete = 0
      xtitle = 0
      if 'completeDisappearance' in Region.Name:
        xtitle_complete = 'Invariant Mass of Tag and Probe (GeV)'
      else:
        xtitle = Variable.TLatexAlias

        ## Save hists
        ## For legend later..
      HistsToDraw_complete = dict()
      HistsToDraw = dict()

        ## Get data first
      if self.DoDebug:
          print '[DEBUG] Trying to get data histogram..'
        #h_Data = f_Data.Get(Region.Name+'/'+Variable.Name+'_'+Region.Name)
        ## Grabbing background temporarily

      print(self.Filename_prefix+Region.Name+'/'+self.Filename_prefix+Region.Name) 
      if 'completeDisappearance' in Region.Name:
        h_Data_complete = f_Data.Get(self.Filename_prefix+Region.Name+'/'+self.Filename_prefix+Region.Name)
        if not h_Data_complete:
          print f_Data.GetName()
          print 'data_obs'+' ==> No data, skipped'
          continue

        ## Make overflow
        h_Data_complete.GetXaxis().SetRangeUser(xMin_complete,xMax_complete)
        #h_Data = mylib.MakeOverflowBin(h_Data)

        ## Rebin
	print "nBins: ", h_Data_complete.GetXaxis().GetNbins()
        h_Data_complete = self.Rebin(h_Data_complete, Region.Name, Variable.Name, nRebin)
        ## Copy data axis
        dataAxis = h_Data_complete.GetXaxis()
        nBin = dataAxis.GetNbins()
        xBins = [dataAxis.GetBinLowEdge(1)]
        if ("partialDisappearanceLargeStaDPhi" in Region.Name):
	  xBins = [0,1,2,2.5,3.,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.5,10,20]
	  nBin = len(xBins) - 1
	else:
          for ix in range(0,nBin):
	    print "dataAxis.GetBinUpEdge(ix+1): ", dataAxis.GetBinUpEdge(ix+1)
            xBins.append( dataAxis.GetBinUpEdge(ix+1) )
        xBins = array("d",xBins)

        ### Att data histogram
        h_Data_complete.SetMarkerStyle(20)
        h_Data_complete.SetMarkerSize(1.2)
        h_Data_complete.SetMarkerColor(ROOT.kBlack)
        h_Data_complete.SetLineColor(ROOT.kBlack)

        if self.DoDebug:
          print '[DEBUG] data histogram finished'
          print 'Data:'
          for z in range(0,h_Data_complete.GetXaxis().GetNbins()):
            x_l = h_Data_complete.GetXaxis().GetBinLowEdge(z+1)
            x_r = h_Data_complete.GetXaxis().GetBinUpEdge(z+1)
            y = h_Data_complete.GetBinContent(z+1)
            print '[%f,%f] : %f +- %f'%(x_l,x_r,y,h_Data_complete.GetBinError(z+1))

        ## Loop over samples
        ## For Legend, save 
        HistsForLegend = []
        AliasForLegend = [] ## Prevent double-counting
        stack_Bkgd_complete = ROOT.THStack("stack_Bkgd_complete", "")
        h_Bkgd_complete = 0
        ## Save systematic
        SystematicUps_complete = dict()
        SystematicDowns_complete = dict()
        ## If we take errors fro shapes
        h_TotalBackgroundFromShape_complete = 0

        ## Loop over systematics
        for Syst in self.Systematics:

          if self.DoDebug:
            print '[DEBUG] Trying to make histogram for Syst = ',
            Syst.Print()

          h_Bkgd_ForSyst_complete = 0
          dirName = Region.Name

          if Syst.Name!="Central":

            if Syst.Direction>0:
              dirName = "Syst_"+Syst.Name+"Up_"+Region.Name
            else:
              dirName = "Syst_"+Syst.Name+"Down_"+Region.Name

          for SampleGroup in self.SampleGroups:
            Color = SampleGroup.Color
            LegendAdded = False

            if self.DoDebug:
	      print "[DEBUG] SampleGroup.Samples: ", SampleGroup.Samples
            for Sample in SampleGroup.Samples:
              if self.DoDebug:
	        print "[DEBUG] Sample: ", Sample
                print '[DEBUG] Trying to make histogram for Sample = '+Sample

              f_Sample = ROOT.TFile(Indir+'/'+self.Filename_prefix+self.Filename_skim+self.Filename_suffix+Sample+'_'+Signal+'.root')
              h_Sample_complete = 0

              ## Uncorrelated sources has Syst.Year = 2016 or 2017 or 2018
              ## For this cases, SampleGroup.Year should be matched
	      if (Syst.Name != 'Central'):
		print "Histo name: ", Sample+'_'+Syst.Name+str(Syst.Direction)
		h_Sample_complete = f_Sample.Get(Sample+'_'+Syst.Name+Syst.Direction)
              ## 2) Lumi, MC normalizaion
              ## Use central and scale them later
              elif (Syst.Name in ["Lumi", "DYNorm", "NonPromptNorm", "OthersNorm"]):
                tmp_dirName = Region.Name
                h_Sample_complete = f_Sample.Get(tmp_dirName+'/'+Variable.Name+'_'+tmp_dirName)
              ## For all other cases
              else:
		print "Histo name: ", self.Filename_prefix+Region.Name+'/'+self.Filename_prefix+Region.Name 
                h_Sample_complete = f_Sample.Get(self.Filename_prefix+Region.Name+'/'+self.Filename_prefix+Region.Name)

              if not h_Sample_complete:
                print 'No hist : %s %s'%(Syst.Name,Sample)
                continue

              print "h_Sample_complete.GetBinContent(1): ", h_Sample_complete.GetBinContent(1)
	      print "h_Sample_complete.GetXaxis.GetNbins(): ", h_Sample_complete.GetXaxis().GetNbins()
	      print "h_Sample_complete.GetXaxis().GetBinUpEdge(15): ", h_Sample_complete.GetXaxis().GetBinUpEdge(15)
              ## Make overflow
              h_Sample_complete.GetXaxis().SetRangeUser(xMin_complete,xMax_complete)
              #h_Sample = mylib.MakeOverflowBin(h_Sample)

              h_Sample_complete = self.Rebin(h_Sample_complete, Region.Name, Variable.Name, nRebin)
              h_Sample_complete.SetLineColor(Color)
              h_Sample_complete.SetLineWidth(1)
              h_Sample_complete.SetFillColor(Color)

              ## Scale
              MCSF, MCSFerr = 1., 0.
              if self.ScaleMC:
                ## now, only for DY
                if "DYJets" in Sample:
                  MCSF, MCSFerr = mylib.GetDYNormSF(SampleGroup.Year, Region.Name)
              h_Sample_complete.Scale( MCSF )

              ## Manual systematic
              ## 1) [Lumi] Uncorrelated
              if (Syst.Name=="Lumi") and (Syst.Year==SampleGroup.Year):
                lumierr = mylib.LumiError(SampleGroup.Year)
                for ix in range(0,h_Sample_complete.GetXaxis().GetNbins()):
                  y = h_Sample_complete.GetBinContent(ix+1)
                  y_new = y + y*float(Syst.Direction)*lumierr
                  h_Sample_complete.SetBinContent(ix+1, y_new)
              ## 2) [DYNorm] Uncorrelated, only for DY
              if (Syst.Name=="DYNorm") and ("DYJets" in Sample) and (Syst.Year==SampleGroup.Year):
                for ix in range(0,h_Sample_complete.GetXaxis().GetNbins()):
                  y = h_Sample_complete.GetBinContent(ix+1) ## h_Sample is already scaled by MCSF above
                  y_new = y * ( MCSF + float(Syst.Direction)*MCSFerr ) / MCSF
                  h_Sample_complete.SetBinContent(ix+1, y_new)

              ## AddError option
              AddErrorOption = ''
              if self.AddErrorLinear:
                AddErrorOption = 'L'

              ## If central, add to h_Bkgd
              if self.DoDebug:
	        print "[DEBUG] Syst.Name: ", Syst.Name
	        print "[DEBUG] Sample: ", Sample
	      print "[DEBUG] Syst.Name: ", Syst.Name
	      print "[DEBUG] Sample: ", Sample
		
              if Syst.Name=="Central" and Sample=='total_background':
                if not h_TotalBackgroundFromShape_complete:
                  h_TotalBackgroundFromShape_complete = h_Sample_complete.Clone()
                else:
                  h_TotalBackgroundFromShape_complete = mylib.AddHistograms( h_TotalBackgroundFromShape_complete, h_Sample_complete, AddErrorOption )
                HistsToDraw_complete[Sample] = h_Sample_complete.Clone()

              elif Syst.Name=="Central":

                print "ADDING Sample:", Sample
                stack_Bkgd_complete.Add( h_Sample_complete )
                if not h_Bkgd_complete:
                  print "Cloning sample"
                  h_Bkgd_complete = h_Sample_complete.Clone()
                else:
                  print "Adding sample"
                  print "h_Sample_complete.GetBinContent(4): ", h_Sample_complete.GetBinContent(4)
                  h_Bkgd_complete = mylib.AddHistograms( h_Bkgd_complete, h_Sample_complete, AddErrorOption)

                print "h_Bkgd.GetBinContent(4): ", h_Bkgd_complete.GetBinContent(4)
                HistsToDraw_complete[Sample] = h_Sample_complete.Clone()
                if (not LegendAdded) and (SampleGroup.TLatexAlias not in AliasForLegend):
                  HistsForLegend.append( [HistsToDraw_complete[Sample],SampleGroup.TLatexAlias] )
                  AliasForLegend.append(SampleGroup.TLatexAlias)
                  LegendAdded = True
              ## else (i.e., systematic), add to h_Bkgd_ForSyst
              else:

                if not h_Bkgd_ForSyst_complete:
                  h_Bkgd_ForSyst_complete = h_Sample_complete.Clone()
                else:
                  h_Bkgd_ForSyst_complete = mylib.AddHistograms(h_Bkgd_ForSyst_complete, h_Sample_complete, AddErrorOption)

              ## Close file
              f_Sample.Close()

            ##==>End Sample loop
          ##==>End SampleGroup loop

          if Syst.Name!="Central":
            if Syst.Direction=="Up":
              SystematicUps_complete[Syst.FullName()] = h_Bkgd_ForSyst_complete
            else:
              SystematicDowns_complete[Syst.FullName()] = h_Bkgd_ForSyst_complete

        ##==>End Syst loop
        #print SystematicUps
        #print SystematicDowns

        ## Syst Up/Down . Max/Min
        h_Bkgd_TotErr_Max_complete = h_Bkgd_complete.Clone()
        h_Bkgd_TotErr_Min_complete = h_Bkgd_complete.Clone()
        for SystKey in SystematicUps_complete:
          h_Up = SystematicUps_complete[SystKey]
          h_Down = SystematicDowns_complete[SystKey]
	  print "SystKey: ", SystKey
          for ix in range(0,h_Bkgd_complete.GetXaxis().GetNbins()):
            x_l = h_Bkgd_complete.GetXaxis().GetBinLowEdge(ix+1)
            x_r = h_Bkgd.GetXaxis().GetBinUpEdge(ix+1)
            y_Central = h_Bkgd_complete.GetBinContent(ix+1)
            y_Up = h_Up.GetBinContent(ix+1)
            y_Down = h_Down.GetBinContent(ix+1)
            ## -.
            y_Max = max( max(y_Central,y_Up), y_Down)
            y_Min = min( min(y_Central,y_Up), y_Down)
            #print '[%d,%d] : %f, (Max,Min) = (%f,%f)'%(x_l,x_r,y_Central,y_Up,y_Down)

            ## Update Max
            err_Max_Current = h_Bkgd_TotErr_Max_complete.GetBinError(ix+1)
            err_Max_ToAdd = y_Max-y_Central
            err_Max_Updated = math.sqrt( err_Max_Current*err_Max_Current + err_Max_ToAdd*err_Max_ToAdd )
            h_Bkgd_TotErr_Max_complete.SetBinError(ix+1, err_Max_Updated)
            ## Update Min
            err_Min_Current = h_Bkgd_TotErr_Min_complete.GetBinError(ix+1)
            err_Min_ToAdd = y_Central-y_Min
            err_Min_Updated = math.sqrt( err_Min_Current*err_Min_Current + err_Min_ToAdd*err_Min_ToAdd )
            h_Bkgd_TotErr_Min_complete.SetBinError(ix+1, err_Min_Updated)
        ##==>End Systematic loop

        ## Blind mode
        if not Region.UnblindData:
          print "h_Bkgd_complete.GetBinContent(1): ", h_Bkgd_complete.GetBinContent(1)
          print "h_Bkgd_complete.GetBinError(1): ", h_Bkgd_complete.GetBinError(1)
          h_Data_complete = h_Bkgd_complete.Clone(h_Data_complete.GetName())
          for ibin in range(1,h_Data_complete.GetNbinsX()+1):
              print "h_Data_complete.GetBinError(ibin): ", h_Data_complete.GetBinError(ibin)
              h_Data_complete.SetBinError(ibin,h_Bkgd_TotErr_Max_complete.GetBinError(ibin))
              print "h_Data_complete.GetBinError(ibin): ", h_Data_complete.GetBinError(ibin)
          print "h_Data.GetBinContent(1): ", h_Data_complete.GetBinContent(1)
          h_Data_complete.SetMarkerStyle(20)
          h_Data_complete.SetMarkerSize(1.2)
          h_Data_complete.SetMarkerColor(ROOT.kBlack)
          h_Data_complete.SetLineColor(ROOT.kBlack)
          h_Data_complete.SetLineWidth(2)

        ## hist => gr
	print "h_Bkgd_TotErr_Min_complete.GetBinError(1): ", h_Bkgd_TotErr_Min_complete.GetBinError(1)
	print "h_Bkgd_TotErr_Max_complete.GetBinError(1): ", h_Bkgd_TotErr_Max_complete.GetBinError(1)
        print "h_Data_complete.GetBinContent(1): ", h_Data_complete.GetBinContent(1)
        gr_Bkgd_TotErr_complete = mylib.GetAsymmError(h_Bkgd_TotErr_Max_complete,h_Bkgd_TotErr_Min_complete)
        gr_Data_complete = ROOT.TGraphAsymmErrors(h_Data_complete)

        ## Error from shape
        if self.ErrorFromShape:
          if self.DoDebug:
            print 'Total background :'
            for z in range(0,h_TotalBackgroundFromShape_complete.GetXaxis().GetNbins()):
              x_l = h_TotalBackgroundFromShape_complete.GetXaxis().GetBinLowEdge(z+1)
              x_r = h_TotalBackgroundFromShape_complete.GetXaxis().GetBinUpEdge(z+1)
              y = h_TotalBackgroundFromShape.GetBinContent(z+1)
              print '[%f,%f] : %f +- %f'%(x_l,x_r,y,h_TotalBackgroundFromShape_complete.GetBinError(z+1))
	  #print "h_TotalBackgroundFromShape.GetBinError(1): ", h_TotalBackgroundFromShape.GetBinError(1)
          gr_Bkgd_TotErr_complete = ROOT.TGraphAsymmErrors(h_TotalBackgroundFromShape_complete)

        err_up_tmp = []
        err_down_tmp = []
        alpha = 1. - 0.6827
        for i in range(0, gr_Data_complete.GetN()):
          N = gr_Data_complete.GetY()[i]

          L = 0.                                          if (N==0.) else (ROOT.Math.gamma_quantile(alpha/2.,N,1.))
          U = ( ROOT.Math.gamma_quantile_c(alpha,N+1,1) ) if (N==0.) else (ROOT.Math.gamma_quantile_c(alpha/2.,N+1.,1.))

          #print '%d - %f + %f'%(N, N-L, U-N)

          if N!=0:
            gr_Data_complete.SetPointEYlow(i, N-L )
            gr_Data_complete.SetPointEYhigh(i, U-N )
            err_down_tmp.append(N-L)
            err_up_tmp.append(U-N)
          else:
            zerodata_err_low = 0.1
            zerodata_err_high = 1.8

            xlow = gr_Data_complete.GetX()[i]-gr_Data_complete.GetEXlow()[i]
            xhigh = gr_Data_complete.GetX()[i]+gr_Data_complete.GetEXhigh()[i]
            if self.ZeroDataCheckCut(Variable.Name,xlow,xhigh):
              zerodata_err_low = 0.
              zerodata_err_high = 0.

            gr_Data_complete.SetPointEYlow(i, zerodata_err_low)
            gr_Data_complete.SetPointEYhigh(i, zerodata_err_high)
            err_down_tmp.append(zerodata_err_low)
            err_up_tmp.append(zerodata_err_high)

        ## Legend
        lg_complete = 0
        ## No signal
        if len(self.SignalsToDraw)==0:
          lg_complete = ROOT.TLegend(0.55, 0.45, 0.92, 0.90)
        ## With Signal
        else:
          if Region.DrawRatio:
            lg_complete = ROOT.TLegend(0.55, 0.56, 0.92, 0.90)
          else:
            lg_complete = ROOT.TLegend(0.55, 0.56, 0.92, 0.90)
        lg_complete.SetBorderSize(0)
        lg_complete.SetFillStyle(0)

        if not self.NoErrorBand:
          lg_complete.AddEntry(gr_Bkgd_TotErr_complete, "Stat. #oplus syst. uncert.", "f")
        ## dummy graph for legend..
        ## this is because h_Data does not have horizontal error bars,
        ## and gr_data does not have points
        gr_Data_dummy_complete = ROOT.TGraphAsymmErrors(gr_Data_complete)
        gr_Data_dummy_complete.SetMarkerStyle(20)
        gr_Data_dummy_complete.SetMarkerSize(1.2)
        dataLegendGOption="ep"
        if Variable.Name=="WRCand_Mass":
          dataLegendGOption="lpe"
        if Region.DrawData:
          if Region.UnblindData:
            lg_complete.AddEntry(gr_Data_dummy_complete, "Data", dataLegendGOption)
          else:
            lg_complete.AddEntry(gr_Data_dummy_complete, "Total background", dataLegendGOption)
        for i_lg in range(0,len(HistsForLegend)):
          h_lg = HistsForLegend[ len(HistsForLegend)-1-i_lg ][0]
          tlatexaliax = HistsForLegend[ len(HistsForLegend)-1-i_lg ][1]
          lg_complete.AddEntry(h_lg,tlatexaliax,"f")

        ## Prepare canvas

        if self.DoDebug:
          print '[DEBUG] Preparing canvas..'

        c1_up_left.cd()
        if Region.Logy>0:
          c1_up_left.SetLogy(True)

        c1.cd()

        #### axis histograms

        h_dummy_up_complete = ROOT.TH1D('h_dummy_up_complete', '', nBin, xBins)
        h_dummy_up_complete.GetXaxis().SetRangeUser(xMin_complete, xMax_complete)
        if nRebin>0:
          binsize = h_dummy_up_complete.GetXaxis().GetBinUpEdge(1)-h_dummy_up_complete.GetXaxis().GetBinLowEdge(1)
          str_binsize = '%d'%(binsize)
          if binsize!=int(binsize):
            str_binsize = '%1.2f'%(binsize)
          h_dummy_up_complete.GetYaxis().SetTitle('Events / '+str_binsize+' '+Variable.Unit)
        else:
          h_dummy_up_complete.GetYaxis().SetTitle('Events / bin')
        if Variable.Name=='WRCand_Mass':
          h_dummy_up_complete.GetYaxis().SetTitle('Events / bin')

        h_dummy_down_complete = ROOT.TH1D('h_dummy_down_complete', '', nBin, xBins)
        h_dummy_down_complete.GetYaxis().SetRangeUser(0.,2.0)

        #if ('completeDisappearance' in Region.Name):
        #  h_dummy_down_complete.GetYaxis().SetRangeUser(0.40,1.90)

        h_dummy_down_complete.SetNdivisions(504,"Y")
        h_dummy_down_complete.SetTickLength(0)
        #h_dummy_down_complete.SetNdivisions(5,"X")
        #if ('completeDisappearance' in Region.Name):
          #h_dummy_down.GetXaxis().SetNdivisions(1)
          #h_dummy_down.GetXaxis().SetLabelOffset(100)
          #h_dummy_down.GetXaxis().CenterTitle(True)
        h_dummy_down_complete.GetXaxis().SetRangeUser(xMin_complete, xMax_complete)
        h_dummy_down_complete.GetXaxis().SetTitle(xtitle_complete)
        h_dummy_down_complete.GetYaxis().SetTitle("#frac{Data}{Pred.}")
        h_dummy_down_complete.SetFillColor(0)
        h_dummy_down_complete.SetMarkerSize(0)
        h_dummy_down_complete.SetMarkerStyle(0)
        h_dummy_down_complete.SetLineColor(ROOT.kWhite)

        if Region.DrawRatio:
          h_dummy_up_complete, h_dummy_down_complete = canvas_margin.hist_axis(h_dummy_up_complete, h_dummy_down_complete)
        else:
          h_dummy_up_complete.SetTitle("")
          h_dummy_up_complete.GetYaxis().SetLabelSize(0.04)
          h_dummy_up_complete.GetYaxis().SetTitleSize(0.054)
          #h_dummy_up_complete.GetYaxis().SetTitleOffset(1.30)
          h_dummy_up_complete.GetXaxis().SetLabelSize(0.035)
          h_dummy_up_complete.GetXaxis().SetTitleSize(0.055)
         #h_dummy_up_complete.GetXaxis().SetTitle(xtitle_complete)
          h_dummy_up_complete.GetXaxis().SetTitle("Comp. Dis. Region")

        h_dummy_up_complete.SetTickLength(0)
        h_dummy_up_complete.GetXaxis().CenterTitle(True)
        #h_dummy_up_complete.SetNdivisions(0,"X")
        h_dummy_up_complete.GetXaxis().SetTitleOffset(0.5)
        h_dummy_up_complete.GetXaxis().SetTitle("Complete Dis. Region")
        h_dummy_down_complete.GetXaxis().SetLabelOffset(0.03)
        #h_dummy_down_complete.GetXaxis().SetTitle(xtitle_complete)

        ## Get Y maximum
        yMax = max( yMax, mylib.GetMaximum(gr_Data_complete) )
        yMax = max( yMax, mylib.GetMaximum(gr_Bkgd_TotErr_complete) )
        ## Yaxis range
        yMin = 0.
        yMaxScale = 1.2
        if Region.Logy>0:
          yMaxScale = 10
          yMin = 0.2
        #print "!!!! yMin: ", yMin
        #h_dummy_up.GetYaxis().SetRangeUser( yMin, yMaxScale*yMax )
        h_dummy_up_complete.GetYaxis().SetRangeUser( yMin, 1000 )

        ## Exception control
	print "Variable.Name: ", Variable.Name

        if self.DoDebug:
          print '[DEBUG] Canvas is ready'

        ## Draw up
        c1_up_left.cd()

        print "[DEBUG] h_dummy_up_complete.Integral(): ", h_dummy_up_complete.Integral()
        h_dummy_up_complete.Draw("hist")
        stack_Bkgd_complete.Draw("histsame")

        gr_Bkgd_TotErr_complete.SetMarkerColor(0)
        gr_Bkgd_TotErr_complete.SetMarkerSize(0)
        gr_Bkgd_TotErr_complete.SetFillStyle(3013)
        gr_Bkgd_TotErr_complete.SetFillColor(ROOT.kBlack)
        gr_Bkgd_TotErr_complete.SetLineColor(0)
	print "self.NoErrorBand: ", self.NoErrorBand
        if not self.NoErrorBand:
          gr_Bkgd_TotErr_complete.Draw("sameE2")

        gr_Data_complete.SetLineWidth(2)
        gr_Data_complete.SetMarkerSize(0.)
        gr_Data_complete.SetMarkerColor(ROOT.kBlack)
        gr_Data_complete.SetLineColor(ROOT.kBlack)
        if Region.DrawData:
	  print "!!!!!!!!!!!!!!!!! h_Data.Integral(): ", h_Data_complete.Integral()
          print "h_Data_complete.GetBinError(1): ", h_Data_complete.GetBinError(1)
          h_Data_complete.Draw("p EX0 same")
          #gr_Data_complete.Draw("p0same")


        ## Signal
        h_Sigs = []
        print "self.SignalsToDraw: ", self.SignalsToDraw
        for Sig in self.SignalsToDraw:
	  print "Sig.mDBrem_FileString: ", Sig.mDBrem_FileString
	  #if Signal not in Sig.mDBrem_FileString: continue
	  print "Sig.mDBrem_FileString: ", Sig.mDBrem_FileString

          #fpullpath_Sig = Indir+'/'+self.Filename_prefix+self.Filename_skim+self.Filename_suffix+'DBrem_'+Sig.mDBrem_FileString+'_'+Signal+'.root'
          fpullpath_Sig = Indir+'/'+self.Filename_prefix+self.Filename_skim+self.Filename_suffix+'DBrem_'+Sig.mDBrem_FileString+'_'+Sig.mDBrem_FileString+'.root'
          f_Sig = ROOT.TFile(fpullpath_Sig)
          print fpullpath_Sig
          print self.Filename_prefix+Region.Name+'/DBrem_'+Region.Name
          h_Sig = f_Sig.Get(self.Filename_prefix+Region.Name+'/DBrem_'+Region.Name)
          if not h_Sig:
            print fpullpath_Sig
            print self.Filename_prefix+Region.Name+'/DBrem_'+Region.Name
            continue

          print "h_Sig.GetBinContent(1): ", h_Sig.GetBinContent(1)
          ## Make overflow
          h_Sig.GetXaxis().SetRangeUser(xMin_complete,xMax_complete)
          #h_Sig = mylib.MakeOverflowBin(h_Sig)

          ## Rebin
          h_Sig = self.Rebin(h_Sig, Region.Name, Variable.Name, nRebin)
          print "2nd check h_Sig.GetBinContent(1): ", h_Sig.GetBinContent(1)
          print "h_Sig.Integral(): ", h_Sig.Integral()

          ## Scale
          print "Scaling h_Sig"
          h_Sig.Scale( Sig.xsec * Sig.xsecScale )
          print "3rd check h_Sig.GetBinContent(1): ", h_Sig.GetBinContent(1)

          ## Att
          h_Sig.SetLineWidth(3)
          h_Sig.SetLineColor(Sig.Color)
          h_Sig.SetLineStyle(Sig.Style)

          ## Save
          h_Sigs.append(h_Sig)

          ## legend
          lg_complete.AddEntry(h_Sig, Sig.GetTLatexAlias(), 'l')
          lg_complete.SetTextFont(42)

          ## Draw
	  print "!!!!!!!!!!!!!!!!! h_Sig.Integral(): ", h_Sig.Integral()
          h_Sig.Draw("histsame")

        h_dummy_up_complete.Draw("axissame")

        ## Legend
        lg_complete.Draw()

        ## Draw down
        c1_down_left.cd()
        #h_dummy_down_complete.GetYaxis().SetTitle("#frac{Data}{Pred.}")
        print "h_dummy_down_complete.GetXaxis().GetTitle(): ", h_dummy_down_complete.GetXaxis().GetTitle()
        h_dummy_down_complete.Draw("hist axis")
        #h_dummy_down_complete.Draw("axissame")

        ## values must be set later
        h_Data_Ratio_complete = h_Data_complete.Clone('h_Data_Ratio_complete')
        ## BinContent set by divide here, but errors must be set later
        tmp_h_Data_Ratio = h_Data_complete.Clone('tmp_h_Data_Ratio')
        tmp_h_Data_Ratio.Divide(h_Bkgd_complete)
        gr_Data_Ratio_complete = ROOT.TGraphAsymmErrors(tmp_h_Data_Ratio)
        gr_Data_Ratio_complete.SetName('gr_Data_Ratio_complete')
        gr_Data_Ratio_complete.SetLineWidth(2)
        gr_Data_Ratio_complete.SetMarkerSize(0.)
        gr_Data_Ratio_complete.SetLineColor(ROOT.kBlack)
        ## values must be set later, but BinContent will be simply 1
        gr_Bkgd_Ratio_complete = gr_Bkgd_TotErr_complete.Clone('gr_Bkgd_Ratio_complete')

        for i in range(1,h_Data_Ratio_complete.GetXaxis().GetNbins()+1):

          ## FIXME for zero? how?
          if h_Bkgd_complete.GetBinContent(i)!=0:

            ## ratio point
            ## BinContent = Data/Bkgd
            ## BinError = DataError/Bkgd
            h_Data_Ratio_complete.SetBinContent( i, h_Data_Ratio_complete.GetBinContent(i) / h_Bkgd_complete.GetBinContent(i) )
            h_Data_Ratio_complete.SetBinError ( i, h_Data_Ratio_complete.GetBinError(i) / h_Bkgd_complete.GetBinContent(i) )

            if err_down_tmp[i-1]!=0.:
              gr_Data_Ratio_complete.SetPointEYlow(i-1, err_down_tmp[i-1] / h_Bkgd_complete.GetBinContent(i) )
              gr_Data_Ratio_complete.SetPointEYhigh(i-1, err_up_tmp[i-1] / h_Bkgd_complete.GetBinContent(i))
              if Variable.Name!="WRCand_Mass":
                gr_Data_Ratio_complete.SetPointEXlow(i-1, 0)
                gr_Data_Ratio_complete.SetPointEXhigh(i-1, 0)
            else:
              gr_Data_Ratio_complete.SetPointEYlow(i-1, 0)
              gr_Data_Ratio_complete.SetPointEYhigh(i-1, 1.8 / h_Bkgd_complete.GetBinContent(i))
              if Variable.Name!="WRCand_Mass":
                gr_Data_Ratio_complete.SetPointEXlow(i-1, 0)
                gr_Data_Ratio_complete.SetPointEXhigh(i-1, 0)

            ## ratio allerr
            ## BinContent = 1
            ## BinError = Bkgd(Stat+Syst)Error/Bkgd

            gr_Bkgd_Ratio_complete.SetPoint(i-1,h_Bkgd_complete.GetXaxis().GetBinCenter(i), 1.)
            gr_Bkgd_Ratio_complete.SetPointEYhigh( i-1, gr_Bkgd_Ratio_complete.GetErrorYhigh(i-1) / h_Bkgd_complete.GetBinContent(i) )
            gr_Bkgd_Ratio_complete.SetPointEYlow( i-1,  gr_Bkgd_Ratio_complete.GetErrorYlow(i-1) / h_Bkgd_complete.GetBinContent(i) )

          elif h_Bkgd_complete.GetBinContent(i)==0. and h_Data_Ratio_complete.GetBinContent(i)==0.:

            h_Data_Ratio_complete.SetBinContent( i, 0 )
            h_Data_Ratio_complete.SetBinError ( i, 0 )
            gr_Data_Ratio_complete.SetPoint(i-1, 0, 0)
            gr_Data_Ratio_complete.SetPointEYlow(i-1, 0)
            gr_Data_Ratio_complete.SetPointEYhigh(i-1, 0)

            gr_Bkgd_Ratio_complete.SetPoint(i-1,h_Bkgd_complete.GetXaxis().GetBinCenter(i), 1.)
            gr_Bkgd_Ratio_complete.SetPointEYhigh( i-1, 0. )
            gr_Bkgd_Ratio_complete.SetPointEYlow( i-1, 0. )

            if Variable.Name!="WRCand_Mass":
              gr_Data_Ratio_complete.SetPointEXlow(i-1, 0)
              gr_Data_Ratio_complete.SetPointEXhigh(i-1, 0)

          ## If bkgd <= 0
          else:
            this_max_ratio = 20.0
            this_data = h_Data_Ratio_complete.GetBinContent(i)
            this_data_err = h_Data_Ratio_complete.GetBinError(i)

            h_Data_Ratio_complete.SetBinContent( i, this_max_ratio )
            h_Data_Ratio_complete.SetBinError ( i, this_data_err*this_max_ratio/this_data )

            tmp_x = ROOT.Double(0.)
            tmp_y = ROOT.Double(0.)
            gr_Data_Ratio_complete.GetPoint(i-1, tmp_x, tmp_y)
            gr_Data_Ratio_complete.SetPoint(i-1, tmp_x, this_max_ratio)
            gr_Data_Ratio_complete.SetPointEYlow(i-1, err_down_tmp[i-1]*this_max_ratio/this_data)
            gr_Data_Ratio_complete.SetPointEYhigh(i-1, err_up_tmp[i-1]*this_max_ratio/this_data)

            gr_Bkgd_Ratio_complete.SetPoint(i-1,h_Bkgd_complete.GetXaxis().GetBinCenter(i), 1.)
            gr_Bkgd_Ratio_complete.SetPointEYhigh( i-1, 0. )
            gr_Bkgd_Ratio_complete.SetPointEYlow( i-1, 0. )

            if Variable.Name!="WRCand_Mass":
              gr_Data_Ratio_complete.SetPointEXlow(i-1, 0)
              gr_Data_Ratio_complete.SetPointEXhigh(i-1, 0)
        ##==>End bin loop

        gr_Bkgd_Ratio_complete.SetMarkerColor(0)
        gr_Bkgd_Ratio_complete.SetMarkerSize(0)
        gr_Bkgd_Ratio_complete.SetFillStyle(3013)
        gr_Bkgd_Ratio_complete.SetFillColor(ROOT.kBlack)
        gr_Bkgd_Ratio_complete.SetLineColor(0)
        gr_Bkgd_Ratio_complete.Draw("sameE2")

        h_Data_Ratio_complete.Draw("p EX0 same")
        #gr_Data_Ratio_complete.Draw("p0same")

        ## y=1 graph
        g1_x = [-9000, 9000]
        g1_y = [1, 1]
        g1_complete = ROOT.TGraph(2, array("d", g1_x ), array("d", g1_y ))
        g1_complete.Draw("same")
        for ibin in range(1,h_Data_complete.GetNbinsX()+1):
            print "FINAL h_Data_complete.GetBinError(ibin): ", h_Data_complete.GetBinError(ibin)

        c1.cd()
        c1.SaveAs(Outdir+'Test_complete.pdf')

        ## TLatex
        c1_up_left.cd()
        channelname = ROOT.TLatex()
        channelname.SetNDC()
        channelname.SetTextSize(0.03)
        #channelname.DrawLatex(0.25, 0.88, "#font[42]{"+Region.TLatexAlias+"}")

        ## Extra lines
        #exec(self.ExtraLines)
        print "Draw self.ExtraLines"

      #################################################################################################################################
      #### NEED TO COPY ALL OF ABOVE
      elif 'partialDisappearance' in Region.Name:
        h_Data = f_Data.Get(self.Filename_prefix+Region.Name+'/'+self.Filename_prefix+Region.Name)
        if not h_Data:
          print f_Data.GetName()
          print 'data_obs'+' ==> No data, skipped'
          continue

        ## Make overflow
        h_Data.GetXaxis().SetRangeUser(xMin,xMax)
        #h_Data = mylib.MakeOverflowBin(h_Data)

        ## Rebin
	print "nBins: ", h_Data.GetXaxis().GetNbins()
        h_Data = self.Rebin(h_Data, Region.Name, Variable.Name, nRebin)
        ## Copy data axis
        dataAxis = h_Data.GetXaxis()
        nBin = dataAxis.GetNbins()
        xBins = [dataAxis.GetBinLowEdge(1)]
        if ("partialDisappearanceLargeStaDPhi" in Region.Name):
	  xBins = [0,1,2,2.5,3.,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.5,10,20]
	  nBin = len(xBins) - 1
	else:
          for ix in range(0,nBin):
	    print "dataAxis.GetBinUpEdge(ix+1): ", dataAxis.GetBinUpEdge(ix+1)
            xBins.append( dataAxis.GetBinUpEdge(ix+1) )
        xBins = array("d",xBins)

        ### Att data histogram
        h_Data.SetMarkerStyle(20)
        h_Data.SetMarkerSize(1.2)
        h_Data.SetMarkerColor(ROOT.kBlack)
        h_Data.SetLineColor(ROOT.kBlack)

        if self.DoDebug:
          print '[DEBUG] data histogram finished'
          print 'Data:'
          for z in range(0,h_Data.GetXaxis().GetNbins()):
            x_l = h_Data.GetXaxis().GetBinLowEdge(z+1)
            x_r = h_Data.GetXaxis().GetBinUpEdge(z+1)
            y = h_Data.GetBinContent(z+1)
            print '[%f,%f] : %f +- %f'%(x_l,x_r,y,h_Data.GetBinError(z+1))

        ## Loop over samples
        ## For Legend, save 
        HistsForLegend = []
        AliasForLegend = [] ## Prevent double-counting
        stack_Bkgd = ROOT.THStack("stack_Bkgd", "")
        h_Bkgd = 0
        ## Save systematic
        SystematicUps = dict()
        SystematicDowns = dict()
        ## If we take errors fro shapes
        h_TotalBackgroundFromShape = 0

        ## Loop over systematics
        for Syst in self.Systematics:

          if self.DoDebug:
            print '[DEBUG] Trying to make histogram for Syst = ',
            Syst.Print()

          h_Bkgd_ForSyst = 0
          dirName = Region.Name

          if Syst.Name!="Central":

            if Syst.Direction>0:
              dirName = "Syst_"+Syst.Name+"Up_"+Region.Name
            else:
              dirName = "Syst_"+Syst.Name+"Down_"+Region.Name

          for SampleGroup in self.SampleGroups:
            Color = SampleGroup.Color
            LegendAdded = False

            if self.DoDebug:
	      print "[DEBUG] SampleGroup.Samples: ", SampleGroup.Samples
            for Sample in SampleGroup.Samples:
              if self.DoDebug:
	        print "[DEBUG] Sample: ", Sample
                print '[DEBUG] Trying to make histogram for Sample = '+Sample

              f_Sample = ROOT.TFile(Indir+'/'+self.Filename_prefix+self.Filename_skim+self.Filename_suffix+Sample+'_'+Signal+'.root')
              h_Sample = 0

              ## Uncorrelated sources has Syst.Year = 2016 or 2017 or 2018
              ## For this cases, SampleGroup.Year should be matched
	      if (Syst.Name != 'Central'):
		print "Histo name: ", Sample+'_'+Syst.Name+str(Syst.Direction)
		h_Sample = f_Sample.Get(Sample+'_'+Syst.Name+Syst.Direction)
              ## 2) Lumi, MC normalizaion
              ## Use central and scale them later
              elif (Syst.Name in ["Lumi", "DYNorm", "NonPromptNorm", "OthersNorm"]):
                tmp_dirName = Region.Name
                h_Sample = f_Sample.Get(tmp_dirName+'/'+Variable.Name+'_'+tmp_dirName)
              ## For all other cases
              else:
		print "Histo name: ", self.Filename_prefix+Region.Name+'/'+self.Filename_prefix+Region.Name 
                h_Sample = f_Sample.Get(self.Filename_prefix+Region.Name+'/'+self.Filename_prefix+Region.Name)

              if not h_Sample:
                print 'No hist : %s %s'%(Syst.Name,Sample)
                continue

              print "h_Sample.GetBinContent(1): ", h_Sample.GetBinContent(1)
	      print "h_Sample.GetXaxis.GetNbins(): ", h_Sample.GetXaxis().GetNbins()
	      print "h_Sample.GetXaxis().GetBinUpEdge(15): ", h_Sample.GetXaxis().GetBinUpEdge(15)
              ## Make overflow
              h_Sample.GetXaxis().SetRangeUser(xMin,xMax)
              #h_Sample = mylib.MakeOverflowBin(h_Sample)

              h_Sample = self.Rebin(h_Sample, Region.Name, Variable.Name, nRebin)
              h_Sample.SetLineColor(Color)
              h_Sample.SetLineWidth(1)
              h_Sample.SetFillColor(Color)

              ## Scale
              MCSF, MCSFerr = 1., 0.
              if self.ScaleMC:
                ## now, only for DY
                if "DYJets" in Sample:
                  MCSF, MCSFerr = mylib.GetDYNormSF(SampleGroup.Year, Region.Name)
              h_Sample.Scale( MCSF )

              ## Manual systematic
              ## 1) [Lumi] Uncorrelated
              if (Syst.Name=="Lumi") and (Syst.Year==SampleGroup.Year):
                lumierr = mylib.LumiError(SampleGroup.Year)
                for ix in range(0,h_Sample.GetXaxis().GetNbins()):
                  y = h_Sample.GetBinContent(ix+1)
                  y_new = y + y*float(Syst.Direction)*lumierr
                  h_Sample.SetBinContent(ix+1, y_new)
              ## 2) [DYNorm] Uncorrelated, only for DY
              if (Syst.Name=="DYNorm") and ("DYJets" in Sample) and (Syst.Year==SampleGroup.Year):
                for ix in range(0,h_Sample.GetXaxis().GetNbins()):
                  y = h_Sample.GetBinContent(ix+1) ## h_Sample is already scaled by MCSF above
                  y_new = y * ( MCSF + float(Syst.Direction)*MCSFerr ) / MCSF
                  h_Sample.SetBinContent(ix+1, y_new)

              ## AddError option
              AddErrorOption = ''
              if self.AddErrorLinear:
                AddErrorOption = 'L'

              ## If central, add to h_Bkgd
              if self.DoDebug:
	        print "[DEBUG] Syst.Name: ", Syst.Name
	        print "[DEBUG] Sample: ", Sample
	      print "[DEBUG] Syst.Name: ", Syst.Name
	      print "[DEBUG] Sample: ", Sample
		
              if Syst.Name=="Central" and Sample=='total_background':
                if not h_TotalBackgroundFromShape:
                  h_TotalBackgroundFromShape = h_Sample.Clone()
                else:
                  h_TotalBackgroundFromShape = mylib.AddHistograms( h_TotalBackgroundFromShape, h_Sample, AddErrorOption )
                HistsToDraw[Sample] = h_Sample.Clone()

              elif Syst.Name=="Central":

                print "ADDING Sample:", Sample
                stack_Bkgd.Add( h_Sample )
                if not h_Bkgd:
                  print "Cloning sample"
                  h_Bkgd = h_Sample.Clone()
                else:
                  print "Adding sample"
                  print "h_Sample.GetBinContent(4): ", h_Sample.GetBinContent(4)
                  h_Bkgd = mylib.AddHistograms( h_Bkgd, h_Sample, AddErrorOption)

                print "h_Bkgd.GetBinContent(4): ", h_Bkgd.GetBinContent(4)
                HistsToDraw[Sample] = h_Sample.Clone()
                if (not LegendAdded) and (SampleGroup.TLatexAlias not in AliasForLegend):
                  HistsForLegend.append( [HistsToDraw[Sample],SampleGroup.TLatexAlias] )
                  AliasForLegend.append(SampleGroup.TLatexAlias)
                  LegendAdded = True
              ## else (i.e., systematic), add to h_Bkgd_ForSyst
              else:

                if not h_Bkgd_ForSyst:
                  h_Bkgd_ForSyst = h_Sample.Clone()
                else:
                  h_Bkgd_ForSyst = mylib.AddHistograms(h_Bkgd_ForSyst, h_Sample, AddErrorOption)

              ## Close file
              f_Sample.Close()

            ##==>End Sample loop
          ##==>End SampleGroup loop

          if Syst.Name!="Central":
            if Syst.Direction=="Up":
              SystematicUps[Syst.FullName()] = h_Bkgd_ForSyst
            else:
              SystematicDowns[Syst.FullName()] = h_Bkgd_ForSyst

        ##==>End Syst loop
        #print SystematicUps
        #print SystematicDowns

        ## Syst Up/Down . Max/Min
        h_Bkgd_TotErr_Max = h_Bkgd.Clone()
        h_Bkgd_TotErr_Min = h_Bkgd.Clone()
        for SystKey in SystematicUps:
          h_Up = SystematicUps[SystKey]
          h_Down = SystematicDowns[SystKey]
	  print "SystKey: ", SystKey
          for ix in range(0,h_Bkgd.GetXaxis().GetNbins()):
            x_l = h_Bkgd.GetXaxis().GetBinLowEdge(ix+1)
            x_r = h_Bkgd.GetXaxis().GetBinUpEdge(ix+1)
            y_Central = h_Bkgd.GetBinContent(ix+1)
            y_Up = h_Up.GetBinContent(ix+1)
            y_Down = h_Down.GetBinContent(ix+1)
            ## -.
            y_Max = max( max(y_Central,y_Up), y_Down)
            y_Min = min( min(y_Central,y_Up), y_Down)
            #print '[%d,%d] : %f, (Max,Min) = (%f,%f)'%(x_l,x_r,y_Central,y_Up,y_Down)

            ## Update Max
            err_Max_Current = h_Bkgd_TotErr_Max.GetBinError(ix+1)
            err_Max_ToAdd = y_Max-y_Central
            err_Max_Updated = math.sqrt( err_Max_Current*err_Max_Current + err_Max_ToAdd*err_Max_ToAdd )
            h_Bkgd_TotErr_Max.SetBinError(ix+1, err_Max_Updated)
            ## Update Min
            err_Min_Current = h_Bkgd_TotErr_Min.GetBinError(ix+1)
            err_Min_ToAdd = y_Central-y_Min
            err_Min_Updated = math.sqrt( err_Min_Current*err_Min_Current + err_Min_ToAdd*err_Min_ToAdd )
            h_Bkgd_TotErr_Min.SetBinError(ix+1, err_Min_Updated)
        ##==>End Systematic loop

        ## Blind mode
        if not Region.UnblindData:
          print "h_Bkgd.GetBinContent(1): ", h_Bkgd.GetBinContent(1)
          h_Data = h_Bkgd.Clone(h_Data.GetName())
          print "h_Data.GetNbinsX(): ", h_Data.GetNbinsX()
          for ibin in range(1,h_Data.GetNbinsX()+1):
              print "h_Data.GetBinError(ibin): ", h_Data.GetBinError(ibin)
              h_Data.SetBinError(ibin,h_Bkgd_TotErr_Max.GetBinError(ibin))
              print "h_Data.GetBinError(ibin): ", h_Data.GetBinError(ibin)
          print "h_Data.GetBinContent(1): ", h_Data.GetBinContent(1)
          h_Data.SetMarkerStyle(20)
          h_Data.SetMarkerSize(1.2)
          h_Data.SetMarkerColor(ROOT.kBlack)
          h_Data.SetLineColor(ROOT.kBlack)
          h_Data.SetLineWidth(2)

        ## hist => gr
	print "h_Bkgd_TotErr_Min.GetBinError(1): ", h_Bkgd_TotErr_Min.GetBinError(1)
	print "h_Bkgd_TotErr_Max.GetBinError(1): ", h_Bkgd_TotErr_Max.GetBinError(1)
        print "h_Data.GetBinContent(1): ", h_Data.GetBinContent(1)
        gr_Bkgd_TotErr = mylib.GetAsymmError(h_Bkgd_TotErr_Max,h_Bkgd_TotErr_Min)
        gr_Data = ROOT.TGraphAsymmErrors(h_Data)

        ## Error from shape
        if self.ErrorFromShape:
          if self.DoDebug:
            print 'Total background :'
            for z in range(0,h_TotalBackgroundFromShape.GetXaxis().GetNbins()):
              x_l = h_TotalBackgroundFromShape.GetXaxis().GetBinLowEdge(z+1)
              x_r = h_TotalBackgroundFromShape.GetXaxis().GetBinUpEdge(z+1)
              y = h_TotalBackgroundFromShape.GetBinContent(z+1)
              print '[%f,%f] : %f +- %f'%(x_l,x_r,y,h_TotalBackgroundFromShape.GetBinError(z+1))
	  #print "h_TotalBackgroundFromShape.GetBinError(1): ", h_TotalBackgroundFromShape.GetBinError(1)
          gr_Bkgd_TotErr = ROOT.TGraphAsymmErrors(h_TotalBackgroundFromShape)

        err_up_tmp = []
        err_down_tmp = []
        alpha = 1. - 0.6827
        for i in range(0, gr_Data.GetN()):
          N = gr_Data.GetY()[i]

          L = 0.                                          if (N==0.) else (ROOT.Math.gamma_quantile(alpha/2.,N,1.))
          U = ( ROOT.Math.gamma_quantile_c(alpha,N+1,1) ) if (N==0.) else (ROOT.Math.gamma_quantile_c(alpha/2.,N+1.,1.))

          #print '%d - %f + %f'%(N, N-L, U-N)

          if N!=0:
            gr_Data.SetPointEYlow(i, N-L )
            gr_Data.SetPointEYhigh(i, U-N )
            err_down_tmp.append(N-L)
            err_up_tmp.append(U-N)
          else:
            zerodata_err_low = 0.1
            zerodata_err_high = 1.8

            xlow = gr_Data.GetX()[i]-gr_Data.GetEXlow()[i]
            xhigh = gr_Data.GetX()[i]+gr_Data.GetEXhigh()[i]
            if self.ZeroDataCheckCut(Variable.Name,xlow,xhigh):
              zerodata_err_low = 0.
              zerodata_err_high = 0.

            gr_Data.SetPointEYlow(i, zerodata_err_low)
            gr_Data.SetPointEYhigh(i, zerodata_err_high)
            err_down_tmp.append(zerodata_err_low)
            err_up_tmp.append(zerodata_err_high)

        ## Legend
        lg = 0
        ## No signal
        if len(self.SignalsToDraw)==0:
          lg = ROOT.TLegend(0.55, 0.45, 0.92, 0.90)
        ## With Signal
        else:
          if Region.DrawRatio:
            lg = ROOT.TLegend(0.55, 0.56, 0.92, 0.90)
          else:
            lg = ROOT.TLegend(0.55, 0.56, 0.92, 0.90)
        lg.SetBorderSize(0)
        lg.SetFillStyle(0)

        if not self.NoErrorBand:
          lg.AddEntry(gr_Bkgd_TotErr, "Stat. #oplus syst. uncert.", "f")
        ## dummy graph for legend..
        ## this is because h_Data does not have horizontal error bars,
        ## and gr_data does not have points
        gr_Data_dummy = ROOT.TGraphAsymmErrors(gr_Data)
        gr_Data_dummy.SetMarkerStyle(20)
        gr_Data_dummy.SetMarkerSize(1.2)
        dataLegendGOption="ep"
        if Variable.Name=="WRCand_Mass":
          dataLegendGOption="lpe"
        if Region.DrawData:
          if Region.UnblindData:
            lg.AddEntry(gr_Data_dummy, "Data", dataLegendGOption)
          else:
            lg.AddEntry(gr_Data_dummy, "Total background", dataLegendGOption)
        for i_lg in range(0,len(HistsForLegend)):
          h_lg = HistsForLegend[ len(HistsForLegend)-1-i_lg ][0]
          tlatexaliax = HistsForLegend[ len(HistsForLegend)-1-i_lg ][1]
          print "tlatexaliax: ", tlatexaliax
          lg.AddEntry(h_lg,tlatexaliax,"f")

        ## Prepare canvas

        if self.DoDebug:
          print '[DEBUG] Preparing canvas..'

        c1_up_right.cd()
        if Region.Logy>0:
          c1_up_right.SetLogy(True)

        c1.cd()

        #### axis histograms

        h_dummy_up = ROOT.TH1D('h_dummy_up', '', nBin, xBins)
        h_dummy_up.GetXaxis().SetRangeUser(xMin, xMax)
        if nRebin>0:
          binsize = h_dummy_up.GetXaxis().GetBinUpEdge(1)-h_dummy_up.GetXaxis().GetBinLowEdge(1)
          str_binsize = '%d'%(binsize)
          if binsize!=int(binsize):
            str_binsize = '%1.2f'%(binsize)
          h_dummy_up.GetYaxis().SetTitle('Events / '+str_binsize+' '+Variable.Unit)
        else:
          h_dummy_up.GetYaxis().SetTitle('Events / bin')
        if Variable.Name=='WRCand_Mass':
          h_dummy_up.GetYaxis().SetTitle('Events / bin')

        h_dummy_down = ROOT.TH1D('h_dummy_down', '', nBin, xBins)
        h_dummy_down.GetYaxis().SetRangeUser(0.,2.0)

        if ('completeDisappearance' in Region.Name):
          h_dummy_down.GetYaxis().SetRangeUser(0.40,1.90)

        h_dummy_down.SetNdivisions(504,"Y")
        h_dummy_down.SetNdivisions(0,"X")
        #if ('completeDisappearance' in Region.Name):
          #h_dummy_down.GetXaxis().SetNdivisions(1)
          #h_dummy_down.GetXaxis().SetLabelOffset(100)
          #h_dummy_down.GetXaxis().CenterTitle(True)
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! xMin: ", xMin
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! xMax: ", xMax

        h_dummy_down.GetXaxis().SetRangeUser(xMin, xMax)
        h_dummy_down.GetXaxis().SetTitle(xtitle)
        h_dummy_down.GetXaxis().SetTitleSize(0.1)        
        h_dummy_down.GetXaxis().SetTitleOffset(0.1)        
        h_dummy_down.GetYaxis().SetTitle("#frac{Data}{Pred.}")
        h_dummy_down.SetFillColor(0)
        h_dummy_down.SetMarkerSize(0)
        h_dummy_down.SetMarkerStyle(0)
        h_dummy_down.SetLineColor(ROOT.kWhite)
        h_dummy_down.GetXaxis().SetLabelOffset(100)
        h_dummy_down.GetXaxis().CenterTitle(True)        

        print "Region.DrawRatio: ", Region.DrawRatio
        if Region.DrawRatio:
          h_dummy_up, h_dummy_down = canvas_margin.hist_axis(h_dummy_up, h_dummy_down)
        else:
          h_dummy_up.SetTitle("")
          h_dummy_up.GetYaxis().SetLabelSize(0.04)
          h_dummy_up.GetYaxis().SetTitleSize(0.054)
          h_dummy_up.GetYaxis().SetTitleOffset(1.30)
          h_dummy_up.GetXaxis().SetLabelSize(0.035)
          h_dummy_up.GetXaxis().SetTitleSize(0.055)
          h_dummy_up.GetXaxis().SetTitle(xtitle)

        h_dummy_up.GetYaxis().SetLabelOffset(100)
        h_dummy_up.GetYaxis().SetTitleOffset(100)
        h_dummy_up.GetYaxis().SetTickLength(0.09)
        h_dummy_up.GetXaxis().SetTickLength(0)
        h_dummy_up.GetXaxis().CenterTitle(True)        
        h_dummy_up.GetXaxis().SetTitle("Partial Dis. Region")        
        h_dummy_up.GetXaxis().SetTitleOffset(0.27) 
        h_dummy_up.GetXaxis().SetTitleSize(0.105)
        h_dummy_down.GetYaxis().SetLabelOffset(100)
        h_dummy_down.GetYaxis().SetTitleOffset(100)
        h_dummy_down.GetYaxis().SetTickSize(0.09)


        ## Get Y maximum
        yMax = max( yMax, mylib.GetMaximum(gr_Data) )
        yMax = max( yMax, mylib.GetMaximum(gr_Bkgd_TotErr) )
        ## Yaxis range
        yMin = 0.
        yMaxScale = 1.2
        if Region.Logy>0:
          yMaxScale = 10
          yMin = 0.2
        #print "!!!! yMin: ", yMin
        #h_dummy_up.GetYaxis().SetRangeUser( yMin, yMaxScale*yMax )
        h_dummy_up.GetYaxis().SetRangeUser( yMin, 1000 )

        ## Exception control

        if (Variable.Name=="staDEoverE"):
          if ("partialDisappearance" in Region.Name):
            h_dummy_up.GetYaxis().SetRangeUser( 1, 10000*yMax )

        if self.DoDebug:
          print '[DEBUG] Canvas is ready'

        ## Draw up
        c1_up_right.cd()

        print "[DEBUG] h_dummy_up.Integral(): ", h_dummy_up.Integral()
        h_dummy_up.Draw("histsame")
        stack_Bkgd.Draw("histsame")

        gr_Bkgd_TotErr.SetMarkerColor(0)
        gr_Bkgd_TotErr.SetMarkerSize(0)
        gr_Bkgd_TotErr.SetFillStyle(3013)
        gr_Bkgd_TotErr.SetFillColor(ROOT.kBlack)
        gr_Bkgd_TotErr.SetLineColor(0)
	print "self.NoErrorBand: ", self.NoErrorBand
        if not self.NoErrorBand:
          gr_Bkgd_TotErr.Draw("sameE2")

        gr_Data.SetLineWidth(2)
        gr_Data.SetMarkerSize(0.)
        gr_Data.SetMarkerColor(ROOT.kBlack)
        gr_Data.SetLineColor(ROOT.kBlack)
        if Region.DrawData:
	  print "!!!!!!!!!!!!!!!!! h_Data.Integral(): ", h_Data.Integral()
          h_Data.Draw("p EX0 histsame")
          #gr_Data.Draw("p0same")


        ## Signal
        h_Sigs = []
        print "self.SignalsToDraw: ", self.SignalsToDraw
        for Sig in self.SignalsToDraw:
	  print "Sig.mDBrem_FileString: ", Sig.mDBrem_FileString
	  #if Signal not in Sig.mDBrem_FileString: continue
	  print "Sig.mDBrem_FileString: ", Sig.mDBrem_FileString

          #fpullpath_Sig = Indir+'/'+self.Filename_prefix+self.Filename_skim+self.Filename_suffix+'DBrem_'+Sig.mDBrem_FileString+'_'+Signal+'.root'
          fpullpath_Sig = Indir+'/'+self.Filename_prefix+self.Filename_skim+self.Filename_suffix+'DBrem_'+Sig.mDBrem_FileString+'_'+Sig.mDBrem_FileString+'.root'
          f_Sig = ROOT.TFile(fpullpath_Sig)
          print fpullpath_Sig
          print self.Filename_prefix+Region.Name+'/DBrem_'+Region.Name
          h_Sig = f_Sig.Get(self.Filename_prefix+Region.Name+'/DBrem_'+Region.Name)
          if not h_Sig:
            print fpullpath_Sig
            print self.Filename_prefix+Region.Name+'/DBrem_'+Region.Name
            continue

          print "h_Sig.GetBinContent(1): ", h_Sig.GetBinContent(1)
          ## Make overflow
          h_Sig.GetXaxis().SetRangeUser(xMin,xMax)
          #h_Sig = mylib.MakeOverflowBin(h_Sig)

          ## Rebin
          h_Sig = self.Rebin(h_Sig, Region.Name, Variable.Name, nRebin)
          print "2nd check h_Sig.GetBinContent(1): ", h_Sig.GetBinContent(1)
          print "h_Sig.Integral(): ", h_Sig.Integral()

          ## Scale
          h_Sig.Scale( Sig.xsec * Sig.xsecScale )
          print "3rd check h_Sig.GetBinContent(1): ", h_Sig.GetBinContent(1)

          ## Att
          h_Sig.SetLineWidth(3)
          h_Sig.SetLineColor(Sig.Color)
          h_Sig.SetLineStyle(Sig.Style)

          ## Save
          h_Sigs.append(h_Sig)

          ## legend
          lg.AddEntry(h_Sig, Sig.GetTLatexAlias(), 'l')
          lg.SetTextFont(42)

          ## Draw
	  print "!!!!!!!!!!!!!!!!! h_Sig.Integral(): ", h_Sig.Integral()
          h_Sig.Draw("histsame")

        h_dummy_up.Draw("axissame")

        ## Legend
        #lg.Draw()

        ## Draw down
        c1_down_right.cd()
        h_dummy_down.Draw("axis hist")
        print "h_dummy_down.GetXaxis().GetTitle(): ", h_dummy_down.GetXaxis().GetTitle()


        ## values must be set later
        h_Data_Ratio = h_Data.Clone('h_Data_Ratio')
        ## BinContent set by divide here, but errors must be set later
        tmp_h_Data_Ratio = h_Data.Clone()
        tmp_h_Data_Ratio.Divide(h_Bkgd)
        gr_Data_Ratio = ROOT.TGraphAsymmErrors(tmp_h_Data_Ratio)
        gr_Data_Ratio.SetName('gr_Data_Ratio')
        gr_Data_Ratio.SetLineWidth(2)
        gr_Data_Ratio.SetMarkerSize(0.)
        gr_Data_Ratio.SetLineColor(ROOT.kBlack)
        ## values must be set later, but BinContent will be simply 1
        gr_Bkgd_Ratio = gr_Bkgd_TotErr.Clone('gr_Bkgd_Ratio')

        for i in range(1,h_Data_Ratio.GetXaxis().GetNbins()+1):

          ## FIXME for zero? how?
          if h_Bkgd.GetBinContent(i)!=0:

            ## ratio point
            ## BinContent = Data/Bkgd
            ## BinError = DataError/Bkgd
            h_Data_Ratio.SetBinContent( i, h_Data_Ratio.GetBinContent(i) / h_Bkgd.GetBinContent(i) )
            h_Data_Ratio.SetBinError ( i, h_Data_Ratio.GetBinError(i) / h_Bkgd.GetBinContent(i) )

            if err_down_tmp[i-1]!=0.:
              gr_Data_Ratio.SetPointEYlow(i-1, err_down_tmp[i-1] / h_Bkgd.GetBinContent(i) )
              gr_Data_Ratio.SetPointEYhigh(i-1, err_up_tmp[i-1] / h_Bkgd.GetBinContent(i))
              if Variable.Name!="WRCand_Mass":
                gr_Data_Ratio.SetPointEXlow(i-1, 0)
                gr_Data_Ratio.SetPointEXhigh(i-1, 0)
            else:
              gr_Data_Ratio.SetPointEYlow(i-1, 0)
              gr_Data_Ratio.SetPointEYhigh(i-1, 1.8 / h_Bkgd.GetBinContent(i))
              if Variable.Name!="WRCand_Mass":
                gr_Data_Ratio.SetPointEXlow(i-1, 0)
                gr_Data_Ratio.SetPointEXhigh(i-1, 0)

            ## ratio allerr
            ## BinContent = 1
            ## BinError = Bkgd(Stat+Syst)Error/Bkgd

            gr_Bkgd_Ratio.SetPoint(i-1,h_Bkgd.GetXaxis().GetBinCenter(i), 1.)
            gr_Bkgd_Ratio.SetPointEYhigh( i-1, gr_Bkgd_Ratio.GetErrorYhigh(i-1) / h_Bkgd.GetBinContent(i) )
            gr_Bkgd_Ratio.SetPointEYlow( i-1,  gr_Bkgd_Ratio.GetErrorYlow(i-1) / h_Bkgd.GetBinContent(i) )

          elif h_Bkgd.GetBinContent(i)==0. and h_Data_Ratio.GetBinContent(i)==0.:

            h_Data_Ratio.SetBinContent( i, 0 )
            h_Data_Ratio.SetBinError ( i, 0 )
            gr_Data_Ratio.SetPoint(i-1, 0, 0)
            gr_Data_Ratio.SetPointEYlow(i-1, 0)
            gr_Data_Ratio.SetPointEYhigh(i-1, 0)

            gr_Bkgd_Ratio.SetPoint(i-1,h_Bkgd.GetXaxis().GetBinCenter(i), 1.)
            gr_Bkgd_Ratio.SetPointEYhigh( i-1, 0. )
            gr_Bkgd_Ratio.SetPointEYlow( i-1, 0. )

            if Variable.Name!="WRCand_Mass":
              gr_Data_Ratio.SetPointEXlow(i-1, 0)
              gr_Data_Ratio.SetPointEXhigh(i-1, 0)

          ## If bkgd <= 0
          else:
            this_max_ratio = 20.0
            this_data = h_Data_Ratio.GetBinContent(i)
            this_data_err = h_Data_Ratio.GetBinError(i)

            h_Data_Ratio.SetBinContent( i, this_max_ratio )
            h_Data_Ratio.SetBinError ( i, this_data_err*this_max_ratio/this_data )

            tmp_x = ROOT.Double(0.)
            tmp_y = ROOT.Double(0.)
            gr_Data_Ratio.GetPoint(i-1, tmp_x, tmp_y)
            gr_Data_Ratio.SetPoint(i-1, tmp_x, this_max_ratio)
            gr_Data_Ratio.SetPointEYlow(i-1, err_down_tmp[i-1]*this_max_ratio/this_data)
            gr_Data_Ratio.SetPointEYhigh(i-1, err_up_tmp[i-1]*this_max_ratio/this_data)

            gr_Bkgd_Ratio.SetPoint(i-1,h_Bkgd.GetXaxis().GetBinCenter(i), 1.)
            gr_Bkgd_Ratio.SetPointEYhigh( i-1, 0. )
            gr_Bkgd_Ratio.SetPointEYlow( i-1, 0. )

            if Variable.Name!="WRCand_Mass":
              gr_Data_Ratio.SetPointEXlow(i-1, 0)
              gr_Data_Ratio.SetPointEXhigh(i-1, 0)
        ##==>End bin loop

        gr_Bkgd_Ratio.SetMarkerColor(0)
        gr_Bkgd_Ratio.SetMarkerSize(0)
        gr_Bkgd_Ratio.SetFillStyle(3013)
        gr_Bkgd_Ratio.SetFillColor(ROOT.kBlack)
        gr_Bkgd_Ratio.SetLineColor(0)
        gr_Bkgd_Ratio.Draw("sameE2")

        h_Data_Ratio.Draw("p EX0 same")
        #gr_Data_Ratio.Draw("p0same")

        ## y=1 graph
        g1_x = [-9000, 9000]
        g1_y = [1, 1]
        g1 = ROOT.TGraph(2, array("d", g1_x ), array("d", g1_y ))
        g1.Draw("same")

        c1.cd()
        c1.SaveAs(Outdir+'Test_partial.pdf')

        ## TLatex
        c1_up_right.cd()
        channelname = ROOT.TLatex()
        channelname.SetNDC()
        channelname.SetTextSize(0.03)
        #channelname.DrawLatex(0.25, 0.88, "#font[42]{"+Region.TLatexAlias+"}")

        ## Extra lines
        #exec(self.ExtraLines)
        print "Draw self.ExtraLines"


    print "All variables are done for this region, closing data TFile"
    f_Data.Close()
    

        ## Save
   print "SAVING c1"
   c1.cd()
   latex_CMS = ROOT.TLatex()
   latex_CMSPriliminary = ROOT.TLatex()
   latex_Lumi = ROOT.TLatex()

   latex_CMS.SetNDC()
   latex_CMSPriliminary.SetNDC()
   latex_Lumi.SetNDC()
   latex_CMS.SetTextSize(0.055)
   latex_CMS.DrawLatex(0.13, 0.89, "#font[62]{CMS}")
   latex_CMSPriliminary.SetTextSize(0.055)
   latex_CMSPriliminary.DrawLatex(0.13, 0.85, "#font[42]{#it{#scale[0.8]{Preliminary}}}")
   #latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS}")

   latex_Lumi.SetTextSize(0.035)
   latex_Lumi.SetTextFont(42)
   latex_Lumi.DrawLatex(0.83, 0.96, mylib.TotalLumi(float(self.DataYear))+" fb^{-1} (13 TeV)")
   
   c1.SaveAs(Outdir+Variable.Name+'_'+Region.Name+'.pdf')
   print Variable.Name+'_'+Region.Name+'.pdf ==> Saved'

   c1.Close()

      ##==>End Variable loop

#      if self.DoDebug:
#        print '[DEBUG] All variables are done for this region, closing data TFile'
#      f_Data.Close()

   ##==>End Region loop
   if self.DoDebug:
     print '[DEBUG] All regions are done.'
