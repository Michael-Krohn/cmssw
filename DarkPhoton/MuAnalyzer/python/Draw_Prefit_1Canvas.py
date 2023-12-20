import os,ROOT
from Plotter_1Canvas import SampleGroup, DBremSignalInfo, Variable, Region, Systematic
from Plotter_1Canvas import Plotter
import mylib
import argparse

ROOT.gROOT.SetBatch(ROOT.kTRUE)

## Arguments

parser = argparse.ArgumentParser(description='PreFit plot commands')
parser.add_argument('-c', dest='Category', type=int, default=0)
parser.add_argument('-y', dest='Year', type=int, default=2018)
parser.add_argument('--debug',action='store_true')
parser.add_argument('--blind',action='store_true')
args = parser.parse_args()

## Blind mode
UnblindData = not args.blind

## Enviroment
WORKING_DIR = os.environ['PLOTTER_WORKING_DIR']
dataset = os.environ['CATANVERSION']
ENV_PLOT_PATH = os.environ['PLOT_PATH']

## Signal mass
mDBrem = 1.0
DBremSignalInfoToDraw2nd = DBremSignalInfo(mDBrem = mDBrem)
DBremSignalInfoToDraw2nd.Color = ROOT.kBlack
DBremSignalInfoToDraw2nd.useOfficial = True
DBremSignalInfoToDraw2nd.Style = 3
DBremSignalInfoToDraw2nd.xsecScale = 10.0

mDBrem = 0.4
DBremSignalInfoToDraw3rd = DBremSignalInfo(mDBrem = mDBrem)
DBremSignalInfoToDraw3rd.Color = ROOT.kBlack
DBremSignalInfoToDraw3rd.useOfficial = True
DBremSignalInfoToDraw3rd.Style = 3
DBremSignalInfoToDraw3rd.xsecScale = 1.0

mDBrem = 0.2
DBremSignalInfoToDraw4th = DBremSignalInfo(mDBrem = mDBrem)
DBremSignalInfoToDraw4th.Color = ROOT.kBlack
DBremSignalInfoToDraw4th.useOfficial = True
#DBremSignalInfoToDraw4th.Style = 4
DBremSignalInfoToDraw4th.Style = 2
DBremSignalInfoToDraw4th.xsecScale = 1.0
#DBremSignalInfoToDraw.xsec = mylib.GetSignalXsec(WORKING_DIR+'/data/'+dataset+'/xsec_190705_GenXsecAN_eeANDmm.txt', mWR, mN)

m = Plotter()

m.DoDebug = args.debug

#### In/Out
m.DataYear = args.Year
str_Year = str(args.Year)
m.InputDirectory = WORKING_DIR+'/rootfiles/'+dataset+'/PreFit/'
if args.Year<0:
  str_Year = 'YearCombined'
m.DataDirectory = str_Year
m.Filename_prefix = "DBrem_"
m.Filename_suffix = ""
m.Filename_skim = ""
m.OutputDirectory = ENV_PLOT_PATH+"/"+dataset+"/PreFit_2022_InitialFit/"+str_Year+"/"


#### Systematic
#tmp_Systematics = ["Trig","ID","PU"]
tmp_Systematics = []
m.ErrorFromShape = True
m.AddErrorLinear = False
if args.Year<0:
  m.AddErrorLinear = True

m.Systematics = [ Systematic(Name="Central", Direction="Central", Year=-1) ]
for s in tmp_Systematics:
    m.Systematics.append( Systematic(Name=s, Direction="Up", Year=-1) )
    m.Systematics.append( Systematic(Name=s, Direction="Down", Year=-1) )
m.PrintSystematics()

#### Define Samples
from PredefinedSamples import *
m.SampleGroups = [
    SampleGroup_DY, SampleGroup_WJets, SampleGroup_TTBar, SampleGroup_muPX, SampleGroup_total_background
  ]
#### Signal
#m.SignalsToDraw = [DBremSignalInfoToDraw, DBremSignalInfoToDraw2nd, DBremSignalInfoToDraw3rd, DBremSignalInfoToDraw4th]
#m.SignalsToDraw = [DBremSignalInfoToDraw2nd, DBremSignalInfoToDraw3rd, DBremSignalInfoToDraw4th]
m.SignalsToDraw = [DBremSignalInfoToDraw2nd, DBremSignalInfoToDraw4th]

#### Print
m.PrintSamples()

#### Define reiongs
print "UnblindData: ", UnblindData
if args.Category==0:
  m.RegionsToDraw = [
    Region('completeDisappearance', UnblindData=UnblindData, Logy=1, TLatexAlias='Complete Disappearance'),
    #Region('partialDisappearance', UnblindData=UnblindData, Logy=1, TLatexAlias='Partial Disappearance'),
    #Region('partialDisappearanceLargeStaDPhi', UnblindData=UnblindData, Logy=1, TLatexAlias='Partial Disappearance'),
    Region('partialDisappearanceBDT', UnblindData=UnblindData, Logy=1, TLatexAlias='Partial Disappearance'),
  ]

m.PrintRegions()

#### Define Variables
m.VariablesToDraw = [
  Variable('nEvents', 'Invariant Mass of Tag and Probe (GeV)', 'nEvents'),
  #Variable('nEvents', 'Complete Disappearance Region', 'nEvents'),
  Variable('bdtScore', 'Part. Dis. Region', ''),
]
m.PrintVariables()

#### Extra lines
m.ExtraLines='''tl = ROOT.TLatex()
tl.SetNDC()
tl.SetTextSize(0.03)
tl.DrawLatex(0.2, 0.815, '#font[42]{Prefit}')
'''

#### Draw
m.Draw()
