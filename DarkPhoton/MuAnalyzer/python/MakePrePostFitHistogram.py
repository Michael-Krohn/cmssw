import os,ROOT
from array import array
import numpy as np

WORKING_DIR = os.environ['PLOTTER_WORKING_DIR']
dataset = os.environ['CATANVERSION']
FILE_PATH = os.environ['FILE_PATH']
PLOT_PATH = os.environ['PLOT_PATH']

Signals = [
#"2p0",
"1p0",
"0p4",
"0p2",
]

PreORPosts = [
"Pre",
"Post",
]

for Signal in Signals:

  for PreORPost in PreORPosts:
  
    basedir = FILE_PATH+'/'+dataset+'/FitDiagnostics/'
    outdirBase = FILE_PATH+'/'+dataset+'/'+PreORPost+'Fit/'
  
    Years = [
    "2018",
    ]
  
    Regions = [
    "completeDisappearance",
    #"partialDisappearance",
    #"partialDisappearanceLargeStaDPhi",
    "partialDisappearanceBDT",
    ]
  
    Samples = [
    'TTBar',
    'DY',
    'muPX',
  #  'WJets',
    'total_background',
    'DBrem_'+Signal,
    ]
  
    for Sample in Samples:
  
  
        print '-------------------------------------------'
        print 'Making '+Sample+' for '+'...'
        this_filename = 'DBrem_'+Sample+'_'+Signal+'.root'
        print '- filename = '+this_filename
  
        this_outdirBase = outdirBase+'/'
        os.system('mkdir -p '+this_outdirBase)
        print this_outdirBase
  
        out = ROOT.TFile(this_outdirBase+this_filename,'RECREATE')
  
  
        shapedirName = 'shapes_fit_b'
        if 'Pre' in PreORPost or 'DBrem' in Sample:
            shapedirName = 'shapes_prefit'
  
        for Region in Regions:
	    if "completeDisappearance" in Region and "TTBar" in Sample: continue
	    #if "partialDisappearance" in Region and "muPX" in Sample: continue
	    if "partialDisappearance" in Region and "TTBar" in Sample: continue
  
            shapefilename = 'fitDiagnostics_DBrem_'+Signal+'.root'
            f = ROOT.TFile(basedir+'/'+shapefilename)
  
            dir_shapes_fit_b = f.Get(shapedirName)
  
            print '  - Region = '+Region
  
            #plotbins = [0, 1]
            #plotbins = [80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]
            plotbins = [80,82,84,86,88,90,92,94,96,98,100]
            if "LargeStaDPhi" in Region:
#              plotbins = np.linspace(0,20,51)
	       plotbins = [0,1,2,2.5,3.,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.5,10,20]
	    elif "BDT" in Region:
               #plotbins = [0,0.05,0.1,0.2,0.3,0.4,0.55,0.7,1.0]
	       #plotbins = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
	       #plotbins = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.84,0.88,0.9,0.92,0.94,0.96,0.98,1.0]
	       plotbins = [0.98,1.0]

#	       plotbins = [0,0.05,0.1,1.0]
	       print "!!!!!!!!!!!!!Signal: ", Signal
#	       if "2p0" in Signal:
#        	 plotbins = [0,0.05,0.1,1.0]
#      	       elif "1p0" in Signal:
 #       	 plotbins = [0,0.05,0.1,0.2,0.3,0.4,0.55,0.7,1.0]
 #     	       elif "0p8" in Signal:
 #       	 plotbins = [0,0.05,0.1,0.2,0.3,0.4,0.5,1.0]
 #     	       elif "0p6" in Signal:
 #       	 plotbins = [0,0.05,0.1,0.2,0.3,0.45,0.6,1.0]
 #     	       elif "0p4" in Signal:
 #       	 plotbins = [0,0.05,0.1,0.2,0.3,0.4,0.55,1.0]
 #     	       elif "0p2" in Signal:
 #       	 plotbins = [0,0.05,0.1,0.2,0.3,0.4,0.55,0.7,1.0]

            elif "partial" in Region:
              plotbins = np.linspace(-1,1,51)
  
            dirName = 'DBrem_'+Region
            #### Exception for CRs
  
            print '    - bin =',
            print plotbins
  
            thisShapeDirName = Region
            print thisShapeDirName
            dir_fitDiag_thisShapes = dir_shapes_fit_b.Get(thisShapeDirName)
            ## hist_binned starts from 800 GeV
            hist_binned = dir_fitDiag_thisShapes.Get(Sample)

            hist_massbinned = ROOT.TH1D(dirName, '', len(plotbins)-1, array('d', plotbins))
            for i in range(0, hist_massbinned.GetXaxis().GetNbins()):
              binWidth = plotbins[i+1] - plotbins[i]
              ix = i+1
              #hist_massbinned.SetBinContent( ix, ( hist_binned.GetBinContent(ix) / binWidth ) * 0.05 )
              #hist_massbinned.SetBinError( ix, ( hist_binned.GetBinError(ix) / binWidth ) * 0.05 )
              hist_massbinned.SetBinContent( ix, ( hist_binned.GetBinContent(ix) ))
              hist_massbinned.SetBinError( ix, ( hist_binned.GetBinError(ix) ))
              #print hist_binned.GetBinContent(ix)
  
            out.mkdir(dirName)
            out.cd(dirName)
            hist_massbinned.Write()
            out.cd()
  
  
            #### Closing fitDiag file
            f.Close()
  
        #### Closing output file
        out.Close()
