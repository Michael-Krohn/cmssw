import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'runLocally',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if running locally")

options.parseArguments()

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

if(options.runLocally):
	process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
	import FWCore.Utilities.FileUtils as FileUtils
	mylist = FileUtils.loadListFromFile ('datafiles/DYJets.txt')
	readFiles = cms.untracked.vstring( *mylist)

process.options = cms.untracked.PSet(
 SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
 IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
 #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/10003/05215DB9-A1EB-D449-9332-2005BDA65284.root")
    fileNames = cms.untracked.vstring()
)

process.MCWeights = cms.EDAnalyzer('WeightsFinder',
    edmPileupInfo = cms.InputTag("mixData"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("MCWeights.root"),
				   closeFileFast = cms.untracked.bool(False)
)

process.p = cms.Path(process.MCWeights)
