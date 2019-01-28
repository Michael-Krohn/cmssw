import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  True,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.bool,
				  "True if is MC dataset")

options.register( 'runRandomTrack',
                                  False,
                                  VarParsing.multiplicity.singleton,
				  VarParsing.varType.bool,
				  "True if is studying random track efficiency")

options.register( 'runLocally',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if running locally")

options.parseArguments()


process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.isMC):
        process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v10')
else:
        process.GlobalTag.globaltag = cms.string('94X_dataRun2_ReReco_EOY17_v6')

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

if(options.runLocally):
	process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

if(options.runLocally):
	import FWCore.Utilities.FileUtils as FileUtils
	mylist = FileUtils.loadListFromFile ('Filtered_Files_DY_temp.txt')
	readFiles = cms.untracked.vstring( *mylist)

process.options = cms.untracked.PSet(
 SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
 IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
 #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    fileNames = readFiles
)

process.demo = cms.EDAnalyzer('MuAnalyzer',
    recoMuons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    CSCSegmentLabel = cms.InputTag("cscSegments"),
    trigResults = cms.InputTag("TriggerResults","","HLT"),
    muonPathsToPass = cms.vstring("HLT_IsoMu24_v","HLT_IsoMu27_v"),
    HBHERecHits = cms.InputTag("reducedHcalRecHits","hbhereco"),
    isMC = cms.untracked.bool(options.isMC),
    runRandomTrackEfficiency = cms.untracked.bool(options.runRandomTrack)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)


process.p = cms.Path(process.demo)
