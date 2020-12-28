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

options.register( 'isSig',
				  True,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.bool,
				  "True if using signal injected events")

options.register( 'hasDpho',
                                  True,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if dark brem particle is present")

options.parseArguments()


process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.isMC):
        process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v17')
else:
        process.GlobalTag.globaltag = cms.string('94X_dataRun2_ReReco_EOY17_v6')

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

if(options.runLocally):
	process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
	import FWCore.Utilities.FileUtils as FileUtils
        #mylist = FileUtils.loadListFromFile ('file_temp_RECO.txt')
	#mylist = FileUtils.loadListFromFile ('Filtered_Files_SingleMuon2017E-PromptReco_temp.txt')
	#mylist = FileUtils.loadListFromFile ('file_temp.txt')
        #mylist = FileUtils.loadListFromFile ('Filtered_Files_DY_2017.txt')
	if(options.isMC):
	   mylist = FileUtils.loadListFromFile ('datafiles/map1p0_13TeV.txt')
           #mylist = FileUtils.loadListFromFile ('datafiles/DY_MC_Files.txt')
        else:
	   mylist = FileUtils.loadListFromFile('datafiles/Filtered_Files_Doublemuon.txt')
	readFiles = cms.untracked.vstring( *mylist)

process.options = cms.untracked.PSet(
 SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
 IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
 #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    #fileNames = cms.untracked.vstring("file:RECOdata_Test.root")
    fileNames = readFiles
)

process.demo = cms.EDAnalyzer('SigAnalyzer',
    recoMuons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    g4SimHits = cms.InputTag("g4SimHits"),
    CSCSegmentLabel = cms.InputTag("cscSegments"),
    trigResults = cms.InputTag("TriggerResults","","HLT"),
    muonPathsToPass = cms.vstring("HLT_IsoMu24_v","HLT_IsoMu27_v"),
    HBHERecHits = cms.InputTag("hbhereco"),
    StandAloneTracks = cms.InputTag("standAloneMuons"),
#    ReducedHBHERecHits = cms.InputTag("reducedHcalRecHits"),
    isMC = cms.untracked.bool(options.isMC),
    isSig = cms.untracked.bool(options.isSig),
    hasDpho = cms.untracked.bool(options.hasDpho),
    runrandomTrackEfficiency = cms.untracked.bool(options.runRandomTrack)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)


process.p = cms.Path(process.demo)
