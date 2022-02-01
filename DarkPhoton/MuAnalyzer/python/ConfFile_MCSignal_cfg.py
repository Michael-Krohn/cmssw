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

options.register( 'scaleFactor',
                                  2000,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.float,
                                  "Reweighting factor for dark brem process")
				

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
	   #mylist = FileUtils.loadListFromFile ('datafiles/DYJetsToMM_dbrem.txt')
           mylist = FileUtils.loadListFromFile ('datafiles/DYJets_PrivateProduction_mumuFilter.txt')
           #mylist = FileUtils.loadListFromFile ('filetest.txt')
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
    electrons = cms.InputTag("gedGsfElectrons"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    edmPileupInfo = cms.InputTag("mixData"),
    g4SimHits = cms.InputTag("g4SimHits"),
    dbremweighttag = cms.InputTag("g4SimHits","DBremEventWeight"),
    dbremInitialETag = cms.InputTag("g4SimHits","DBremInitialEnergy"),
    dbremFinalETag = cms.InputTag("g4SimHits","DBremFinalEnergy"),
    dbremAngleTag = cms.InputTag("g4SimHits","DBremAngle"),
    dbremXTag = cms.InputTag("g4SimHits","DBremLocationX"),
    dbremYTag = cms.InputTag("g4SimHits","DBremLocationY"),
    dbremZTag = cms.InputTag("g4SimHits","DBremLocationZ"),
    CSCSegmentLabel = cms.InputTag("cscSegments"),
    trigResults = cms.InputTag("TriggerResults","","HLT"),
    muonPathsToPass = cms.vstring("HLT_IsoMu24_v","HLT_IsoMu27_v"),
    HBHERecHits = cms.InputTag("hbhereco"),
    PFJets = cms.InputTag("ak4PFJetsCHS"),
    CaloJetSource = cms.InputTag("ak4CaloJets"),
    DisplacedStandaloneTracks = cms.InputTag("displacedStandAloneMuons"),
    StandAloneTracks = cms.InputTag("standAloneMuons"),
    EERecHits = cms.InputTag("reducedEcalRecHitsEE"),
    EBRecHits = cms.InputTag("reducedEcalRecHitsEB"),
#    ReducedHBHERecHits = cms.InputTag("reducedHcalRecHits"),
    isMC = cms.untracked.bool(options.isMC),
    isSig = cms.untracked.bool(options.isSig),
    hasDpho = cms.untracked.bool(options.hasDpho),
    scaleFactor = cms.untracked.double(options.scaleFactor),
    runrandomTrackEfficiency = cms.untracked.bool(options.runRandomTrack)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)


process.p = cms.Path(process.demo)
