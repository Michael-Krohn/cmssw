import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  False,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.bool,
				  "True if is MC dataset")

options.register( 'runLocally',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if running locally")

options.register( 'signalRegion',
                                  True,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if running on MC in isolated region")

options.parseArguments()


process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.isMC):
        process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
else:
        process.GlobalTag.globaltag = cms.string('106X_dataRun2_v32')

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
	   mylist = FileUtils.loadListFromFile ('datafiles/HardBrem_RECO/dataFiles_DYJetsToLL_2nd_0000_test.txt')
	else:
	   mylist = FileUtils.loadListFromFile('datafiles/dataFiles_2018D.txt')
	readFiles = cms.untracked.vstring( *mylist)

process.options = cms.untracked.PSet(
 SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
 IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
 #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring("file:RECO_RAW2DIGI_L1Reco_RECO_ALCA_EI_PAT_DQM_99.root")
    fileNames = cms.untracked.vstring("file:ConfFile_MuPXAnalyzer_HardBrem_Run2018D_RAWRECO_BadTrackIso_018.root")
#    fileNames = readFiles
)
if(options.isMC):
   process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
else:
   process.source.lumisToProcess = LumiList.LumiList(filename = "/local/cms/user/revering/dphoton/slc7/CMSSW_10_6_17_patch1/src/DarkPhoton/MuAnalyzer/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt").getVLuminosityBlockRange()

process.demo = cms.EDAnalyzer('MuPXAnalyzer',
    recoMuons = cms.InputTag("muons"),
    photons = cms.InputTag("photons"),
    tracks = cms.InputTag("generalTracks"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    edmPileupInfo = cms.InputTag("mixData"),
    CSCSegmentLabel = cms.InputTag("cscSegments"),
    trigResults = cms.InputTag("TriggerResults","","HLT"),
    muonPathsToPass = cms.vstring("HLT_IsoMu24_v","HLT_IsoMu27_v"),
    #HBHERecHits = cms.InputTag("hbhereco"),
    EERecHits = cms.InputTag("reducedEcalRecHitsEE"),
    EBRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    #HBHERecHits = cms.InputTag("reducedHcalRecHits","hbhereco"),
    HBHERecHits = cms.InputTag("hbhereco"),
    PFJets = cms.InputTag("ak4PFJetsCHS"),
    TowerSource = cms.InputTag("caloStage2Digis","CaloTower"),
    CaloJetSource = cms.InputTag("ak4CaloJets"),
    isMC = cms.untracked.bool(options.isMC),
    signalRegion = cms.untracked.bool(options.signalRegion)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)


process.p = cms.Path(process.demo)
