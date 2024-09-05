import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  False,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.bool,
				  "True if is MC dataset")

options.register( 'sameSign',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if same sign control region")

options.register( 'offPeakCR',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Require inv. Mass <70 or >110")

options.register( 'ecalBrem',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Require brem to be in ECAL")

options.register( 'studyBarrel',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Allow tag aligned muons in the barrel")

options.register( 'tagAlignedCr',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Use tag track instead of probe")

options.register( 'ttbarCR',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if ttbar control region")

options.register( 'metFailName',
                                  "",
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.string,
                                  "Name of met filter results file")

options.register( 'missingCSC',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if missing CSC filtered events")

options.register( 'highHcalCr',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if high HCAL iso CR")

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
				
options.register( 'binSFsfname',
                                  "/local/cms/user/revering/datafiles/dbremBinningSFs/dbremBinningSFs_map0p2.root",
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.string,
                                  "Dark Photon Energy binning filename")

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
	   mylist = FileUtils.loadListFromFile ('datafiles/ZmmSim_EndcapMuFilter_part12.txt')
           #mylist = FileUtils.loadListFromFile ('datafiles/DYJets_PrivateProduction_mumuFilter.txt')
           #mylist = FileUtils.loadListFromFile ('datafiles/testfiles.txt')
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
    met = cms.InputTag("pfMet"),
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
    CSCRecHitsLabel = cms.InputTag("csc2DRecHits"),
    DTSegmentLabel = cms.InputTag("dt4DSegments"),
    cscSimHits = cms.InputTag("g4SimHits","MuonCSCHits"),
    RPCRecHitTag = cms.InputTag("rpcRecHits"),
    trigResults = cms.InputTag("TriggerResults","","HLT"),
    muonPathsToPass = cms.vstring("HLT_IsoMu24_v","HLT_IsoMu27_v"),
    HBHERecHits = cms.InputTag("hbhereco"),
    PFJets = cms.InputTag("ak4PFJetsCHS"),
    CaloJetSource = cms.InputTag("ak4CaloJets"),
    btagAlgos = cms.VInputTag(cms.InputTag("pfDeepCSVJetTags:probb"), cms.InputTag("pfDeepCSVJetTags:probbb")),
    DisplacedStandaloneTracks = cms.InputTag("displacedStandAloneMuons"),
    StandAloneTracks = cms.InputTag("standAloneMuons"),
    RefittedStandaloneTracks = cms.InputTag("refittedStandAloneMuons"),
    EERecHits = cms.InputTag("reducedEcalRecHitsEE"),
    EBRecHits = cms.InputTag("reducedEcalRecHitsEB"),
#    ReducedHBHERecHits = cms.InputTag("reducedHcalRecHits"),
    isMC = cms.untracked.bool(options.isMC),
    isSig = cms.untracked.bool(options.isSig),
    hasDpho = cms.untracked.bool(options.hasDpho),
    ecalBrem = cms.untracked.bool(options.ecalBrem),
    sameSign = cms.untracked.bool(options.sameSign),
    offPeakCR = cms.untracked.bool(options.offPeakCR),
    ttbarCR = cms.untracked.bool(options.ttbarCR),
    missingCSC = cms.untracked.bool(options.missingCSC),
    highHcalCr = cms.untracked.bool(options.highHcalCr),
    studyBarrelMuons = cms.untracked.bool(options.studyBarrel),
    tagAlignedCr = cms.untracked.bool(options.tagAlignedCr),
    scaleFactor = cms.untracked.double(options.scaleFactor),
    binSFsfname = cms.untracked.string(options.binSFsfname),
    metFailName = cms.untracked.string(options.metFailName),
    runrandomTrackEfficiency = cms.untracked.bool(options.runRandomTrack)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)


process.p = cms.Path(process.demo)
