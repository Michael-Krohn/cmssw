import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v17')

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

process.options = cms.untracked.PSet(
 SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
 IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
 #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    #fileNames = cms.untracked.vstring("root://cmseos.fnal.gov//store/user/revering/DYJetsToMM_EndcapDBrem_mAp1p5_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_revering/RunIISummer20UL18_SIM/220429_162528/0000/SIM_1.root")
    fileNames = readFiles
)

process.demo = cms.EDAnalyzer('DBremMatChecker',
    dbremweighttag = cms.InputTag("g4SimHits","DBremEventWeight"),
    dbremInitialETag = cms.InputTag("g4SimHits","DBremInitialEnergy"),
    dbremXTag = cms.InputTag("g4SimHits","DBremLocationX"),
    dbremYTag = cms.InputTag("g4SimHits","DBremLocationY"),
    dbremZTag = cms.InputTag("g4SimHits","DBremLocationZ"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)

process.p = cms.Path(process.demo)
