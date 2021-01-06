import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if is MC dataset")
options.parseArguments()
																	

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.isMC):
        process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v17')
else:
        process.GlobalTag.globaltag = cms.string('102X_dataRun2_Prompt_v11')

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('Run2018D_UL2018_v4.txt')
readFiles = cms.untracked.vstring( *mylist )

#Uncomment when running over condor
process.source = cms.Source("PoolSource",
    # duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    # replace 'myfile.root' with the source file you want to use
    fileNames = readFiles
    
)

#if not(options.isMC):
#        process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt').getVLuminosityBlockRange()


process.DiMuonFilter = cms.EDFilter('MuPlusXFilter_AOD',
    recoMuons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    HBHERecHits = cms.InputTag("reducedHcalRecHits","hbhereco"),
    EERecHits = cms.InputTag("reducedEcalRecHitsEE"),
    EBRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    genParticles = cms.InputTag("genParticles"),
    isMC = cms.untracked.bool(False)
)

process.USER = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('p1')
    ),
    fileName = cms.untracked.string('mupx_filtering.root')
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("skim_monitor_hists.root"),
				   closeFileFast = cms.untracked.bool(False)
)

process.p1 = cms.Path(process.DiMuonFilter)
process.outpath = cms.EndPath(process.USER)
