import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  True,
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
        process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v10')
else:
        process.GlobalTag.globaltag = cms.string('94X_dataRun2_ReReco_EOY17_v6')

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#Uncomment when running over condor
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = readFiles
)

if not(options.isMC):
        process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt').getVLuminosityBlockRange()


process.DiMuonFilter = cms.EDFilter('MuFilter_AOD',
    recoMuons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    isMC = cms.untracked.bool(False)
)

process.USER = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('p1')
    ),
    fileName = cms.untracked.string('test_filtering.root')
)

#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("DiMuon_Histos.root"),
#				   closeFileFast = cms.untracked.bool(False)
#)

process.p1 = cms.Path(process.DiMuonFilter)
process.outpath = cms.EndPath(process.USER)
