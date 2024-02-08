import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
options.register( 'runLocally',
                            False,
                            VarParsing.multiplicity.singleton,
                            VarParsing.varType.bool,
                            "True if running locally")

options.register( 'isMC',
                            True,
                            VarParsing.multiplicity.singleton,
                            VarParsing.varType.bool,
                            "True if using MC")

options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff') 
process.load("Configuration.StandardSequences.MagneticField_cff") 
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.isMC):
        process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
else:
        process.GlobalTag.globaltag = cms.string('106X_dataRun2_v32')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Uncomment when running over condor
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/local/cms/user/revering/ZmmSim/ZmmSim_miniAOD/ZmmSim_miniAOD_000.root'))

    # duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    # replace 'myfile.root' with the source file you want to use

process.myMetFilter = cms.EDFilter('METFilter_AOD',
    TriggerResultsTag = cms.InputTag("TriggerResults","","PAT")
)

from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
        muons = cms.InputTag("slimmedMuons"),
        vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
        PFCandidates = cms.InputTag("packedPFCandidates"),
        minDzBestTrack = cms.double(0.5),
        taggingMode    = cms.bool(True)
)

process.USER = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('p1')
    ),
    fileName = cms.untracked.string('mupx_filtering.root')
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("metFilterResults.root"),
				   closeFileFast = cms.untracked.bool(False)
)

process.p1 = cms.Path(process.BadPFMuonFilterUpdateDz*process.myMetFilter)
process.outpath = cms.EndPath(process.USER)
