# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: RECO -s RAW2DIGI,L1Reco,RECO,ALCA:TkAlZMuMu+MuAlCalIsolatedMu+MuAlOverlaps+MuAlZMuMu+DtCalib+HcalCalLowPUHBHEMuonFilter,EI,PAT,DQM:@rerecoCommon --runUnscheduled --nThreads 8 --data --era Run2_2018,pf_badHcalMitigation --scenario pp --conditions 106X_dataRun2_v26 --eventcontent AOD,MINIAOD,DQM --datatier AOD,MINIAOD,DQMIO --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2018,Configuration/DataProcessing/Utils.addMonitoring --filein file:test.root -n 100 --python_filename=recoAOD_Run2018D_SingleMuon.py --no_exec
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if is MC dataset")

options.register( 'signalRegion',
                                  True,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if running on MC in isolated region")

options.parseArguments()


from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Modifier_pf_badHcalMitigation_cff import pf_badHcalMitigation

process = cms.Process('RECO',Run2_2018,pf_badHcalMitigation)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PAT_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile ('SingleMuon_Run2018D_RAW_1000files.txt')
readFiles = cms.untracked.vstring( *mylist)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:82651E35-F896-E811-8A0E-02163E01A133.root'),
#    fileNames = readFiles,
#    skipEvents = cms.untracked.uint32(773)
)

#process.source.lumisToProcess = LumiList.LumiList(filename = '/uscms_data/d3/mkrohn/DarkPhoton/RAWskim/CMSSW_10_6_17_patch1/src/DarkPhoton/MuFilter/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt').getVLuminosityBlockRange()


process.options = cms.untracked.PSet(

)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECO nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.DiMuonFilter = cms.EDFilter('MuPlusXFilter_HardBrem_RAW_RECO',
    recoMuons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    CSCSegmentLabel = cms.InputTag("cscSegments"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    HBHERecHits = cms.InputTag("hbhereco"),
#    HBHERecHits = cms.InputTag("reducedHcalRecHits","hbhereco"),
    EERecHits = cms.InputTag("reducedEcalRecHitsEE"),
    EBRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    CaloJetSource = cms.InputTag("ak4CaloJets"),
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    genParticles = cms.InputTag("genParticles"),
    isMC = cms.untracked.bool(options.isMC),
    signalRegion = cms.untracked.bool(options.signalRegion)
)

process.USER = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('p1')
    ),
    fileName = cms.untracked.string('mupx_RAW_filtering.root')
)

process.p1 = cms.Path(process.DiMuonFilter)
process.outpath = cms.EndPath(process.USER)

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('p1')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO_ALCA_EI_PAT_DQM.root'),
    outputCommands = process.AODEventContent.outputCommands +
        cms.untracked.vstring(
        'keep *_hbhereco_*_*'
    )
)


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v26', '')


# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.DiMuonFilterPath = cms.Path(process.DiMuonFilter)
process.dqmoffline_step = cms.EndPath(process.DQMOfflineCommon)
process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineMuon)
process.dqmoffline_2_step = cms.EndPath(process.DQMOfflineHcal)
process.dqmoffline_3_step = cms.EndPath(process.DQMOfflineJetMET)
process.dqmoffline_4_step = cms.EndPath(process.DQMOfflineEcal)
process.dqmoffline_5_step = cms.EndPath(process.DQMOfflineEGamma)
process.dqmoffline_6_step = cms.EndPath(process.DQMOfflineL1TMuon)
process.dqmoffline_7_step = cms.EndPath(process.DQMOfflineL1TEgamma)
process.dqmoffline_8_step = cms.EndPath(process.DQMOfflineCTPPS)
process.dqmoffline_9_step = cms.EndPath(process.DQMOfflineL1TMonitoring)
process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
process.AODoutput_step = cms.EndPath(process.AODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.pathALCARECOMuAlZMuMu,process.pathALCARECOMuAlZMuMuGeneralTracks,process.pathALCARECOMuAlOverlaps,process.pathALCARECOMuAlOverlapsGeneralTracks,process.pathALCARECOTkAlZMuMu,process.pathALCARECODtCalib,process.pathALCARECOHcalCalLowPUHBHEMuonFilter,process.pathALCARECOMuAlCalIsolatedMu,process.pathALCARECOMuAlCalIsolatedMuGeneralTracks,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.dqmoffline_step,process.dqmoffline_1_step,process.dqmoffline_2_step,process.dqmoffline_3_step,process.dqmoffline_4_step,process.dqmoffline_5_step,process.dqmoffline_6_step,process.dqmoffline_7_step,process.dqmoffline_8_step,process.dqmoffline_9_step,process.dqmofflineOnPAT_step,process.p1,process.AODoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2018

#call to customisation function customisePostEra_Run2_2018 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2018(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllData(process)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


