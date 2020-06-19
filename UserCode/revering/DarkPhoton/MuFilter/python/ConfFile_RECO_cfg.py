import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if is MC dataset")

options.register( 'is2017',
                                  True,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if using 2017 dataset")
											
options.register( 'fileName',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "File name to fiter")

options.parseArguments()
																	

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.isMC):
        process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v15')
else:
    if(options.is2017):
        process.GlobalTag.globaltag = cms.string('94X_dataRun2_ReReco_EOY17_v6')
    else:
        process.GlobalTag.globaltag = cms.string('101X_dataRun2_Prompt_v11')


process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#import FWCore.Utilities.FileUtils as FileUtils
#mylist = FileUtils.loadListFromFile ('files_SingleMuon2017EE-PromptReco-v1_RECO_temp.txt')
#readFiles = cms.untracked.vstring( *mylist)

#Uncomment when running over condor
process.source = cms.Source("PoolSource",
    inputCommands = cms.untracked.vstring('keep *',
    					  'drop recoPFJets_*_*_*',
					  'drop recoPFMETs_*_*_*',
					  'drop recoJPTJets_*_*_*',
					  'drop recoCaloMETs_*_*_*',
					  'drop recoPreshowerClusterShapes_*_*_*',
					  'drop recoPreshowerClusters_*_*_*',
					  'drop doubles_ak*_*_*',
					  'drop recoBasicJets_*_*_*',
					  'drop doubles_cmsTopTagPFJetsCHS_*_*',
					  'drop DcsStatuss_*_*_*',
					  'drop L1Gct*_*_*_*',
					  'drop long_tcdsDigis_*_*',
					  'drop int_tcdsDigis_*_*',
					  'drop recoJetIDedmValueMap_*_*_*',
					  'drop double_ak*_*_*',
					  'drop *_fixedGridRho*_*_*',
					  'drop *_reducedEcalRecHits*_*_*',
					  'drop booledmValueMap_Photon*_*_*',
					  'drop floatedmValueMap_eid*_*_*',
					  'drop floatedmValueMap_*Photon*_*_*',
					  'drop Totem*_*_*_*',
					  'drop recoJetedm*_*_*_*',
					  'drop *HaloData_*_*_*',
					  'drop double_cmsTopTag*_*_*',
					  'drop CTPPS*_*_*_*'),
					  #'drop RPCDet*_*_*_*'),
					  #'drop EBDigiCollection_*_*_*',
					  #'drop EEDigiCollection_*_*_*'),
					  #'drop recoCaloClusters_*_*_*'),
					  #'drop DetIdedm*_*_*_*'),
#					  'drop *_trackerDrivenElectronSeeds_*_*'),
					  #'drop *_ecal*_*_*'),
					  #'drop recoPhotonCores_*_*_*'),
					  #'drop recoGsfElectrons_*_*_*'),
					  #'drop recoSuperClusters_*_*_*'),
    					  #'drop l1tEGammaBXVector_*_*_*',
					  #'drop l1tEtSumBXVector_*_*_*',
					  #'drop l1tJetBXVector_*_*_*',
					  #'drop li1MuonBXVector_*_*_*',
					  #'drop double_*_*_*',
					  #'drop recoJetedmRefToBaseProdTofloatsAssociationVector_*_*_*',
					  #'drop EcalRecHitsSorted_*_*_*'),
					  #'drop booledmValueMap_*_*_*'),
#    					  'drop *',
#					  'keep *_*_*_HLT',
#					  'keep *_generalTracks_*_*',
#					  'keep recoMuons_muons_*_*',
#					  'keep *_offlinePrimaryVertices_*_*',
#					  'keep *_cscSegments_*_*',
#					  'keep *_hbhereco_*_*'),

    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring("file:"+options.fileName)
    #fileNames = cms.untracked.vstring("file:/local/cms/user/revering/dphoton/CMSSW_10_2_11_patch1/src/DarkPhoton/MuGenerator/test.root")
)

if not(options.isMC):
    if(options.is2017):
        process.source.lumisToProcess = LumiList.LumiList(filename = '/home/krohn045/MuMu/CMSSW_10_2_X/CMSSW_10_2_6/src/DarkPhoton/MuFilter/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt').getVLuminosityBlockRange()
    else:
        process.source.lumisToProcess = LumiList.LumiList(filename = '/home/krohn045/MuMu/CMSSW_10_2_X/CMSSW_10_2_6/src/DarkPhoton/MuFilter/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON_MuonPhys.txt').getVLuminosityBlockRange()

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
