from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Dphoton_MuPlusXSkim_Run2018B_lumisplit'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ConfFile_MuPlusX_AOD_cfg.py'
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/mireveri/MuPlusX/'
config.Data.inputDataset = '/SingleMuon/Run2018B-12Nov2019_UL2018-v2/AOD'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
#config.Data.userInputFiles = open('/uscms_data/d3/revering/MuPxFilter/CMSSW_10_6_17_patch1/src/DarkPhoton/MuFilter/Run2018C.txt').readlines()
config.Data.publication = False
config.Data.outputDatasetTag = 'Run2018B'
config.Data.allowNonValidInputDataset = True

config.Site.blacklist = ['T0_*']
config.Site.storageSite = "T3_US_FNALLPC"

