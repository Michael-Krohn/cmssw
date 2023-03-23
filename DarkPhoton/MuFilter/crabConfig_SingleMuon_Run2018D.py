from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Dphoton_MuPlusXSkim_SingleMuon_Run2018D_RAW_RECO_4thTry'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ConfFile_MuPlusX_RAW_RECO_cfg.py'
config.JobType.maxMemoryMB = 14000
config.JobType.maxJobRuntimeMin = 2000
config.JobType.allowUndistributedCMSSW = True
config.JobType.numCores = 8
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 2
#config.Data.totalUnits = 1000
config.Data.outLFNDirBase = '/store/user/mkrohn/MuPlusX/'
config.Data.inputDataset = '/SingleMuon/Run2018D-v1/RAW'
#config.Data.userInputFiles = open('SingleMuon_Run2018D_RAW_1000files.txt').readlines()
config.Data.publication = False
config.Data.outputDatasetTag = 'SingleMuon_Run2018D_RAW'
config.Data.allowNonValidInputDataset = True
config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

config.Site.blacklist = ['T0_*']
config.Site.storageSite = "T3_US_FNALLPC"

