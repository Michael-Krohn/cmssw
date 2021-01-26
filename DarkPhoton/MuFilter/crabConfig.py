from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Dphoton_MuPlusXSkim_WJets'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ConfFile_MuPlusX_AOD_cfg.py'
config.JobType.maxMemoryMB = 5000
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/mireveri/MuPlusX/'
#config.Data.inputDataset  '/SingleMuon/Run2018D-22Jan2019-v2/AOD'
config.Data.userInputFiles = open('/uscms_data/d3/revering/CMSSW_10_6_17_patch1/src/DarkPhoton/MuFilter/WJetsToLNu.txt').readlines()
config.Data.publication = False
config.Data.outputDatasetTag = 'WJets'
config.Data.allowNonValidInputDataset = True

config.Site.blacklist = ['T0_*']
config.Site.storageSite = "T3_US_FNALLPC"

