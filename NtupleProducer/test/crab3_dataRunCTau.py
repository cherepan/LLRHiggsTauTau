# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Production_09_2017Tau'
config.General.workArea = 'crab3_dataRunCTau'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py' # to produce LLR ntuples or EnrichedMiniAOD according to the RunNtuplizer bool
config.JobType.sendExternalFolder = True #Needed until the PR including the Spring16 ele MVA ID is integrated in CMSSW/cms-data.

config.section_("Data")
config.Data.inputDataset = '/Tau/Run2016C-23Sep2016-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40 #number of events per jobs # 18K FOR SINGLE ELE, 10k for others
config.Data.totalUnits = -1 #number of event
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.runRange = '271036-284044' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/cherepan'
config.Data.publication = True
config.Data.outputDatasetTag = 'dataC_Production_09_2017'

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'
