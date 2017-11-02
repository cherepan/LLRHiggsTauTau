# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'SimFitTest_v2'
config.General.workArea = 'crab3_dy_ll_SimFitTest'
config.General.transferLogs = True


config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py' # to produce LLR ntuples or EnrichedMiniAOD according to the RunNtuplizer bool
config.JobType.sendExternalFolder = True #Needed until the PR including the Spring16 ele MVA ID is integrated in CMSSW/cms-data.

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 100000 #number of events per jobs # 18K FOR SINGLE ELE, 10k for others
config.Data.totalUnits = 5000000 #number of event
config.Data.outLFNDirBase = '/store/user/cherepan'
config.Data.publication = True
config.Data.outputDatasetTag = 'SimFitTest_v2'

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'
