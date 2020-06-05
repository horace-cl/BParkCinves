from CRABClient.UserUtilities import config

config=config()

#A name the user gives to it's request/task. 
#In particular, it is used by CRAB to create a project directory 
#(named crab_<requestName>) 
config.General.requestName = 'Bplus_parkedData_20-05-01'

#The area (full or relative path) where to create the CRAB project directory.
config.General.workArea = 'crab_projects'

#Whether or not to transfer the output files to the storage site.
config.General.transferOutputs = True

# Whether or not to copy the jobs log files to the storage site
config.General.transferLogs = True





# Specifies if this task is running an analysis ('Analysis') on an existing dataset
# or is running MC event generation ('PrivateMC').
config.JobType.pluginName = 'Analysis'

#The name of the CMSSW parameter-set configuration file that should be run via cmsRun.
config.JobType.psetName = 'run_BToKmumu_cfg.py'

#List of parameters to pass to the CMSSW parameter-set configuration file
#['myOption','param1=value1','param2=value2']
#config.JobType.pyCfgParams = ''

config.JobType.allowUndistributedCMSSW = True

#config.JobType.maxMemoryMB = 500



config.Data.inputDataset = '/ParkingBPH4/Run2018B-05May2019-v2/MINIAOD'
#config.Data.inputDataset = '/ParkingBPH5/Run2018D-05May2019promptD-v1/MINIAOD' 
config.Data.inputDBS = 'global'
#config.Data.splitting ='Automatic'# 'LumiBased'
config.Data.splitting ='LumiBased'
config.Data.unitsPerJob = 25
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '317392,317340,317320' # '193093-194075'
#config.Data.publication = True
#config.Data.outputDatasetTag = 'CRAB3_tutorial_May2015_Data_analysis'


#config.Data.inputDataset = '/ParkingBPH5/Run2018D-05May2019promptD-v1/MINIAOD' 
#config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outputDatasetTag = 'Bplus-Kplus-mumu-all'

#config.Site.storageSite = 'T3_MX_Cinvestav'
config.Site.storageSite = 'T2_CH_CERNBOX'



