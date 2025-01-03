from CRABClient.UserUtilities import config
import sys

config = config()


#**************************submit function***********************                                                      
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
'''
from httplib import HTTPException
def submit(config):
        try:
                crabCommand('submit', config = config)
        except HTTPException as hte:
                print("Failed submitting task: %s" % (hte.headers))
        except ClientException as cle:
                print("Failed submitting task: %s" % (cle))
'''
#****************************************************************                                




# Common configuration

config.General.workArea     = 'crab_projects_ntuples'
config.General.transferLogs = False
config.JobType.pluginName   = 'Analysis' # PrivateMC
#    config.JobType.inputFiles   = ['Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016V4_MC.db']
config.JobType.sendExternalFolder = True
config.Data.inputDBS        = 'global'    
#config.Data.splitting       = 'LumiBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
config.Data.splitting       = 'FileBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
config.Data.totalUnits      = -1
config.Data.allowNonValidInputDataset      = True
config.Data.publication     = False
config.Site.storageSite     = 'T2_CH_CERN'
config.JobType.allowUndistributedCMSSW = True
#config.Data.useParent                   = True
config.Data.useParent                   = False
#config.Data.inputBlocks = '/DoubleElectron_FlatPt-1To300/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/GEN-SIM-RAW#d5be51ad-1f40-4565-a8cd-fe0de08338f3'
# dataset dependent configuration

config.JobType.psetName     = 'ConfFile_cfg.py'
config.General.requestName = 'mc'
config.Data.unitsPerJob    = 5
config.Data.useParent                   = True
#config.Data.inputDataset   = '/DoubleElectron_FlatPt-1To300/RunIIWinter19PFCalibDR-2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1/AODSIM'
#config.Data.inputDataset   = '/DoubleElectron_FlatPt-1To300/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/AODSIM'
#config.Data.inputDataset   = '/DoubleElectron_FlatPt-1To300/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/GEN-SIM-RAW'
config.Data.inputDataset   = '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter22DRPremix-EcalUncalZElectron-ALCARECO_122X_mcRun3_2021_realistic_v9_ext1-v3/ALCARECO'
config.Data.outLFNDirBase  = '/store/user/shilpi/ECAL_EM_noise/with5x5/'

### used 
#config.Data.inputDataset   = '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter22DRPremix-EcalUncalZElectron-ALCARECO_122X_mcRun3_2021_realistic_v9_ext1-v3/ALCARECO'

### 11th nov, 2023 - to initiate the transfer
#config.Data.inputDataset   = '/AlCaPhiSym/Run2022E-v1/RAW'

### to initiate the data transfer - 11th Dec, 2023
#config.Data.inputDataset   = '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter22DRPremix-ALCARECO_122X_mcRun3_2021_realistic_v9_ext1-v3/GEN-SIM-RAW'

'''
config.JobType.psetName     = 'ConfFile_cfg_minBias.py'
#config.Data.inputDataset   = '/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter22DR-TkAlMinBias-PUForTRKAPE_TRK_122X_mcRun3_2021_design_v9-v2/ALCARECO'
config.Data.inputDataset   = '/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter22DR-PUForTRKAPE_TRK_122X_mcRun3_2021_design_v9-v2/GEN-SIM-RAW'
 
#config.Data.outLFNDirBase  = '/store/user/shilpi' ### first version here
config.Data.outLFNDirBase  = '/store/user/shilpi/ECAL_EM_noise/with5x5/'
#submit(config)
'''
