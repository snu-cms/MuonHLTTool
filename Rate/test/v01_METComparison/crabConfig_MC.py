from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.workArea = 'CRABDir'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = ''
config.JobType.numCores = 4
# config.JobType.maxMemoryMB = 2500
# config.JobType.maxJobRuntimeMin = 2000

config.Data.inputDataset = ''
# config.Data.useParent = True

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Site.storageSite = 'T3_KR_KISTI'

# -- for CMSSW_10_1_X
# -- Ref: https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/3731/1/1/1/1/1/2/1/4.html
config.section_("General")
config.General.instance = 'preprod'

version = '_v20180424_'
# 'MultiCRAB' part
if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    config.General.requestName = 'METPlots'+version+'ttbar92XForTSG_OldMuonReco'
    config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/GEN-SIM-DIGI-RAW'
    config.JobType.psetName = 'HLTOpenCfgMET_OldMuonReco.py'
    crabCommand('submit', config = config)

    config.General.requestName = 'METPlots'+version+'ttbar92XForTSG_NewMuonReco'
    config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/GEN-SIM-DIGI-RAW'
    config.JobType.psetName = 'HLTOpenCfgMET_NewMuonReco.py'
    crabCommand('submit', config = config)

    config.General.requestName = 'METPlots'+version+'ttbar92XForTSG_IDonhltMuons'
    config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/GEN-SIM-DIGI-RAW'
    config.JobType.psetName = 'HLTOpenCfgMET_IDonhltMuons.py'
    crabCommand('submit', config = config)

    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)