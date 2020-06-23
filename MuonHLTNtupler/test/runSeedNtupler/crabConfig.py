from CRABClient.UserUtilities import config
import sys, os
import datetime
now = datetime.datetime.now()
date = now.strftime('%Y%m%d')

def doSkim(menu, flag):

    fin = open(menu, "r")
    fData = fin.read()
    fin.close()

    if flag:
        fData = fData.replace('customizerFuncForMuonHLTNtupler(process, "MYHLT", False)', 'customizerFuncForMuonHLTNtupler(process, "MYHLT", True)')
        fData = fData.replace('customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", False)', 'customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", True)')
    else:
        fData = fData.replace('customizerFuncForMuonHLTNtupler(process, "MYHLT", True)', 'customizerFuncForMuonHLTNtupler(process, "MYHLT", False)')        
        fData = fData.replace('customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", True)', 'customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", False)')        

    fout = open(menu, "w")
    fout.write(fData)
    fout.close()


workDir = 'PhaseII'
submitVersion = 'BDTImp'
mainOutputDir = '/store/user/hkwon/%s/v%s/%s' % (workDir,submitVersion,date)

config = config()

config.General.requestName = ''
# config.General.workArea = 'CRABDir'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = ''

config.Data.inputDataset = ''
# config.Data.useParent = True

config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.JobType.maxMemoryMB = 4000
config.Data.outLFNDirBase = mainOutputDir
config.Data.publication = False
config.Data.ignoreLocality = True
config.Site.storageSite = 'T2_KR_KNU'
config.Site.whitelist = ['T2_CH_CERN','T2_FR_*']

# 'MultiCRAB' part
if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    HLT_menu = 'HLT_Phase2C9_trigger_D49.py'

    ## DY 110X ##

    doSkim(HLT_menu, True)

    config.General.workArea = 'crab_%s_%s_%s' % (workDir, submitVersion, date)
    config.General.requestName = '%s_%s_%s' % (workDir, 'DYMuMu_PU200_110X_pilot', date)
    config.JobType.psetName = HLT_menu
    config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    crabCommand('submit', config = config)
    print "\n\n", config

    config.General.workArea = 'crab_%s_%s_%s' % (workDir, submitVersion, date)
    config.General.requestName = '%s_%s_%s' % (workDir, 'DYMuMu_PU200_110X_pilot2', date)
    config.JobType.psetName = HLT_menu
    config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-PU200_pilot2_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    crabCommand('submit', config = config)
    print "\n\n", config

    config.General.workArea = 'crab_%s_%s_%s' % (workDir, submitVersion, date)
    config.General.requestName = '%s_%s_%s' % (workDir, 'DYMuMu_FlatPU0To200_110X_pilot', date)
    config.JobType.psetName = HLT_menu
    config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-FlatPU0To200_pilot_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    crabCommand('submit', config = config)
    print "\n\n", config

    ## ttbar ##

    doSkim(HLT_menu, False)

    config.General.workArea = 'crab_%s_%s_%s' % (workDir, submitVersion, date)
    config.General.requestName = '%s_%s_%s' % (workDir, 'TTToSemiLepton_PU200_110X', date)
    config.JobType.psetName = 'HLT_Phase2C9_trigger_D49.py'
    config.Data.inputDataset = '/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    crabCommand('submit', config = config)
    print "\n\n", config

    config.General.workArea = 'crab_%s_%s_%s' % (workDir, submitVersion, date)
    config.General.requestName = '%s_%s_%s' % (workDir, 'TTTo2L2Nu_PU200_110X', date)
    config.JobType.psetName = 'HLT_Phase2C9_trigger_D49.py'
    config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    crabCommand('submit', config = config)
    print "\n\n", config

    ## 106X ##

    # config.General.workArea = 'crab_%s_%s_%s' % (submitVersion, date, NtupleType)
    # config.General.requestName = '%s_%s_%s' % (submitVersion, 'DYMuMu_PU200_106X', date)
    # config.JobType.psetName = 'HLT.py'
    # config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/PhaseIITDRSpring19DR-PU200_pilot_106X_upgrade2023_realistic_v3-v1/GEN-SIM-DIGI-RAW'
    # crabCommand('submit', config = config)
    # print "\n\n", config

    # config.General.workArea = 'crab_%s_%s_%s' % (submitVersion, date, NtupleType)
    # config.General.requestName = '%s_%s_%s' % (submitVersion, 'DYMuMu_FlatPU0To200_106X', date)
    # config.JobType.psetName = 'HLT.py'
    # config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/PhaseIITDRSpring19DR-FlatPU0To200_pilot_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    # crabCommand('submit', config = config)
    # print "\n\n", config

    # config.General.workArea = 'crab_%s_%s_%s' % (submitVersion, date, NtupleType)
    # config.General.requestName = '%s_%s_%s' % (submitVersion, 'Mu_FlatPt2to100_PU200', date)
    # config.JobType.psetName = 'HLT.py'
    # config.Data.inputDataset = '/Mu_FlatPt2to100-pythia8-gun/PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    # crabCommand('submit', config = config)
    # print "\n\n", config

    # config.General.workArea = 'crab_%s_%s_%s' % (submitVersion, date, NtupleType)
    # config.General.requestName = '%s_%s_%s' % (submitVersion, 'Mu_FlatPt2to100_noPU', date)
    # config.JobType.psetName = 'HLT.py'
    # config.Data.inputDataset = '/Mu_FlatPt2to100-pythia8-gun/PhaseIITDRSpring19DR-NoPU_106X_upgrade2023_realistic_v3-v1/GEN-SIM-DIGI-RAW'
    # crabCommand('submit', config = config)
    # print "\n\n", config

    # config.General.workArea = 'crab_%s_%s_%s' % (submitVersion, date, NtupleType)
    # config.General.requestName = '%s_%s_%s' % (submitVersion, 'WToLNu_PU200', date)
    # config.JobType.psetName = 'HLT.py'
    # config.Data.inputDataset = '/WToLNu_14TeV_TuneCP5_pythia8/PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v1/GEN-SIM-DIGI-RAW'
    # crabCommand('submit', config = config)
    # print "\n\n", config

    # config.General.workArea = 'crab_%s_%s_%s' % (submitVersion, date, NtupleType)
    # config.General.requestName = '%s_%s_%s' % (submitVersion, 'WToLNu_PU140', date)
    # config.JobType.psetName = 'HLT.py'
    # config.Data.inputDataset = '/WToLNu_14TeV_TuneCP5_pythia8/PhaseIITDRSpring19DR-PU140_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    # crabCommand('submit', config = config)
    # print "\n\n", config

    # config.General.workArea = 'crab_%s_%s_%s' % (submitVersion, date, NtupleType)
    # config.General.requestName = '%s_%s_%s' % (submitVersion, 'WToLNu_NoPU', date)
    # config.JobType.psetName = 'HLT.py'
    # config.Data.inputDataset = '/WToLNu_14TeV_TuneCP5_pythia8/PhaseIITDRSpring19DR-NoPU_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW'
    # crabCommand('submit', config = config)
    # print "\n\n", config
