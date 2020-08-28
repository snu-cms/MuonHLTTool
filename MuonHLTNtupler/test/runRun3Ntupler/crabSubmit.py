from CRABClient.UserUtilities import config, getUsernameFromCRIC
import sys, os
import gc
import datetime
now = datetime.datetime.now()
date = now.strftime('%Y%m%d')

submitVersion = 'MuonHLTRun3'
mainOutputDir = '/store/user/%s/%s/%s' % (getUsernameFromCRIC(), submitVersion, date)



# 'MultiCRAB' part
if __name__ == '__main__':
    # from CRABAPI.RawCommand import crabCommand

    crab_cfg = """
from CRABClient.UserUtilities import config, getUsernameFromCRIC

config = config()

config.General.requestName = '%(datasetTag)s_%(menuTag)s_%(date)s'
config.General.workArea = 'crab_%(submitVersion)s_%(date)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '%(menu)s'

config.Data.inputDataset = '%(datasetPath)s'
# config.Data.useParent = True
config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.JobType.maxMemoryMB = 4000
config.Data.outLFNDirBase = '%(mainOutputDir)s'
config.Data.publication = False
config.Data.ignoreLocality = True
config.Site.storageSite = 'T2_KR_KNU'
config.Site.whitelist = ['T2_CH_CERN','T2_FR_*']
    """

    datasets = [
        ("DYToLL_M50_110X", "/DYToLL_M-50_TuneCP5_14TeV-pythia8_HCAL/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW")
    ]

    HLT_menus = [
        "HLT_MC_Run3.py",
    ]

    # proxy = '"/tmp/x509up_u95096"'

    for menu in HLT_menus:
        menuTag = 'HLTRun3'  # menu.replace(".py", "").replace("HLT_MC_", "")

        for datasetTag, datasetPath in datasets:

            Crab_Config = 'crabConfig_'+datasetTag+'.py'
            print "\n\n", crab_cfg % locals()
            sys.stdout.flush()
            gc.collect()

            open(Crab_Config, 'wt').write(crab_cfg % locals())

            cmd = 'crab submit -c '+Crab_Config  # +' --proxy='+proxy

            print cmd
            sys.stdout.flush()
            gc.collect()

            os.system(cmd)


