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
        # ("DYToLL_M50_PU200_110X",     "/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("DYToLL_M10to50_PU200_110X", "/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("WToLNu",                    "/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),

        # ("TTToSemiLep",               "/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("TTTo2L2Nu",                 "/TTTo2L2Nu_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("TT",                        "/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),

        # ("QCD_Pt30to50",   "/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt50to80",   "/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt80to120",  "/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt120to170", "/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt170to300", "/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt300to470", "/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt470to600", "/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt600toInf", "/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),

        # ("QCD_Pt15to20_MuEn",   "/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v3/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt20to30_MuEn",   "/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v3/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt30to50_MuEn",   "/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt50to80_MuEn",   "/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt80to120_MuEn",  "/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt120to170_MuEn", "/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt170to300_MuEn", "/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW"),
        # ("QCD_Pt300toInf_MuEn", "/QCD_Pt-300toInf_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v3/GEN-SIM-DIGI-RAW")
    ]

    HLT_menus = [
        "HLT_MC_Mu.py",
    ]

    # proxy = '"/tmp/x509up_u95096"'

    for menu in HLT_menus:
        # menuTag = menu.replace(".py", "").replace("HLT_Phase2D49_", "")
        menuTag = "HLT2018"

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


