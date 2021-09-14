# MuonHLT Ntupler

## Run3 113X Recipe
```
cmsrel CMSSW_11_3_2
cd CMSSW_11_3_2/src
cmsenv
git cms-init

# Muon HLT customizers for Run 3
git cms-addpkg HLTrigger/Configuration
git clone -b dev https://github.com/khaosmos93/MuonHLTForRun3.git HLTrigger/Configuration/python/MuonHLTForRun3

# Simple Muon HLT menu
hltGetConfiguration /dev/CMSSW_11_3_0/GRun/V14 --type GRun \
--path HLTriggerFirstPath,HLT_IsoMu24_v*,HLT_Mu50_v*,HLTriggerFinalPath,HLTAnalyzerEndpath \
--unprescale --cff >$CMSSW_BASE/src/HLTrigger/Configuration/python/HLT_MuonHLT_cff.py

scram b -j 8

# cmsDriver
cmsDriver.py hlt_muon \
--python_filename=hlt_muon_Run3_mc.py \
--step HLT:MuonHLT \
--process MYHLT --era=Run3 \
--mc --conditions=113X_mcRun3_2021_realistic_v10 \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForDoubletRemoval \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForCscSegment \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForGEM \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizerFuncForMuonHLTSeeding \
--filein=root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI\
-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120003/e786c41e-21ba-489f-880c-42d0a248e59e.root \
-n 100 --no_output --no_exec

# in order to change parameters of the new muon seed classifier,
# e.g. to choose top 20 seeds with the highest quality, modify
# the following line in hlt_muon_Run3_mc.py
process = customizerFuncForMuonHLTSeeding(process)
 ->
process = customizerFuncForMuonHLTSeeding(process, newProcessName='MYHLT', doSort=True, nSeedsMaxBs = (20, 20), nSeedsMaxEs = (20, 20), mvaCutBs = (0.0, 0.0), mvaCutEs = (0.0, 0.0))

# optionally, add the following lines at the end of the configuration file to print out full trigger reports
process.options.wantSummary = cms.untracked.bool( True )
if 'MessageLogger' in process.__dict__:
     process.MessageLogger.TriggerSummaryProducerAOD = cms.untracked.PSet()
     process.MessageLogger.L1GtTrigReport = cms.untracked.PSet()
     process.MessageLogger.L1TGlobalSummary = cms.untracked.PSet()
     process.MessageLogger.HLTrigReport = cms.untracked.PSet()
     process.MessageLogger.FastReport = cms.untracked.PSet()
     process.MessageLogger.ThroughputService = cms.untracked.PSet()
     process.MessageLogger.cerr.FwkReport.reportEvery = 10

## Test run
cmsRun hlt_muon_Run3_mc.py
```

## Ntupler
```
git clone -b Run3 git@github.com:snu-cms/MuonHLTTool.git
scram b -j8
```
In hlt_muon_Run3_mc.py, Add these line
```
process.load("Configuration.StandardSequences.Reconstruction_cff")
# -- Ignoring filters -- #
process.HLT_Mu50_v13 = cms.Path(
    process.HLTBeginSequence +
    cms.ignore(process.hltL1sSingleMu22or25) +
    cms.ignore(process.hltPreMu50) +
    cms.ignore(process.hltL1fL1sMu22or25L1Filtered0) +
    process.HLTL2muonrecoSequence +
    cms.ignore(process.hltL2fL1sMu22or25L1f0L2Filtered10Q) +
    process.HLTL3muonrecoSequence +
    cms.ignore(process.hltL1fForIterL3L1fL1sMu22or25L1Filtered0) +
    cms.ignore(process.hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q) +
    process.HLTEndSequence
)

skimDY = True  # set True (False) for DY (otherwise)
isDIGI = True  # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW)
from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", skimDY, isDIGI, "Run3v6")
process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeedsFiltered",       "", "MYHLT")
process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeedsFiltered", "", "MYHLT")
from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", skimDY, isDIGI, "Run3v6")
process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeedsFiltered",       "", "MYHLT")
process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeedsFiltered", "", "MYHLT")

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v13,
     process.HLT_Mu50_v13,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myseedpath
)

## Test run
cmsRun hlt_muon_Run3_mc.py
```

## Timing
In hlt_muon_Run3_mc.py, Add these line (Do not run ntupler? -> use If-Else)
```
# -- Timing -- # (Copied from CMSSW_11_0_0 Menu)
# configure the FastTimerService
process.load( "HLTrigger.Timer.FastTimerService_cfi" )

# print a text summary at the end of the job
process.FastTimerService.printEventSummary         = False
process.FastTimerService.printRunSummary           = False
process.FastTimerService.printJobSummary           = True

# enable DQM plots
process.FastTimerService.enableDQM                 = True

# enable per-path DQM plots (starting with CMSSW 9.2.3-patch2)
process.FastTimerService.enableDQMbyPath           = True

# enable per-module DQM plots
process.FastTimerService.enableDQMbyModule         = True

# enable per-event DQM plots vs lumisection
process.FastTimerService.enableDQMbyLumiSection    = True
process.FastTimerService.dqmLumiSectionsRange      = 2500

# set the time resolution of the DQM plots
process.FastTimerService.dqmTimeRange              = 2000.
process.FastTimerService.dqmTimeResolution         =   10.
process.FastTimerService.dqmPathTimeRange          = 1000.
process.FastTimerService.dqmPathTimeResolution     =    5.
process.FastTimerService.dqmModuleTimeRange        =  200.
process.FastTimerService.dqmModuleTimeResolution   =    1.

# set the base DQM folder for the plots
process.FastTimerService.dqmPath                   = 'HLT/TimerService'
process.FastTimerService.enableDQMbyProcesses      = False

# load the DQMStore and DQMRootOutputModule
process.load( "DQMServices.Core.DQMStore_cfi" )
process.DQMStore.enableMultiThread = True

process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
    fileName = cms.untracked.string("DQMIO.root")
)

process.DQMOutput = cms.EndPath( process.dqmOutput )
```