# MuonHLT Ntupler

## Run3 121X Recipe
```
cmsrel CMSSW_12_1_0
cd CMSSW_12_1_0/src
cmsenv
git cms-init

git cms-addpkg DataFormats/SiPixelDetId
vi CommonTools/RecoAlgos/interface/TrackFullCloneSelectorBase.h
# constexpr static Packing thePacking = {11, 11, 0, 10};
git cms-checkdeps -a -A
scram b -j 8

git cms-addpkg RecoMuon/TrackerSeedGenerator
git clone -b dev https://github.com/wonpoint4/RecoMuon-TrackerSeedGenerator.git RecoMuon/TrackerSeedGenerator/data

git cms-addpkg HLTrigger/Configuration
git clone https://github.com/khaosmos93/MuonHLTForRun3.git HLTrigger/Configuration/python/MuonHLTForRun3

## Use cmsDriver for ntupler (due to sim hit matching...)

hltGetConfiguration /dev/CMSSW_12_1_0/GRun/V12 --type GRun \
--path HLTriggerFirstPath,HLT_IsoMu24_v*,HLT_Mu50_v*,HLTriggerFinalPath,HLTAnalyzerEndpath \
--unprescale --cff >$CMSSW_BASE/src/HLTrigger/Configuration/python/HLT_MuonHLT_cff.py

cmsDriver.py hlt_muon \
--python_filename=hlt_muon_Run3_mc.py \
--step HLT:MuonHLT \
--process MYHLT --era=Run3 \
--mc --conditions=auto:phase1_2021_realistic \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeDoubleMuIsoFix \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForDoubletRemoval \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForCscSegment \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForGEM \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForPatatrackWithIsoAndTriplets \
--customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeIOSeedingPatatrack \
--filein=/store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120003/e786c41e-21ba-489f-880c-42d0a248e59e.root \
-n 100 --no_output --no_exec
#--customise=RecoMuon/TrackerSeedGenerator/customizeOIseeding.customizeOIseeding \

## direct hltGetConfiguration

hltGetConfiguration /dev/CMSSW_12_1_0/GRun/V12 --type GRun \
 --process MYHLT \
 --mc --globaltag auto:phase1_2021_realistic \
 --unprescale \
 --paths HLTriggerFirstPath,HLT_IsoMu24_v*,HLT_Mu50_v*,HLTriggerFinalPath,HLTAnalyzerEndpath \
 --eras Run3 \
 --input /store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120003/e786c41e-21ba-489f-880c-42d0a248e59e.root \
 --full --offline --no-output >hlt_muon_run3_mc.py

from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3 import customizeDoubleMuIsoFix,customizeMuonHLTForDoubletRemoval,customizeMuonHLTForCscSegment,customizeMuonHLTForGEM,customizeMuonHLTForPatatrackWithIsoAndTriplets,customizeIOSeedingPatatrack
process = customizeDoubleMuIsoFix(process)
process = customizeMuonHLTForDoubletRemoval(process)
process = customizeMuonHLTForCscSegment(process)
process = customizeMuonHLTForGEM(process)
process = customizeMuonHLTForPatatrackWithIsoAndTriplets(process)
process = customizeIOSeedingPatatrack(process)  # update here for other WPs

#from RecoMuon.TrackerSeedGenerator.customizeOIseeding import customizeOIseeding
#process = customizeOIseeding(process)

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
# -- L2 seed stat recovery -- #
#process.hltIterL3MuonPixelTracksTrackingRegions.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.hltIter3IterL3MuonL2Candidates.src = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.hltL3MuonsIterL3IO.L3TrajBuilderParameters.MuonTrackingRegionBuilder.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence = cms.Sequence(
#    process.HLTIterL3OImuonTkCandidateSequence +
#    process.hltIterL3OIL3MuonsLinksCombination +
#    process.hltIterL3OIL3Muons +
#    process.hltIterL3OIL3MuonCandidates +
#    #process.hltL2SelectorForL3IO +
#    process.HLTIterL3IOmuonTkCandidateSequence +
#    process.hltIterL3MuonsFromL2LinksCombination
#)

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

skimDY = True         # set True (False) for DY (otherwise)
isDIGI = True         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW)
MvaVersion = "v8Pre"  # set v7Fast or v8Pre
WP = 0.0
doOI = False

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", skimDY, isDIGI, MvaVersion)
process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")
from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", skimDY, isDIGI, MvaVersion)
process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v13,
     process.HLT_Mu50_v13,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myseedpath
)
# Automatic addition of the customisation function from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3 import customizeDoubleMuIsoFix,customizeMuonHLTForDoubletRemoval,customizeMuonHLTForCscSegment,customizeMuonHLTForGEM,customizeMuonHLTForPatatrackWithIsoAndTriplets,customizeIOSeedingPatatrack
process = customizeDoubleMuIsoFix(process)
process = customizeMuonHLTForDoubletRemoval(process)
process = customizeMuonHLTForCscSegment(process)
process = customizeMuonHLTForGEM(process)
process = customizeMuonHLTForPatatrackWithIsoAndTriplets(process)
process = customizeIOSeedingPatatrack(process, mvaCutBs = (WP, WP), mvaCutEs = (WP, WP))

if doOI == True :
  from RecoMuon.TrackerSeedGenerator.customizeOIseeding import customizeOIseeding
  process = customizeOIseeding(process)

## Test run
cmsRun hlt_muon_Run3_mc.py
```
