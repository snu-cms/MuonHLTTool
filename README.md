# MuonHLT Ntupler

## Run3 130X Recipe
```
cmsrel CMSSW_13_0_6
cd CMSSW_13_0_6/src
cmsenv
git cms-init

git cms-addpkg HLTrigger/Configuration
git clone -b Run2023 https://github.com/wonpoint4/MuonHLTForRun3.git HLTrigger/Configuration/python/MuonHLTForRun3

### Data (Efficiency)
hltGetConfiguration /dev/CMSSW_13_0_0/GRun \
 --process MYHLT \
 --data --globaltag 130X_dataRun3_HLT_v2 \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2022F/Muon/RAW-RECO/ZMu-PromptReco-v1/000/361/223/00000/7b6a8758-df56-43a7-88cb-e7392b433042.root \
 --eras Run3 \
 --max-events 100 \
 --full --offline --no-output >hlt_muon_data_Run2022.py


### MC (Efficiency)
hltGetConfiguration /dev/CMSSW_13_0_0/GRun \
 --process MYHLT \
 --mc --globaltag 126X_mcRun3_2023_forPU65_v3 \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/mc/Run3Winter23Reco/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/AODSIM/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40000/4b79f858-07ac-4aad-962f-e9473f3141a6.root \
 --eras Run3 \
 --max-events 100 \
 --full --offline --no-output >hlt_muon_mc_Run3.py


### 2018Data (Efficiency)
hltGetConfiguration /dev/CMSSW_13_0_0/GRun \
 --process MYHLT \
 --data --globaltag auto:run3_hlt \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-12Nov2019_UL2018-v6/40014/F7DED4A7-8B3A-574C-8805-2F4061C87ADA.root \
 --customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2018Input \
 --eras Run3 \
 --max-events 100 \
 --full --offline --no-output >hlt_muon_data_Run2018.py


## Test run
cmsRun hlt_muon_data_Run2022.py
cmsRun hlt_muon_mc_Run3.py
cmsRun hlt_muon_data_Run2018.py
```

## Ntupler
```
git clone -b Run3 git@github.com:snu-cms/MuonHLTTool.git
scram b -j8

vi after_menu_data.sh # See below
vi after_menu_mc.sh   #	See below

## Test run
cat after_menu_data.sh >> hlt_muon_data_Run2022.py
sed -i 's/numberOfThreads = 4/numberOfThreads = 1/g' hlt_muon_data_Run2022.py
cmsRun hlt_muon_data_Run2022.py

cat after_menu_mc.sh >> hlt_muon_mc_Run3.py
sed -i 's/numberOfThreads = 4/numberOfThreads = 1/g' hlt_muon_mc_Run3.py
cmsRun hlt_muon_mc_Run3.py
```

# after_menu_data.sh
```
isDIGI = False         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW or Data)

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", isDIGI)
process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

#from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
#process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", isDIGI)
#process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
#process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

# -- L2 seed stat recovery -- #
#process.hltIterL3MuonPixelTracksTrackingRegions.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
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

process.HLT_IsoMu24_v17 = cms.Path(
    process.HLTBeginSequence +
    cms.ignore(process.hltL1sSingleMu22) +
    cms.ignore(process.hltPreIsoMu24) +
    cms.ignore(process.hltL1fL1sMu22L1Filtered0) +
    process.HLTL2muonrecoSequence +
    cms.ignore(process.hltL2fL1sSingleMu22L1f0L2Filtered10Q) +
    process.HLTL3muonrecoSequence +
    cms.ignore(process.hltL1fForIterL3L1fL1sMu22L1Filtered0) +
    cms.ignore(process.hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q) +
    process.HLTMu24IsolationSequence +
    cms.ignore(process.hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered) +
    process.HLTEndSequence
)
process.HLTMu24IsolationSequence = cms.Sequence(
    process.HLTL3muonEcalPFisorecoSequenceNoBoolsForMuons +
    cms.ignore(process.hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfecalIsoRhoFiltered) +
    process.HLTL3muonHcalPFisorecoSequenceNoBoolsForMuons +
    cms.ignore(process.hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfhcalIsoRhoFiltered) +
    process.HLTTrackReconstructionForIsoL3MuonIter02 +
    process.hltMuonTkRelIsolationCut0p08Map
)

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v17,
     process.HLT_Mu50_v17,
     process.HLT_CascadeMu100_v7,
     process.HLT_HighPtTkMu100_v6,
     process.HLT_Mu15_v7,
     process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v9,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myendpath,
     #process.myseedpath
)
```
# after_menu_mc.sh
```
isDIGI = True         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW or Data)

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", isDIGI)
process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", isDIGI)
process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

# -- L2 seed stat recovery -- #
#process.hltIterL3MuonPixelTracksTrackingRegions.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
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

process.HLT_IsoMu24_v17 = cms.Path(
    process.HLTBeginSequence +
    cms.ignore(process.hltL1sSingleMu22) +
    cms.ignore(process.hltPreIsoMu24) +
    cms.ignore(process.hltL1fL1sMu22L1Filtered0) +
    process.HLTL2muonrecoSequence +
    cms.ignore(process.hltL2fL1sSingleMu22L1f0L2Filtered10Q) +
    process.HLTL3muonrecoSequence +
    cms.ignore(process.hltL1fForIterL3L1fL1sMu22L1Filtered0) +
    cms.ignore(process.hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q) +
    process.HLTMu24IsolationSequence +
    cms.ignore(process.hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered) +
    process.HLTEndSequence
)
process.HLTMu24IsolationSequence = cms.Sequence(
    process.HLTL3muonEcalPFisorecoSequenceNoBoolsForMuons +
    cms.ignore(process.hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfecalIsoRhoFiltered) +
    process.HLTL3muonHcalPFisorecoSequenceNoBoolsForMuons +
    cms.ignore(process.hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfhcalIsoRhoFiltered) +
    process.HLTTrackReconstructionForIsoL3MuonIter02 +
    process.hltMuonTkRelIsolationCut0p08Map
)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Run3Winter23Reco/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/AODSIM/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40000/4b79f858-07ac-4aad-962f-e9473f3141a6.root',
    ),
    secondaryFileNames=cms.untracked.vstring(
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/0f51e203-dfd7-4c0e-8f3d-552f1edefa80.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/142265be-53cf-4204-915d-487cd3ae3878.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/172430e4-84f4-4610-812b-b05ef45b6ca2.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/19287029-b952-4d44-9cea-30f0f111ee6c.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/89c74e90-e99d-4bcd-8718-c38b8e39828e.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/a4908eb5-2532-4eef-8832-2149628e8586.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/bb3fb0a6-a20b-43bf-b4b1-8bafca94aebd.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/db7e130d-3b8f-4bee-974d-a399c9d440d1.root',
        '/store/mc/Run3Winter23Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/40006/f682d635-1eec-4d64-96d3-88c4ee85b00a.root',
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v17,
     process.HLT_Mu50_v17,
     process.HLT_CascadeMu100_v7,
     process.HLT_HighPtTkMu100_v6,
     process.HLT_Mu15_v7,
     process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v9,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myendpath,
     process.myseedpath
)
```