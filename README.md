# MuonHLT Ntupler

## Run3 124X Recipe
```
export SCRAM_ARCH=slc7_amd64_gcc10
cmsrel CMSSW_12_4_11
cd CMSSW_12_4_11/src
cmsenv
git cms-init

git cms-addpkg HLTrigger/Configuration
git clone https://github.com/khaosmos93/MuonHLTForRun3.git HLTrigger/Configuration/python/MuonHLTForRun3

### Data
# Muon paths
hltGetConfiguration /dev/CMSSW_12_4_0/GRun \
 --process MYHLT \
 --data --globaltag auto:run3_hlt \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2022F/Muon/RAW-RECO/ZMu-PromptReco-v1/000/361/223/00000/7b6a8758-df56-43a7-88cb-e7392b433042.root \
 --eras Run3 \
 --max-events 100 \
 --full --offline --no-output >hlt_muon_data_Run2022.py

# Full menu
hltGetConfiguration /dev/CMSSW_12_4_0/GRun \
 --process MYHLT \
 --data --globaltag auto:run3_hlt \
 --unprescale \
 --input /store/data/Run2022F/Muon/RAW-RECO/ZMu-PromptReco-v1/000/361/223/00000/7b6a8758-df56-43a7-88cb-e7392b433042.root \
 --eras Run3 \
 --max-events 20 \
 --full --offline --no-output >hlt_muon_data_full_Run2022.py

cp hlt_muon_data_Run2022.py hlt_muon_data_Run2022_NoWP.py
sed -i 's/mvaCutB = cms.double( 0.01 )/mvaCutB = cms.double( 0.0 )/g' hlt_muon_data_Run2022_NoWP.py
sed -i 's/mvaCutE = cms.double( 0.01 )/mvaCutE = cms.double( 0.0 )/g' hlt_muon_data_Run2022_NoWP.py

### MC
# Muon paths
hltGetConfiguration /dev/CMSSW_12_4_0/GRun \
 --process MYHLT \
 --mc --globaltag auto:run3_mc_GRun \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/mc/Run3Winter22DR/DYToLL_M-50_TuneCP5_13p6TeV-pythia8/AODSIM/L1TPU0to99FEVT_122X_mcRun3_2021_realistic_v9-v2/2830001/5ea533ed-d952-4d74-8c48-456405fdc67b.root \
 --eras Run3 --l1-emulator uGT --l1 L1Menu_Collisions2022_v1_3_0-d1_xml \
 --max-events 100 \
 --full --offline --no-output >hlt_muon_mc_Run3.py

# --input /store/mc/Run3Summer21DRPremix/DYToLL_M-50_TuneCP5_14TeV-pythia8/AODSIM/120X_mcRun3_2021_realistic_v6-v2/80003/8e9e1bc0-64c4-427a-9175-26b7290f93d3.root \

# Full menu
hltGetConfiguration /dev/CMSSW_12_4_0/GRun \
 --process MYHLT \
 --mc --globaltag auto:run3_mc_GRun \
 --unprescale \
 --input /store/mc/Run3Winter22DR/DYToLL_M-50_TuneCP5_13p6TeV-pythia8/AODSIM/L1TPU0to99FEVT_122X_mcRun3_2021_realistic_v9-v2/2830001/5ea533ed-d952-4d74-8c48-456405fdc67b.root \
 --eras Run3 --l1-emulator uGT --l1 L1Menu_Collisions2022_v1_3_0-d1_xml \
 --max-events 20 \
 --full --offline --no-output >hlt_muon_mc_full_Run3.py

cp hlt_muon_mc_Run3.py hlt_muon_mc_Run3_NoWP.py
sed -i 's/mvaCutB = cms.double( 0.01 )/mvaCutB = cms.double( 0.0 )/g' hlt_muon_mc_Run3_NoWP.py
sed -i 's/mvaCutE = cms.double( 0.01 )/mvaCutE = cms.double( 0.0 )/g' hlt_muon_mc_Run3_NoWP.py

### 2018Data
# Muon paths
hltGetConfiguration /dev/CMSSW_12_4_0/GRun \
 --process MYHLT \
 --data --globaltag auto:run3_hlt \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-12Nov2019_UL2018-v6/40014/F7DED4A7-8B3A-574C-8805-2F4061C87ADA.root \
 --customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2018Input \
 --eras Run3 \
 --max-events 100 \
 --full --offline --no-output >hlt_muon_data_Run2018.py

# Full menu
hltGetConfiguration /dev/CMSSW_12_4_0/GRun \
 --process MYHLT \
 --data --globaltag auto:run3_hlt \
 --unprescale \
 --input /store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-12Nov2019_UL2018-v6/40014/F7DED4A7-8B3A-574C-8805-2F4061C87ADA.root \
 --customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2018Input \
 --eras Run3 \
 --max-events 20 \
 --full --offline --no-output >hlt_muon_data_full_Run2018.py

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
isMC   = False         # set True (False) for DY (data)
isDIGI = False         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW)

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", isMC, isDIGI, "")
process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")
from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", isMC, isDIGI, "")
process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

process.HLT_IsoMu24_v15 = cms.Path(
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
     process.HLT_IsoMu24_v15,
     process.HLT_Mu50_v15,
     process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v7,
     process.HLT_Mu15_v5,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myendpath,
     process.myseedpath
)
```
# after_menu_mc.sh
```
isMC   = True         # set True (False) for DY (data)
isDIGI = False         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW)

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", isMC, isDIGI, "")
process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")
from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", isMC, isDIGI, "")
process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

process.HLT_IsoMu24_v15 = cms.Path(
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
        '/store/mc/Run3Winter22DR/DYToLL_M-50_TuneCP5_13p6TeV-pythia8/AODSIM/L1TPU0to99FEVT_122X_mcRun3_2021_realistic_v9-v2/2830001/5ea533ed-d952-4d74-8c48-456405fdc67b.root',
    ),
    secondaryFileNames=cms.untracked.vstring(
        '/store/mc/Run3Winter22DR/DYToLL_M-50_TuneCP5_13p6TeV-pythia8/GEN-SIM-DIGI-RAW/L1TPU0to99FEVT_122X_mcRun3_2021_realistic_v9-v2/2830040/e5fa21b1-f07d-4f72-995a-67ac577814da.root',
    ),  
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v15,
     process.HLT_Mu50_v15,
     process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v7,
     process.HLT_Mu15_v5,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myendpath,
     process.myseedpath
)
```