#!/bin/bash

# -- data (the last pp menu for Run2)
hltGetConfiguration orcoff:/cdaq/physics/Run2018/2e34/v3.6.1/HLT/V2 \
--globaltag 101X_dataRun2_HLT_v7 \
--path HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
--input file:/afs/cern.ch/user/k/kplee/work/public/ROOTFile_Test/RAW_SingleMuon_Run2018A_Run316361.root \
--process MYHLT --full --offline \
--l1-emulator uGT \
--l1 L1Menu_Collisions2018_v2_1_0-d1_xml \
--customise MuonHLTTool/MuonHLTNtupler/customizerForMuonHLTNtupler.customizerFuncForMuonHLTNtupler \
--prescale none --max-events 100 --output none > HLTCfgData_Mu.py


# -- MC (92X)
hltGetConfiguration orcoff:/cdaq/physics/Run2018/2e34/v3.6.1/HLT/V2 \
--globaltag 101X_mc2017_realistic_TSG_2018_04_09_20_43_53 \
--path HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
--input file:/afs/cern.ch/user/k/kplee/work/public/ROOTFile_Test/ZMuMu_M50to120_92X_GENSIMRAW.root \
--mc --process MYHLT --full --offline \
--l1-emulator uGT \
--l1 L1Menu_Collisions2018_v2_1_0-d1_xml \
--customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2017DtUnpacking,\
MuonHLTTool/MuonHLTNtupler/customizerForMuonHLTNtupler.customizerFuncForMuonHLTNtupler \
--prescale none --max-events 100 --output none > HLTCfgMC_Mu.py
