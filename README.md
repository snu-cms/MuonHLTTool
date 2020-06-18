# MuonHLT Ntupler

## L1 Phase2
	https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_11_1_0_pre6

	cmsrel CMSSW_11_1_0_pre6
	cd CMSSW_11_1_0_pre6/src
	cmsenv
	git cms-init
	git cms-merge-topic -u cms-L1TK:L1TK-integration-CMSSW_11_1_0_pre4
	git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v3.0.2
	git cms-addpkg L1Trigger/L1TCommon

## Seed classifier
	git cms-addpkg HLTrigger/Muon
	git clone https://github.com/khaosmos93/MuonHLTSeedMVAClassifier.git HLTrigger/MuonHLTSeedMVAClassifier

## Ntupler
	git clone -b SeedingStudy git@github.com:snu-cms/MuonHLTTool.git
	scram b -j8

## Test run
	cd MuonHLTNtupler/test/runSeedNtupler/
	cmsRun HLT_Phase2C9_trigger_D49.py

## Phase2 MuonHLT - not used
	git clone https://username@gitlab.cern.ch/cms-hlt/phase2.git
	mkdir -p HLTrigger/PhaseII/test
	cp -r phase2/cmssw-configs/HLT-v0/pythonFragments/ HLTrigger/PhaseII/python


