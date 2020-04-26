# MuonHLT Ntupler

## L1 Phase2
	cmsrel CMSSW_11_1_0_pre6
	cd CMSSW_11_1_0_pre6/src
	cmsenv
	git cms-init
	git cms-merge-topic -u cms-L1TK:L1TK-integration-CMSSW_11_1_0_pre4
	git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v3.0.0-CMSSW_11_1_0_pre6
	git cms-addpkg L1Trigger/L1TCommon
## Phase2 MuonHLT
	git clone https://username@gitlab.cern.ch/cms-hlt/phase2.git
	mkdir -p HLTrigger/PhaseII/test
	cp -r phase2/cmssw-configs/HLT-v0/pythonFragments/ HLTrigger/PhaseII/python
## Ntupler
	git clone -b SeedingStudy git@github.com:snu-cms/MuonHLTTool.git
## Run
	cmsRun HLTrigger/PhaseII/test/HLT.py



