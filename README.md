# MuonHLT Ntupler

## L1 Track
	cmsrel CMSSW_11_1_0_pre4
	cd CMSSW_11_1_0_pre4/src/
	cmsenv
	git cms-checkout-topic cms-L1TK:L1TK-integration-CMSSW_11_1_0_pre4
## Phase2 MuonHLT
	git clone https://username@gitlab.cern.ch/cms-hlt/phase2.git
	mkdir -p HLTrigger/PhaseII/test
	cp -r phase2/cmssw-configs/HLT-v0/pythonFragments/ HLTrigger/PhaseII/python
## Ntupler
	git clone -b SeedingStudy git@github.com:snu-cms/MuonHLTTool.git
## Run
	cmsRun HLTrigger/PhaseII/test/HLT.py



