import FWCore.ParameterSet.Config as cms

ntuplerBase = cms.EDAnalyzer("MuonHLTNtupler",
	offlineMuon = cms.untracked.InputTag("muons"),
	offlineVertex = cms.untracked.InputTag("offlinePrimaryVertices"),
	triggerResults = cms.untracked.InputTag("TriggerResults::HLT"),
	triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
	myTriggerResults = cms.untracked.InputTag("TriggerResults::MYHLT"), # -- rerun object -- #
	myTriggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"), # -- rerun object -- #
	L3Muon = cms.untracked.InputTag("hltIterL3MuonCandidates"),
	L2Muon = cms.untracked.InputTag("hltL2MuonCandidates"),
	# L1Muon = cms.untracked.InputTag("gmtStage2Digis", "Muon", "RECO"), # -- if not HLT re-run -- #
	L1Muon = cms.untracked.InputTag("hltGmtStage2Digis", "Muon"), # -- if HLT re-run --#

	TkMuon = cms.untracked.InputTag("hltHighPtTkMuonCands"), # -- for TkMu triggers

	lumiScaler = cms.untracked.InputTag("hltScalersRawToDigi"),
	offlineLumiScaler = cms.untracked.InputTag("scalersRawToDigi"),
	PUSummaryInfo = cms.untracked.InputTag("addPileupInfo"),
	genEventInfo = cms.untracked.InputTag("generator"),
	genParticle = cms.untracked.InputTag("genParticles"),

	iterL3OI        = cms.untracked.InputTag("hltL3MuonsIterL3OI", "", "MYHLT"),
	iterL3IOFromL2  = cms.untracked.InputTag("hltL3MuonsIterL3IO", "", "MYHLT"),
	iterL3FromL2    = cms.untracked.InputTag("hltIterL3MuonsFromL2LinksCombination", "", "MYHLT"),
	iterL3IOFromL1  = cms.untracked.InputTag("hltIter3IterL3FromL1MuonMerged", "", "MYHLT"),
	iterL3MuonNoID  = cms.untracked.InputTag("hltIterL3MuonsNoID", "", "MYHLT"),
)
