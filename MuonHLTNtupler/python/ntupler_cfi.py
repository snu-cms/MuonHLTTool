import FWCore.ParameterSet.Config as cms

ntuplerBase = cms.EDAnalyzer("MuonHLTNtupler",
	# -- information stored in edm file
	triggerResults    = cms.untracked.InputTag("TriggerResults::HLT"),
	triggerEvent      = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
	offlineLumiScaler = cms.untracked.InputTag("scalersRawToDigi"),
	offlineVertex     = cms.untracked.InputTag("offlinePrimaryVertices"),
	offlineMuon       = cms.untracked.InputTag("muons"),
	beamSpot          = cms.untracked.InputTag("hltOnlineBeamSpot"),

	# -- newly created objects by HLT rerun
	# -- new process name = "MYHLT"
	myTriggerResults = cms.untracked.InputTag("TriggerResults",       "", "MYHLT"),
	myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "MYHLT"),
	lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi",  "", "MYHLT"),

	L1Muon = cms.untracked.InputTag("hltGmtStage2Digis",       "Muon", "MYHLT"), # -- if L1 emulator is used
	# L1Muon = cms.untracked.InputTag("gmtStage2Digis",          "Muon", "RECO"), # -- if L1 is not emulated
	L2Muon = cms.untracked.InputTag("hltL2MuonCandidates",     "",     "MYHLT"),
	L3Muon = cms.untracked.InputTag("hltIterL3MuonCandidates", "",     "MYHLT"),
	TkMuon = cms.untracked.InputTag("hltHighPtTkMuonCands",    "",     "MYHLT"),

	ECALIsoMap = cms.untracked.InputTag("hltMuonEcalMFPFClusterIsoForMuons", "",               "MYHLT"),
	HCALIsoMap = cms.untracked.InputTag("hltMuonHcalRegPFClusterIsoForMuons", "" ,              "MYHLT"),
	trkIsoMap  = cms.untracked.InputTag("hltMuonTkRelIsolationCut0p08Map",               "trkIsoDeposits", "MYHLT"),

	rho_ECAL = cms.untracked.InputTag("hltFixedGridRhoFastjetECALMFForMuons", "", "MYHLT"),
	rho_HCAL = cms.untracked.InputTag("hltFixedGridRhoFastjetHCAL",           "", "MYHLT"),

	iterL3OI        = cms.untracked.InputTag("hltL3MuonsIterL3OI",                   "", "MYHLT"),
	iterL3IOFromL2  = cms.untracked.InputTag("hltL3MuonsIterL3IO",                   "", "MYHLT"),
	iterL3FromL2    = cms.untracked.InputTag("hltIterL3MuonsFromL2LinksCombination", "", "MYHLT"),
	iterL3IOFromL1  = cms.untracked.InputTag("hltIter3IterL3FromL1MuonMerged",       "", "MYHLT"),
	iterL3MuonNoID  = cms.untracked.InputTag("hltIterL3MuonsNoID",                   "", "MYHLT"),
	iterL3Muon      = cms.untracked.InputTag("hltIterL3Muons",                       "", "MYHLT"),

	hltIterL3MuonTrimmedPixelVertices                 = cms.untracked.InputTag("hltIterL3MuonTrimmedPixelVertices",                   "", "MYHLT"),
	hltIterL3FromL1MuonTrimmedPixelVertices           = cms.untracked.InputTag("hltIterL3FromL1MuonTrimmedPixelVertices",             "", "MYHLT"),

	doMVA  = cms.bool(True),
	doSeed = cms.bool(True),

	hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltIterL3OISeedsFromL2Muons",                         "", "MYHLT"),
	hltIter0IterL3MuonPixelSeedsFromPixelTracks       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",         "", "MYHLT"),
	hltIter2IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeeds",                        "", "MYHLT"),
	hltIter3IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeeds",                        "", "MYHLT"),
	hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",   "", "MYHLT"),
	hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeeds",                  "", "MYHLT"),
	hltIter3IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeeds",                  "", "MYHLT"),

	hltIterL3OIMuonTrack          = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",       "", "MYHLT"),
	hltIter0IterL3MuonTrack       = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",       "", "MYHLT"),
	hltIter2IterL3MuonTrack       = cms.untracked.InputTag("hltIter2IterL3MuonTrackSelectionHighPurity",       "", "MYHLT"),
	hltIter3IterL3MuonTrack       = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",       "", "MYHLT"), 
	hltIter0IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",       "", "MYHLT"),
	hltIter2IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter2IterL3FromL1MuonTrackSelectionHighPurity",       "", "MYHLT"),
	hltIter3IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",       "", "MYHLT"),

	# -- generator information
	PUSummaryInfo = cms.untracked.InputTag("addPileupInfo"),
	genEventInfo = cms.untracked.InputTag("generator"),
	genParticle = cms.untracked.InputTag("genParticles"),
)
