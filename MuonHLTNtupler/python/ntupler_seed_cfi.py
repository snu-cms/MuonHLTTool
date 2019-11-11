import FWCore.ParameterSet.Config as cms

seedNtuplerBase = cms.EDAnalyzer("MuonHLTSeedNtupler",

	L1Muon = cms.untracked.InputTag("hltGmtStage2Digis",       "Muon", "MYHLT"), # -- if L1 emulator is used

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
)
