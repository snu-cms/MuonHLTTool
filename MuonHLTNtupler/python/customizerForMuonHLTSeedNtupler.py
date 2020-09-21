# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
# process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
import HLTrigger.MuonHLTSeedMVAClassifier.mvaScale as _mvaScale

def customizerFuncForMuonHLTSeedNtupler(process, newProcessName = "MYHLT", doDYSkim = False, MvaVersion = "Run3v0"):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput

    from MuonHLTTool.MuonHLTNtupler.ntupler_seed_cfi import seedNtuplerBase
    import SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi
    from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import tpClusterProducer as _tpClusterProducer

    process.hltTPClusterProducer = _tpClusterProducer.clone(
      pixelClusterSrc = "hltSiPixelClusters",
      stripClusterSrc = "hltSiStripRawToClustersFacility"
    )
    # process.hltTPClusterProducer.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel")
    process.hltTrackAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    process.hltTrackAssociatorByHits.cluster2TPSrc            = cms.InputTag("hltTPClusterProducer")
    process.hltTrackAssociatorByHits.UseGrouped               = cms.bool( False )
    process.hltTrackAssociatorByHits.UseSplitting             = cms.bool( False )
    process.hltTrackAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool( False )

    process.seedNtupler = seedNtuplerBase.clone()

    process.seedNtupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", newProcessName)
    # process.seedNtupler.L1Muon           = cms.untracked.InputTag("gmtStage2Digis",        "Muon", newProcessName)
    # process.seedNtupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", "HLT") #for phaseII w/o emulation
    # process.seedNtupler.L1Muon           = cms.untracked.InputTag("simGmtStage2Digis","",newProcessName)  # Phase II sim emul
    process.seedNtupler.L2Muon           = cms.untracked.InputTag("hltL2MuonCandidates",     "",     newProcessName)

    # process.seedNtupler.L1TkMuon                                          = cms.untracked.InputTag("L1TkMuons", "", newProcessName)
    # process.seedNtupler.L1TkPrimaryVertex                                 = cms.untracked.InputTag("L1TkPrimaryVertex", "", newProcessName)

    process.seedNtupler.hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltIterL3OISeedsFromL2Muons",                         "", newProcessName)
    process.seedNtupler.hltIter0IterL3MuonPixelSeedsFromPixelTracks       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",         "", newProcessName)
    process.seedNtupler.hltIter2IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeeds",                        "", newProcessName)
    process.seedNtupler.hltIter3IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeeds",                        "", newProcessName)
    process.seedNtupler.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",   "", newProcessName)
    process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeeds",                  "", newProcessName)
    process.seedNtupler.hltIter3IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeeds",                  "", newProcessName)
    
    process.seedNtupler.hltIterL3OIMuonTrack                              = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",             "", newProcessName)
    process.seedNtupler.hltIter0IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.seedNtupler.hltIter2IterL3MuonTrack                           = cms.untracked.InputTag("hltIter2IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.seedNtupler.hltIter3IterL3MuonTrack                           = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.seedNtupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.seedNtupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter2IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.seedNtupler.hltIter3IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)

    # process.seedNtupler.associatePixel = cms.bool(True)
    # process.seedNtupler.associateRecoTracks = cms.bool(False)
    # process.seedNtupler.associateStrip = cms.bool(True)
    # process.seedNtupler.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel")
    # process.seedNtupler.stripSimLinkSrc = cms.InputTag("simSiStripDigis")
    # process.seedNtupler.ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof', 'g4SimHitsTrackerHitsPixelBarrelHighTof', 'g4SimHitsTrackerHitsPixelEndcapLowTof', 'g4SimHitsTrackerHitsPixelEndcapHighTof')
    # process.seedNtupler.usePhase2Tracker = cms.bool(True)
    # process.seedNtupler.phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker")

    process.seedNtupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.seedNtupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")

    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_0                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIterL3OI_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_1                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIterL3OI_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_2                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIterL3OI_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_3                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIterL3OI_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_0       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_1       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_2       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_3       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0_3.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2_0.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2_1.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2_2.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0FromL1_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0FromL1_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0FromL1_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter0FromL1_3.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2FromL1_0.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2FromL1_1.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2FromL1_2.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter2FromL1_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3FromL1_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3FromL1_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3FromL1_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Barrel_hltIter3FromL1_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_0                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIterL3OI_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_1                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIterL3OI_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_2                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIterL3OI_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_3                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIterL3OI_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_0       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_1       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_2       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_3       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0_3.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2_0.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2_1.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2_2.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0FromL1_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0FromL1_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0FromL1_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter0FromL1_3.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2FromL1_0.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2FromL1_1.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2FromL1_2.xml" % MvaVersion)
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter2FromL1_3.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3FromL1_0.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3FromL1_1.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3FromL1_2.xml" % MvaVersion)
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%s_Endcap_hltIter3FromL1_3.xml" % MvaVersion)


    # process.seedNtupler.mvaScaleMeanHltIterL3OISeedsFromL2Muons_B                       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIterL3OI_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIterL3OISeedsFromL2Muons_B                        = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIterL3OI_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B        = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0_ScaleStd" % MvaVersion) )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2_ScaleMean" % MvaVersion) )
    process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter3IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0FromL1_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B  = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0FromL1_ScaleStd" % MvaVersion) )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2FromL1_ScaleMean" % MvaVersion) )
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2FromL1_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3FromL1_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3FromL1_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIterL3OISeedsFromL2Muons_E                       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIterL3OI_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIterL3OISeedsFromL2Muons_E                        = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIterL3OI_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E        = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0_ScaleStd" % MvaVersion) )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2_ScaleMean" % MvaVersion) )
    process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter3IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0FromL1_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E  = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0FromL1_ScaleStd" % MvaVersion) )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2FromL1_ScaleMean" % MvaVersion) )
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2FromL1_ScaleStd" % MvaVersion) )
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3FromL1_ScaleMean" % MvaVersion) )
    # process.seedNtupler.mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3FromL1_ScaleStd" % MvaVersion) )


    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    if doDYSkim:
        from MuonHLTTool.MuonHLTNtupler.DYmuSkimmer import DYmuSkimmer 
        process.Skimmer = DYmuSkimmer.clone()
        process.myseedpath = cms.Path(process.Skimmer*process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.seedNtupler)
    else:
        process.myseedpath = cms.Path(process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.seedNtupler)

    return process
