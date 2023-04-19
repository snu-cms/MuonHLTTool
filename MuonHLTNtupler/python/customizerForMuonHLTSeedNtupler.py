# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
# process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
import HLTrigger.Configuration.MuonHLTForRun3.mvaScale as _mvaScale

def customizerFuncForMuonHLTSeedNtupler(process, newProcessName = "MYHLT", isDIGI = True):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput

    from MuonHLTTool.MuonHLTNtupler.ntupler_seed_cfi import seedNtuplerBase
    #import SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi
    import SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi
    from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import tpClusterProducer as _tpClusterProducer

    from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import simHitTPAssocProducer as _simHitTPAssocProducer
    process.simHitTPAssocProducer = _simHitTPAssocProducer.clone()

    # process.hltTPClusterProducer = _tpClusterProducer.clone(
    #   pixelClusterSrc = "hltSiPixelClusters",
    #   stripClusterSrc = "hltSiStripRawToClustersFacility"
    # )

    # process.hltTrackAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    # process.hltTrackAssociatorByHits.cluster2TPSrc            = cms.InputTag("hltTPClusterProducer")
    # process.hltTrackAssociatorByHits.UseGrouped               = cms.bool( False )
    # process.hltTrackAssociatorByHits.UseSplitting             = cms.bool( False )
    # process.hltTrackAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool( False )

    # process.hltSeedAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    # process.hltSeedAssociatorByHits.cluster2TPSrc            = cms.InputTag("hltTPClusterProducer")
    # process.hltSeedAssociatorByHits.UseGrouped               = cms.bool( False )
    # process.hltSeedAssociatorByHits.UseSplitting             = cms.bool( False )
    # process.hltSeedAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool( False )
    # process.hltSeedAssociatorByHits.Cut_RecoToSim = cms.double(0.)

    process.hltTrackAssociatorByHits = SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi.trackAssociatorByHits.clone(
        UsePixels = cms.bool(True),
        UseGrouped = cms.bool(True),
        UseSplitting = cms.bool(True),
        ThreeHitTracksAreSpecial = cms.bool(False),
        associatePixel = cms.bool(True),
        associateStrip = cms.bool(True),
        usePhase2Tracker = cms.bool(False),
        pixelSimLinkSrc = cms.InputTag("simSiPixelDigis"),
        stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
        phase2TrackerSimLinkSrc  = cms.InputTag("simSiPixelDigis","Tracker"),
        associateRecoTracks = cms.bool(True)
    )
    process.hltSeedAssociatorByHits = process.hltTrackAssociatorByHits.clone(
        Cut_RecoToSim = cms.double(0.)
    )

    process.seedNtupler = seedNtuplerBase.clone()

    process.seedNtupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", newProcessName)
    process.seedNtupler.L2Muon           = cms.untracked.InputTag("hltL2MuonCandidates",     "",     newProcessName)

    process.seedNtupler.hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltIterL3OISeedsFromL2Muons",                         "", newProcessName)
    process.seedNtupler.hltIter0IterL3MuonPixelSeedsFromPixelTracks       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",         "", newProcessName)
    process.seedNtupler.hltIter2IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",                        "", newProcessName)
    process.seedNtupler.hltIter3IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeeds",                        "", newProcessName)
    process.seedNtupler.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",   "", newProcessName)
    process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",                  "", newProcessName)
    process.seedNtupler.hltIter3IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeeds",                  "", newProcessName)
    
    process.seedNtupler.hltIterL3OIMuonTrack                              = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",             "", newProcessName)
    process.seedNtupler.hltIter0IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.seedNtupler.hltIter2IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.seedNtupler.hltIter3IterL3MuonTrack                           = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.seedNtupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.seedNtupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.seedNtupler.hltIter3IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)

    process.seedNtupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.seedNtupler.seedAssociator = cms.untracked.InputTag("hltSeedAssociatorByHits")
    process.seedNtupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")

    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_barrel_v3.xml")
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v3.xml")
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_endcap_v3.xml")
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v3.xml")

    process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_barrel_v3_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_barrel_v3_ScaleStd") )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v3_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v3_ScaleStd") )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_endcap_v3_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_endcap_v3_ScaleStd") )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v3_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v3_ScaleStd") )

    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    if isDIGI:
        process.myseedpath = cms.Path(process.HLTBeginSequence*
                                      process.HLTL2muonrecoSequence*
                                      process.HLTL3muonrecoSequence*
                                      # process.hltTPClusterProducer*
                                      process.simHitTPAssocProducer*
                                      process.hltTrackAssociatorByHits*
                                      process.hltSeedAssociatorByHits*
                                      process.seedNtupler)
    else:
        process.myseedpath = cms.Path(process.HLTBeginSequence*
                                      process.HLTL2muonrecoSequence*
                                      process.HLTL3muonrecoSequence*
                                      process.seedNtupler)

    return process
