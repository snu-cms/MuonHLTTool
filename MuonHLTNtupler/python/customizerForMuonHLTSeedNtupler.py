# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
# process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
import HLTrigger.Configuration.MuonHLTForRun3.mvaScale as _mvaScale

def customizerFuncForMuonHLTSeedNtupler(process, newProcessName = "MYHLT", doDYSkim = False, isDIGI = True, MvaVersion = "Run3v0"):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput

    from MuonHLTTool.MuonHLTNtupler.ntupler_seed_cfi import seedNtuplerBase
    import SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi
    from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import tpClusterProducer as _tpClusterProducer

    process.hltTPClusterProducer = _tpClusterProducer.clone(
      pixelClusterSrc = "hltSiPixelClusters",
      stripClusterSrc = "hltSiStripRawToClustersFacility"
    )

    process.hltTrackAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    process.hltTrackAssociatorByHits.cluster2TPSrc            = cms.InputTag("hltTPClusterProducer")
    process.hltTrackAssociatorByHits.UseGrouped               = cms.bool( False )
    process.hltTrackAssociatorByHits.UseSplitting             = cms.bool( False )
    process.hltTrackAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool( False )

    process.hltSeedAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    process.hltSeedAssociatorByHits.cluster2TPSrc            = cms.InputTag("hltTPClusterProducer")
    process.hltSeedAssociatorByHits.UseGrouped               = cms.bool( False )
    process.hltSeedAssociatorByHits.UseSplitting             = cms.bool( False )
    process.hltSeedAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool( False )
    process.hltSeedAssociatorByHits.Cut_RecoToSim = cms.double(0.)

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
    process.seedNtupler.hltIter2IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackWithVertexSelector",          "", newProcessName)
    process.seedNtupler.hltIter3IterL3MuonTrack                           = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.seedNtupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.seedNtupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackWithVertexSelector",    "", newProcessName)
    process.seedNtupler.hltIter3IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)

    process.seedNtupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.seedNtupler.seedAssociator = cms.untracked.InputTag("hltSeedAssociatorByHits")
    process.seedNtupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")

    #process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2Seeds_barrel.xml")
    #process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2FromL1Seeds_barrel.xml")
    #process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2Seeds_endcap.xml")
    #process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2FromL1Seeds_endcap.xml")

    #process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_barrel_ScaleMean") )
    #process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_barrel_ScaleStd") )
    #process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_barrel_ScaleMean") )
    #process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_barrel_ScaleStd") )
    #process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_endcap_ScaleMean") )
    #process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_endcap_ScaleStd") )
    #process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_endcap_ScaleMean") )
    #process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_endcap_ScaleStd") )
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Barrel_hltIter2.xml")
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Barrel_hltIter2FromL1.xml")
    process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Endcap_hltIter2.xml")
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Endcap_hltIter2FromL1.xml")

    process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2_ScaleStd") )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2FromL1_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2FromL1_ScaleStd") )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2_ScaleStd") )
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2FromL1_ScaleMean") )
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2FromL1_ScaleStd") )

    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    if doDYSkim:
        from MuonHLTTool.MuonHLTNtupler.DYmuSkimmer import DYmuSkimmer 
        process.Skimmer = DYmuSkimmer.clone()
        if isDIGI:
            process.myseedpath = cms.Path(process.Skimmer*
                                          process.HLTBeginSequence*
                                          process.HLTL2muonrecoSequence*
                                          process.HLTL3muonrecoSequence*
                                          process.hltTPClusterProducer*
                                          process.hltTrackAssociatorByHits*
                                          process.hltSeedAssociatorByHits*
                                          process.seedNtupler)
        else:
            process.myseedpath = cms.Path(process.Skimmer*
                                          process.HLTBeginSequence*
                                          process.HLTL2muonrecoSequence*
                                          process.HLTL3muonrecoSequence*
                                          process.seedNtupler)
    else:
        if isDIGI:
            process.myseedpath = cms.Path(process.HLTBeginSequence*
                                          process.HLTL2muonrecoSequence*
                                          process.HLTL3muonrecoSequence*
                                          process.hltTPClusterProducer*
                                          process.hltTrackAssociatorByHits*
                                          process.hltSeedAssociatorByHits*
                                          process.seedNtupler)
        else:
            process.myseedpath = cms.Path(process.HLTBeginSequence*
                                          process.HLTL2muonrecoSequence*
                                          process.HLTL3muonrecoSequence*
                                          process.seedNtupler)

    return process
