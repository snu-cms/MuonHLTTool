# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
# process = customizerFuncForMuonHLTNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
import RecoMuon.TrackerSeedGenerator.mvaScale as _mvaScale

def customizerFuncForMuonHLTNtupler(process, newProcessName = "MYHLT", doDYSkim = False, isDIGI = True, MvaVersion = "Run3v0"):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput


    from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase
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
    process.ntupler = ntuplerBase.clone()

    # -- set to the new process name
    process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults",          "",     newProcessName)
    process.ntupler.myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD",    "",     newProcessName)
    process.ntupler.lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi",     "",     newProcessName)

    process.ntupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", newProcessName)
    # process.ntupler.L1Muon           = cms.untracked.InputTag("gmtStage2Digis",        "Muon", newProcessName) 
    # process.ntupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", "HLT") #for phaseII w/o emulation
    # process.ntupler.L1Muon           = cms.untracked.InputTag("simGmtStage2Digis","",newProcessName)  # Phase II sim emul
    process.ntupler.L2Muon           = cms.untracked.InputTag("hltL2MuonCandidates",     "",     newProcessName)
    process.ntupler.L3Muon           = cms.untracked.InputTag("hltIterL3MuonCandidates", "",     newProcessName)
    process.ntupler.TkMuon           = cms.untracked.InputTag("hltHighPtTkMuonCands",    "",     newProcessName)

    process.ntupler.iterL3OI         = cms.untracked.InputTag("hltL3MuonsIterL3OI",                   "", newProcessName)
    process.ntupler.iterL3IOFromL2   = cms.untracked.InputTag("hltL3MuonsIterL3IO",                   "", newProcessName)
    process.ntupler.iterL3FromL2     = cms.untracked.InputTag("hltIterL3MuonsFromL2LinksCombination", "", newProcessName)
    process.ntupler.iterL3IOFromL1   = cms.untracked.InputTag("hltIter3IterL3FromL1MuonMerged",       "", newProcessName)
    process.ntupler.iterL3MuonNoID   = cms.untracked.InputTag("hltIterL3MuonsNoID",                   "", newProcessName)
    process.ntupler.iterL3Muon       = cms.untracked.InputTag("hltIterL3Muons",                       "", newProcessName)

    process.ntupler.hltIterL3MuonTrimmedPixelVertices                 = cms.untracked.InputTag("hltIterL3MuonTrimmedPixelVertices",                   "", newProcessName)
    process.ntupler.hltIterL3FromL1MuonTrimmedPixelVertices           = cms.untracked.InputTag("hltIterL3FromL1MuonTrimmedPixelVertices",             "", newProcessName)

    process.ntupler.doMVA  = cms.bool(True)
    process.ntupler.doSeed = cms.bool(False)

    process.ntupler.hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltIterL3OISeedsFromL2Muons",                         "", newProcessName)
    process.ntupler.hltIter0IterL3MuonPixelSeedsFromPixelTracks       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",         "", newProcessName)
    process.ntupler.hltIter2IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeeds",                        "", newProcessName)
    process.ntupler.hltIter3IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeeds",                        "", newProcessName)
    process.ntupler.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",   "", newProcessName)
    process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeeds",                  "", newProcessName)
    process.ntupler.hltIter3IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeeds",                  "", newProcessName)

    process.ntupler.hltIterL3OIMuonTrack                              = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",             "", newProcessName)
    process.ntupler.hltIter0IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.ntupler.hltIter2IterL3MuonTrack                           = cms.untracked.InputTag("hltIter2IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.ntupler.hltIter3IterL3MuonTrack                           = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.ntupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.ntupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter2IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.ntupler.hltIter3IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)

    # process.ntupler.associatePixel = cms.bool(True)
    # process.ntupler.associateRecoTracks = cms.bool(False)
    # process.ntupler.associateStrip = cms.bool(True)
    # process.ntupler.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel")
    # process.ntupler.stripSimLinkSrc = cms.InputTag("simSiStripDigis")
    # process.ntupler.ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof', 'g4SimHitsTrackerHitsPixelBarrelHighTof', 'g4SimHitsTrackerHitsPixelEndcapLowTof', 'g4SimHitsTrackerHitsPixelEndcapHighTof')
    # process.ntupler.usePhase2Tracker = cms.bool(True)
    # process.ntupler.phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker")

    process.ntupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.ntupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")


    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_B_0                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIterL3OI_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_B_1                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIterL3OI_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_B_2                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIterL3OI_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_B_3                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIterL3OI_3.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_0       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0_0.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_1       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0_1.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_2       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0_2.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_3       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0_3.xml" % MvaVersion)
    process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_0                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_1                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_2                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_3                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2_3.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_0                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_1                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_2                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_3                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3_3.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_0 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0FromL1_0.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_1 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0FromL1_1.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_2 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0FromL1_2.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_3 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter0FromL1_3.xml" % MvaVersion)
    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2FromL1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2FromL1_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_1                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2FromL1_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_2                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2FromL1_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_3                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter2FromL1_3.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3FromL1_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_1                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3FromL1_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_2                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3FromL1_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_3                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Barrel_hltIter3FromL1_3.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_E_0                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIterL3OI_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_E_1                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIterL3OI_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_E_2                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIterL3OI_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIterL3OISeedsFromL2Muons_E_3                       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIterL3OI_3.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_0       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0_0.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_1       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0_1.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_2       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0_2.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_3       = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0_3.xml" % MvaVersion)
    process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_0                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_1                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_2                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_3                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2_3.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_0                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_1                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_2                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_3                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3_3.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_0 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0FromL1_0.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_1 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0FromL1_1.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_2 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0FromL1_2.xml" % MvaVersion)
    # process.ntupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_3 = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter0FromL1_3.xml" % MvaVersion)
    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2FromL1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2FromL1_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_1                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2FromL1_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_2                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2FromL1_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_3                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter2FromL1_3.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3FromL1_0.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_1                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3FromL1_1.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_2                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3FromL1_2.xml" % MvaVersion)
    #process.ntupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_3                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/%s_Endcap_hltIter3FromL1_3.xml" % MvaVersion)


    #process.ntupler.mvaScaleMeanHltIterL3OISeedsFromL2Muons_B                       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIterL3OI_ScaleMean" % MvaVersion) )
    #process.ntupler.mvaScaleStdHltIterL3OISeedsFromL2Muons_B                        = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIterL3OI_ScaleStd" % MvaVersion) )
    # process.ntupler.mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0_ScaleMean" % MvaVersion) )
    # process.ntupler.mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B        = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0_ScaleStd" % MvaVersion) )
    process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2_ScaleMean" % MvaVersion) )
    process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2_ScaleStd" % MvaVersion) )
    #process.ntupler.mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3_ScaleMean" % MvaVersion) )
    #process.ntupler.mvaScaleStdHltIter3IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3_ScaleStd" % MvaVersion) )
    # process.ntupler.mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0FromL1_ScaleMean" % MvaVersion) )
    # process.ntupler.mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B  = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter0FromL1_ScaleStd" % MvaVersion) )
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2FromL1_ScaleMean" % MvaVersion) )
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter2FromL1_ScaleStd" % MvaVersion) )
    #process.ntupler.mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3FromL1_ScaleMean" % MvaVersion) )
    #process.ntupler.mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "%s_Barrel_hltIter3FromL1_ScaleStd" % MvaVersion) )
    #process.ntupler.mvaScaleMeanHltIterL3OISeedsFromL2Muons_E                       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIterL3OI_ScaleMean" % MvaVersion) )
    #process.ntupler.mvaScaleStdHltIterL3OISeedsFromL2Muons_E                        = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIterL3OI_ScaleStd" % MvaVersion) )
    # process.ntupler.mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0_ScaleMean" % MvaVersion) )
    # process.ntupler.mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E        = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0_ScaleStd" % MvaVersion) )
    process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2_ScaleMean" % MvaVersion) )
    process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2_ScaleStd" % MvaVersion) )
    #process.ntupler.mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3_ScaleMean" % MvaVersion) )
    #process.ntupler.mvaScaleStdHltIter3IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3_ScaleStd" % MvaVersion) )
    # process.ntupler.mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0FromL1_ScaleMean" % MvaVersion) )
    # process.ntupler.mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E  = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter0FromL1_ScaleStd" % MvaVersion) )
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2FromL1_ScaleMean" % MvaVersion) )
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter2FromL1_ScaleStd" % MvaVersion) )
    #process.ntupler.mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3FromL1_ScaleMean" % MvaVersion) )
    #process.ntupler.mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "%s_Endcap_hltIter3FromL1_ScaleStd" % MvaVersion) )



    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults::HLT") # dummy to avoid ordering error occur in skimming, as it is not used at the moment

    # L1TRK_PROC  =  process.TTTracksFromTrackletEmulation
    # L1TRK_NAME  = "TTTracksFromTrackletEmulation"
    # L1TRK_LABEL = "Level1TTTracks"
    # process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks") )

    process.ntupler.DebugMode = cms.bool(False)
    # process.ntupler.SaveAllTracks = cms.bool(True)
    # process.ntupler.SaveStubs = cms.bool(False)
    # process.ntupler.L1TrackInputTag = cms.InputTag(L1TRK_NAME, L1TRK_LABEL) # TTTrack input 
    # process.ntupler.MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", L1TRK_LABEL)  ## MCTruth input
    # process.ntupler.L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted")
    # process.ntupler.TkMuonToken = cms.InputTag("L1TkMuons","")
    # process.ntupler.l1TkPrimaryVertex = cms.InputTag("L1TkPrimaryVertex","")

    if doDYSkim:
        from MuonHLTTool.MuonHLTNtupler.DYmuSkimmer import DYmuSkimmer 
        process.Skimmer = DYmuSkimmer.clone()
        if isDIGI:
            process.mypath = cms.Path(process.Skimmer*process.HLTBeginSequence*process.HLTL2muonrecoSequence*process.HLTL3muonrecoSequence*process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.ntupler)
        else:
            process.mypath = cms.Path(process.Skimmer*process.ntupler)
    else:
        if isDIGI:
            process.mypath = cms.Path(process.HLTBeginSequence*process.HLTL2muonrecoSequence*process.HLTL3muonrecoSequence*process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.ntupler)
        else:
            process.mypath = cms.Path(process.ntupler)

    return process
