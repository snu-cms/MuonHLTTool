# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
# process = customizerFuncForMuonHLTNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
import HLTrigger.Configuration.MuonHLTForRun3.mvaScale as _mvaScale

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

    process.hltTrackAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    process.hltTrackAssociatorByHits.cluster2TPSrc            = cms.InputTag("hltTPClusterProducer")
    process.hltTrackAssociatorByHits.UseGrouped               = cms.bool( False )
    process.hltTrackAssociatorByHits.UseSplitting             = cms.bool( False )
    process.hltTrackAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool( False )

    # Produce tracks from L3 muons -- to use track hit association
    import SimMuon.MCTruth.MuonTrackProducer_cfi
    process.hltIterL3MuonsNoIDTracks = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
    process.hltIterL3MuonsNoIDTracks.muonsTag                      = cms.InputTag("hltIterL3MuonsNoID", "", newProcessName)
    process.hltIterL3MuonsNoIDTracks.selectionTags                 = ('All',)
    process.hltIterL3MuonsNoIDTracks.trackType                     = "innerTrackPlusSegments"
    process.hltIterL3MuonsNoIDTracks.ignoreMissingMuonCollection   = True
    process.hltIterL3MuonsNoIDTracks.inputCSCSegmentCollection     = cms.InputTag("hltCscSegments", "", newProcessName)
    process.hltIterL3MuonsNoIDTracks.inputDTRecSegment4DCollection = cms.InputTag("hltDt4DSegments", "", newProcessName)

    process.hltIterL3MuonsTracks = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
    process.hltIterL3MuonsTracks.muonsTag                          = cms.InputTag("hltIterL3Muons", "", newProcessName)
    process.hltIterL3MuonsTracks.selectionTags                     = ('All',)
    process.hltIterL3MuonsTracks.trackType                         = "innerTrackPlusSegments"
    process.hltIterL3MuonsTracks.ignoreMissingMuonCollection       = True
    process.hltIterL3MuonsTracks.inputCSCSegmentCollection         = cms.InputTag("hltCscSegments", "", newProcessName)
    process.hltIterL3MuonsTracks.inputDTRecSegment4DCollection     = cms.InputTag("hltDt4DSegments", "", newProcessName)

    # Call the hit associator
    from SimMuon.MCTruth.MuonAssociatorByHits_cfi import muonAssociatorByHits as _muonAssociatorByHits
    hltMuonAssociatorByHits = _muonAssociatorByHits.clone()
    hltMuonAssociatorByHits.PurityCut_track              = 0.75
    hltMuonAssociatorByHits.PurityCut_muon               = 0.75
    hltMuonAssociatorByHits.DTrechitTag                  = cms.InputTag("hltDt1DRecHits", "", newProcessName)
    hltMuonAssociatorByHits.ignoreMissingTrackCollection = True
    hltMuonAssociatorByHits.UseTracker                   = True
    hltMuonAssociatorByHits.UseMuon                      = False  # True

    # Hit associators from each track
    process.AhltIterL3OIMuonTrackSelectionHighPurity          = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIterL3OIMuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter0IterL3MuonTrackSelectionHighPurity       = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter2IterL3MuonTrackSelectionHighPurity       = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3MuonTrackWithVertexSelector", "", newProcessName) )
    process.AhltIter0IterL3FromL1MuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter2IterL3FromL1MuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3FromL1MuonTrackWithVertexSelector", "", newProcessName) )
    process.AhltIter2IterL3MuonMerged                         = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter2IterL3MuonMerged", "", newProcessName) )
    process.AhltIter2IterL3FromL1MuonMerged                   = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter2IterL3FromL1MuonMerged", "", newProcessName) )
    process.AhltIterL3MuonMerged                              = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIterL3MuonMerged", "", newProcessName) )
    process.AhltIterL3MuonAndMuonFromL1Merged                 = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIterL3MuonAndMuonFromL1Merged", "", newProcessName) )
    process.AhltIterL3MuonsNoID                               = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltIterL3MuonsNoIDTracks", "", newProcessName),
        UseMuon = cms.bool(True),
    )
    process.AhltIterL3Muons                                   = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltIterL3MuonsTracks", "", newProcessName),
        UseMuon = cms.bool(True),
    )

    # Names for ntuple variable
    trackNames = [
        'hltIterL3OIMuonTrackAssociated',
        'hltIter0IterL3MuonTrackAssociated',
        'hltIter2IterL3MuonTrackAssociated',
        'hltIter0IterL3FromL1MuonTrackAssociated',
        'hltIter2IterL3FromL1MuonTrackAssociated',
        'hltIter2IterL3MuonMergedAssociated',
        'hltIter2IterL3FromL1MuonMergedAssociated',
        'hltIterL3MuonMergedAssociated',
        'hltIterL3MuonAndMuonFromL1MergedAssociated',
        'iterL3MuonNoIDTrackAssociated',
        'iterL3MuonTrackAssociated'
    ]
    # Labels for calling each track in the Ntupler
    trackLabels = [
        cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",          "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",       "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3MuonTrackWithVertexSelector",       "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackWithVertexSelector", "", newProcessName),
        cms.untracked.InputTag("hltIter2IterL3MuonMerged",                         "", newProcessName),
        cms.untracked.InputTag("hltIter2IterL3FromL1MuonMerged",                   "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonMerged",                              "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonAndMuonFromL1Merged",                 "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonsNoIDTracks",                         "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonsTracks",                             "", newProcessName)
    ]
    # Labels for calling each associator in the Ntupler
    assoLabels = [
        'AhltIterL3OIMuonTrackSelectionHighPurity',
        'AhltIter0IterL3MuonTrackSelectionHighPurity',
        'AhltIter2IterL3MuonTrackSelectionHighPurity',
        'AhltIter0IterL3FromL1MuonTrackSelectionHighPurity',
        'AhltIter2IterL3FromL1MuonTrackSelectionHighPurity',
        'AhltIter2IterL3MuonMerged',
        'AhltIter2IterL3FromL1MuonMerged',
        'AhltIterL3MuonMerged',
        'AhltIterL3MuonAndMuonFromL1Merged',
        'AhltIterL3MuonsNoID',
        'AhltIterL3Muons'
    ]
    process.trackAssoSeq = cms.Sequence(
        process.hltIterL3MuonsNoIDTracks +
        process.hltIterL3MuonsTracks +
        process.AhltIterL3OIMuonTrackSelectionHighPurity +
        process.AhltIter0IterL3MuonTrackSelectionHighPurity +
        process.AhltIter2IterL3MuonTrackSelectionHighPurity +
        process.AhltIter0IterL3FromL1MuonTrackSelectionHighPurity +
        process.AhltIter2IterL3FromL1MuonTrackSelectionHighPurity +
        process.AhltIter2IterL3MuonMerged +
        process.AhltIter2IterL3FromL1MuonMerged +
        process.AhltIterL3MuonMerged +
        process.AhltIterL3MuonAndMuonFromL1Merged +
        process.AhltIterL3MuonsNoID +
        process.AhltIterL3Muons
    )

    # Call the L1 associators
    from MuonAnalysis.MuonAssociators.muonL1Match_cfi import muonL1Match as _muonL1Match
    process.muonL1Info = _muonL1Match.clone(
        src = cms.InputTag("genParticles"),
        useMB2InOverlap = cms.bool(True),
        useStage2L1 = cms.bool(True),
        preselection = cms.string(""),
        matched = cms.InputTag("hltGtStage2Digis:Muon:MYHLT"),
        useTrack = cms.string("none")
    )
    process.muonL1InfoByQ = process.muonL1Info.clone(
        sortBy = cms.string("quality"),
        sortByQual     = cms.bool(True), # see MuonAnalysis/MuonAssociators/src/L1MuonMatcherAlgo.cc
        sortByDeltaPhi = cms.bool(False),
        sortByDeltaEta = cms.bool(False),
        sortByPt       = cms.bool(False)
    )
    process.L1AssoSeq = cms.Sequence(
        process.muonL1Info +
        process.muonL1InfoByQ
    )

    process.ntupler = ntuplerBase.clone()

    # Add track hit association
    process.ntupler.trackCollectionNames  = cms.untracked.vstring(   trackNames )
    process.ntupler.trackCollectionLabels = cms.untracked.VInputTag( trackLabels )
    process.ntupler.associationLabels     = cms.untracked.VInputTag( assoLabels )

    # Add L1 association
    process.ntupler.l1Matches = cms.InputTag("muonL1Info")
    process.ntupler.l1MatchesQuality = cms.InputTag("muonL1Info", "quality")
    process.ntupler.l1MatchesDeltaR = cms.InputTag("muonL1Info", "deltaR")
    process.ntupler.l1MatchesByQ = cms.InputTag("muonL1InfoByQ")
    process.ntupler.l1MatchesByQQuality = cms.InputTag("muonL1InfoByQ", "quality")
    process.ntupler.l1MatchesByQDeltaR = cms.InputTag("muonL1InfoByQ", "deltaR")

    # -- set to the new process name
    process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults",          "",     newProcessName)
    process.ntupler.myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD",    "",     newProcessName)
    process.ntupler.lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi",     "",     newProcessName)

    process.ntupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", newProcessName)
    process.ntupler.L2Muon           = cms.untracked.InputTag("hltL2MuonCandidates",     "",     newProcessName)
    process.ntupler.L3Muon           = cms.untracked.InputTag("hltIterL3MuonCandidates", "",     newProcessName)
    process.ntupler.TkMuon           = cms.untracked.InputTag("hltHighPtTkMuonCands",    "",     newProcessName)

    process.ntupler.iterL3OI         = cms.untracked.InputTag("hltL3MuonsIterL3OI",                   "", newProcessName)
    process.ntupler.iterL3IOFromL2   = cms.untracked.InputTag("hltL3MuonsIterL3IO",                   "", newProcessName)
    process.ntupler.iterL3FromL2     = cms.untracked.InputTag("hltIterL3MuonsFromL2LinksCombination", "", newProcessName)
    process.ntupler.iterL3IOFromL1   = cms.untracked.InputTag("hltIter2IterL3FromL1MuonMerged",       "", newProcessName)
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
    process.ntupler.hltIter2IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackWithVertexSelector",          "", newProcessName)
    process.ntupler.hltIter3IterL3MuonTrack                           = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.ntupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.ntupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackWithVertexSelector",    "", newProcessName)
    process.ntupler.hltIter3IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)

    process.ntupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.ntupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")

    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2Seeds_barrel.xml")
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2FromL1Seeds_barrel.xml")
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2Seeds_endcap.xml")
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter2FromL1Seeds_endcap.xml")

    #process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_barrel_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_barrel_ScaleStd") )
    #process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_barrel_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_barrel_ScaleStd") )
    #process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_endcap_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2Seeds_endcap_ScaleStd") )
    #process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_endcap_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter2FromL1Seeds_endcap_ScaleStd") )

    process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Barrel_hltIter2.xml")
    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Barrel_hltIter2FromL1.xml")
    process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Endcap_hltIter2.xml")
    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/v7Fast_Endcap_hltIter2FromL1.xml")

    process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2_ScaleStd") )
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2FromL1_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "v7Fast_Barrel_hltIter2FromL1_ScaleStd") )
    process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2_ScaleStd") )
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2FromL1_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "v7Fast_Endcap_hltIter2FromL1_ScaleStd") )

    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults::HLT") # dummy to avoid ordering error occur in skimming, as it is not used at the moment
    process.ntupler.DebugMode = cms.bool(False)

    if doDYSkim:
        from MuonHLTTool.MuonHLTNtupler.DYmuSkimmer import DYmuSkimmer 
        process.Skimmer = DYmuSkimmer.clone()
        if isDIGI:
            process.mypath = cms.Path(process.Skimmer*
                                      process.HLTBeginSequence*
                                      process.HLTL2muonrecoSequence*
                                      process.HLTL3muonrecoSequence*
                                      process.hltTPClusterProducer*
                                      process.hltTrackAssociatorByHits*
                                      process.trackAssoSeq*
                                      process.L1AssoSeq*
                                      process.ntupler)
        else:
            process.mypath = cms.Path(process.Skimmer*
                                      process.HLTBeginSequence*
                                      process.HLTL2muonrecoSequence*
                                      process.HLTL3muonrecoSequence*
                                      process.L1AssoSeq*
                                      process.ntupler)
    else:
        if isDIGI:
            process.mypath = cms.Path(process.HLTBeginSequence*
                                      process.HLTL2muonrecoSequence*
                                      process.HLTL3muonrecoSequence*
                                      process.hltTPClusterProducer*
                                      process.hltTrackAssociatorByHits*
                                      process.trackAssoSeq*
                                      process.L1AssoSeq*
                                      process.ntupler)
        else:
            process.mypath = cms.Path(process.HLTBeginSequence*
                                      process.HLTL2muonrecoSequence*
                                      process.HLTL3muonrecoSequence*
                                      process.L1AssoSeq*
                                      process.ntupler)

    return process
