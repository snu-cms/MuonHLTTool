# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
# process = customizerFuncForMuonHLTNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
import HLTrigger.Configuration.MuonHLTForRun3.mvaScale as _mvaScale

def customizerFuncForMuonHLTNtupler(process, newProcessName = "MYHLT", isMC = False, isDIGI = True, MvaVersion = ""):
    process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput


    from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase
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

    process.hltGlbTrkMuonsTracks = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
    process.hltGlbTrkMuonsTracks.muonsTag                          = cms.InputTag("hltGlbTrkMuons", "", newProcessName)
    process.hltGlbTrkMuonsTracks.selectionTags                     = ('All',)
    process.hltGlbTrkMuonsTracks.trackType                         = "innerTrackPlusSegments"
    process.hltGlbTrkMuonsTracks.ignoreMissingMuonCollection       = True
    process.hltGlbTrkMuonsTracks.inputCSCSegmentCollection         = cms.InputTag("hltCscSegments", "", newProcessName)
    process.hltGlbTrkMuonsTracks.inputDTRecSegment4DCollection     = cms.InputTag("hltDt4DSegments", "", newProcessName)

    # Call the hit associator
    # https://github.com/cms-sw/cmssw/blob/CMSSW_12_1_0_pre4/Validation/RecoMuon/python/associators_cff.py
    from SimMuon.MCTruth.trackingParticleMuon_cfi import trackingParticleMuon as _trackingParticleMuon
    process.TPmu = _trackingParticleMuon.clone()

    from SimMuon.MCTruth.MuonAssociatorByHits_cfi import muonAssociatorByHits as _muonAssociatorByHits
    hltMuonAssociatorByHits = _muonAssociatorByHits.clone()
    # HERE
    hltMuonAssociatorByHits.PurityCut_track              = 0.75
    hltMuonAssociatorByHits.PurityCut_muon               = 0.75
    hltMuonAssociatorByHits.DTrechitTag                  = 'hltDt1DRecHits'
    hltMuonAssociatorByHits.ignoreMissingTrackCollection = True
    hltMuonAssociatorByHits.tpTag                        = ("TPmu")
    hltMuonAssociatorByHits.tpRefVector                  = True
    hltMuonAssociatorByHits.UseTracker                   = True
    hltMuonAssociatorByHits.UseMuon                      = False  # True

    # Hit associators from each track
    process.AhltIterL3OIMuonTrackSelectionHighPurity          = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIterL3OIMuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter0IterL3MuonTrackSelectionHighPurity       = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter2IterL3MuonTrackSelectionHighPurity       = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter0IterL3FromL1MuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter2IterL3FromL1MuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "", newProcessName) )
    process.AhltIter2IterL3MuonMerged                         = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter2IterL3MuonMerged", "", newProcessName) )
    process.AhltIter2IterL3FromL1MuonMerged                   = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIter2IterL3FromL1MuonMerged", "", newProcessName) )
    process.AhltIterL3MuonMerged                              = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIterL3MuonMerged", "", newProcessName) )
    process.AhltIterL3MuonAndMuonFromL1Merged                 = hltMuonAssociatorByHits.clone( tracksTag = cms.InputTag("hltIterL3MuonAndMuonFromL1Merged", "", newProcessName) )
    process.AhltIterL3MuonsNoID                               = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltIterL3MuonsNoIDTracks", "", newProcessName),
        UseMuon = cms.bool(True),
        rejectBadGlobal = cms.bool(False),
    )
    process.AhltIterL3Muons                                   = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltIterL3MuonsTracks", "", newProcessName),
        UseMuon = cms.bool(True),
        rejectBadGlobal = cms.bool(False),
    )
    process.AhltPixelTracks = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltPixelTracks", "", newProcessName),
        PurityCut_track = cms.double(0.65),
    )
    process.AhltPixelTracksInRegionL2 = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltPixelTracksInRegionL2", "", newProcessName),
        PurityCut_track = cms.double(0.65),
    )
    process.AhltPixelTracksInRegionL1 = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltPixelTracksInRegionL1", "", newProcessName),
        PurityCut_track = cms.double(0.65),
    )
    process.AhltPixelTracksForSeedsL3Muon = hltMuonAssociatorByHits.clone( # For Run2 Legacy comparison
        tracksTag = cms.InputTag("hltPixelTracksForSeedsL3Muon", "", newProcessName),
        PurityCut_track = cms.double(0.65),
    )
    process.AhltMuCtfTracks= hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltMuCtfTracks", "", newProcessName),
        PurityCut_track = cms.double(0.65),
    )
    process.AhltDiMuonMerging = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltDiMuonMerging", "", newProcessName),
        PurityCut_track = cms.double(0.65),
    )
    process.AhltGlbTrkMuons = hltMuonAssociatorByHits.clone(
        tracksTag = cms.InputTag("hltGlbTrkMuonsTracks", "", newProcessName),
        UseMuon = cms.bool(True),
        rejectBadGlobal = cms.bool(False),
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
        'iterL3MuonTrackAssociated',
        'hltPixelTracksAssociated',
        'hltPixelTracksInRegionL2Associated',
        'hltPixelTracksInRegionL1Associated',
        'hltPixelTracksForSeedsL3MuonAssociated',
        'hltMuCtfTracksAssociated',
        'hltDiMuonMergingAssociated',
        'hltGlbTrkMuonTracksAssociated',
    ]
    # Labels for calling each track in the Ntupler
    trackLabels = [
        cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",          "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",       "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",       "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "", newProcessName),
        cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "", newProcessName),
        cms.untracked.InputTag("hltIter2IterL3MuonMerged",                         "", newProcessName),
        cms.untracked.InputTag("hltIter2IterL3FromL1MuonMerged",                   "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonMerged",                              "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonAndMuonFromL1Merged",                 "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonsNoIDTracks",                         "", newProcessName),
        cms.untracked.InputTag("hltIterL3MuonsTracks",                             "", newProcessName),
        cms.untracked.InputTag("hltPixelTracks",                                   "", newProcessName),
        cms.untracked.InputTag("hltPixelTracksInRegionL2",                         "", newProcessName),
        cms.untracked.InputTag("hltPixelTracksInRegionL1",                         "", newProcessName),
        cms.untracked.InputTag("hltPixelTracksForSeedsL3Muon",                     "", newProcessName),
        cms.untracked.InputTag("hltMuCtfTracks",                                   "", newProcessName),
        cms.untracked.InputTag("hltDiMuonMerging",                                 "", newProcessName),
        cms.untracked.InputTag("hltGlbTrkMuonsTracks",                             "", newProcessName),
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
        'AhltIterL3Muons',
        'AhltPixelTracks',
        'AhltPixelTracksInRegionL2',
        'AhltPixelTracksInRegionL1',
        'AhltPixelTracksForSeedsL3Muon',
        'AhltMuCtfTracks',
        'AhltDiMuonMerging',
        'AhltGlbTrkMuons',
    ]
    process.trackAssoSeq = cms.Sequence(
        process.TPmu +
        process.hltIterL3MuonsNoIDTracks +
        process.hltIterL3MuonsTracks +
        process.hltGlbTrkMuonsTracks +
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
        process.AhltIterL3Muons +
        process.AhltPixelTracks +
        process.AhltPixelTracksInRegionL2 +
        process.AhltPixelTracksInRegionL1 +
        process.AhltPixelTracksForSeedsL3Muon +
        process.AhltMuCtfTracks +
        process.AhltDiMuonMerging +
        process.AhltGlbTrkMuons
    )

    # Call the L1 associators
    from MuonAnalysis.MuonAssociators.muonL1Match_cfi import muonL1Match as _muonL1Match
    process.recomuonL1Info = _muonL1Match.clone(
        src = cms.InputTag("muons"),
        useMB2InOverlap = cms.bool(True),
        useStage2L1 = cms.bool(True),
        preselection = cms.string(""),
        matched = cms.InputTag("hltGtStage2Digis:Muon:MYHLT"),
        useTrack = cms.string("none"),

        useStation2 = cms.bool(True),
        cosmicPropagationHypothesis = cms.bool(False),
        propagatorAlong = cms.ESInputTag("", "SteppingHelixPropagatorAlong"),
        propagatorAny = cms.ESInputTag("", "SteppingHelixPropagatorAny"),
        propagatorOpposite = cms.ESInputTag("", "SteppingHelixPropagatorOpposite"),
        fallbackToME1 = cms.bool(False)
    )
    process.recomuonL1InfoByQ = process.recomuonL1Info.clone(
        sortBy = cms.string("quality"),
        sortByQual     = cms.bool(True), # see MuonAnalysis/MuonAssociators/src/L1MuonMatcherAlgo.cc
        sortByDeltaPhi = cms.bool(False),
        sortByDeltaEta = cms.bool(False),
        sortByPt       = cms.bool(False)
    )
    process.L1AssoSeq = cms.Sequence(
        process.recomuonL1Info +
        process.recomuonL1InfoByQ
    )

    process.genmuonL1Info = _muonL1Match.clone()
    process.genmuonL1InfoByQ = process.genmuonL1Info.clone()
    if isMC:
        process.genmuonL1Info = _muonL1Match.clone(
            src = cms.InputTag("genParticles"),
            useMB2InOverlap = cms.bool(True),
            useStage2L1 = cms.bool(True),
            preselection = cms.string(""),
            matched = cms.InputTag("hltGtStage2Digis:Muon:MYHLT"),
            useTrack = cms.string("none"),

            useStation2 = cms.bool(True),
            cosmicPropagationHypothesis = cms.bool(False),
            propagatorAlong = cms.ESInputTag("", "SteppingHelixPropagatorAlong"),
            propagatorAny = cms.ESInputTag("", "SteppingHelixPropagatorAny"),
            propagatorOpposite = cms.ESInputTag("", "SteppingHelixPropagatorOpposite"),
            fallbackToME1 = cms.bool(False)
        )
        process.genmuonL1InfoByQ = process.genmuonL1Info.clone(
            sortBy = cms.string("quality"),
            sortByQual     = cms.bool(True), # see MuonAnalysis/MuonAssociators/src/L1MuonMatcherAlgo.cc
            sortByDeltaPhi = cms.bool(False),
            sortByDeltaEta = cms.bool(False),
            sortByPt       = cms.bool(False)
        )
        process.L1AssoSeq = cms.Sequence(
            process.recomuonL1Info +
            process.recomuonL1InfoByQ +
            process.genmuonL1Info +
            process.genmuonL1InfoByQ
        )

    process.ntupler = ntuplerBase.clone()

    # Add track hit association
    process.ntupler.trackCollectionNames  = cms.untracked.vstring(   trackNames )
    process.ntupler.trackCollectionLabels = cms.untracked.VInputTag( trackLabels )
    process.ntupler.associationLabels     = cms.untracked.VInputTag( assoLabels )

    # Add L1 association
    process.ntupler.recol1Matches = cms.InputTag("recomuonL1Info")
    process.ntupler.recol1MatchesQuality = cms.InputTag("recomuonL1Info", "quality")
    process.ntupler.recol1MatchesDeltaR = cms.InputTag("recomuonL1Info", "deltaR")
    process.ntupler.recol1MatchesByQ = cms.InputTag("recomuonL1InfoByQ")
    process.ntupler.recol1MatchesByQQuality = cms.InputTag("recomuonL1InfoByQ", "quality")
    process.ntupler.recol1MatchesByQDeltaR = cms.InputTag("recomuonL1InfoByQ", "deltaR")
    process.ntupler.genl1Matches = cms.InputTag("genmuonL1Info")
    process.ntupler.genl1MatchesQuality = cms.InputTag("genmuonL1Info", "quality")
    process.ntupler.genl1MatchesDeltaR = cms.InputTag("genmuonL1Info", "deltaR")
    process.ntupler.genl1MatchesByQ = cms.InputTag("genmuonL1InfoByQ")
    process.ntupler.genl1MatchesByQQuality = cms.InputTag("genmuonL1InfoByQ", "quality")
    process.ntupler.genl1MatchesByQDeltaR = cms.InputTag("genmuonL1InfoByQ", "deltaR")

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
    process.ntupler.hltIter2IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",         "", newProcessName)
    process.ntupler.hltIter3IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeeds",                        "", newProcessName)
    process.ntupler.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",   "", newProcessName)
    process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",   "", newProcessName)
    process.ntupler.hltIter3IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeeds",                  "", newProcessName)

    process.ntupler.hltIterL3OIMuonTrack                              = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",             "", newProcessName)
    process.ntupler.hltIter0IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.ntupler.hltIter2IterL3MuonTrack                           = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.ntupler.hltIter3IterL3MuonTrack                           = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",          "", newProcessName)
    process.ntupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.ntupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)
    process.ntupler.hltIter3IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",    "", newProcessName)

    process.ntupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.ntupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")

    process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_barrel_v2.xml")
    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v2.xml")
    process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_endcap_v2.xml")
    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v2.xml")

    process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_barrel_v2_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_barrel_v2_ScaleStd") )
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v2_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v2_ScaleStd") )
    process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_endcap_v2_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0_PatatrackSeeds_endcap_v2_ScaleStd") )
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v2_ScaleMean") )
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, "xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v2_ScaleStd") )

    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_B                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/"+MvaVersion+"_Barrel_hltIter2.xml")
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/"+MvaVersion+"_Barrel_hltIter2FromL1.xml")
    #process.ntupler.mvaFileHltIter2IterL3MuonPixelSeeds_E                      = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/"+MvaVersion+"_Endcap_hltIter2.xml")
    #process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/"+MvaVersion+"_Endcap_hltIter2FromL1.xml")

    #process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Barrel_hltIter2_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Barrel_hltIter2_ScaleStd") )
    #process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Barrel_hltIter2FromL1_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Barrel_hltIter2FromL1_ScaleStd") )
    #process.ntupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Endcap_hltIter2_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Endcap_hltIter2_ScaleStd") )
    #process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Endcap_hltIter2FromL1_ScaleMean") )
    #process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble( getattr(_mvaScale, MvaVersion+"_Endcap_hltIter2FromL1_ScaleStd") )

    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    # process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults::HLT") # dummy to avoid ordering error occur in skimming, as it is not used at the moment
    process.ntupler.DebugMode = cms.bool(False)

    if isDIGI:
        process.mypath = cms.Path(process.HLTBeginSequence*
                                  process.HLTL2muonrecoSequence*
                                  process.HLTL3muonrecoSequence*
                                  # process.hltTPClusterProducer*
                                  process.simHitTPAssocProducer*
                                  process.hltTrackAssociatorByHits*
                                  process.trackAssoSeq*
                                  process.L1AssoSeq)
        process.myendpath = cms.EndPath(process.ntupler)
    else:
        process.mypath = cms.Path(process.HLTBeginSequence*
                                  process.HLTL2muonrecoSequence*
                                  process.HLTL3muonrecoSequence*
                                  process.L1AssoSeq)
        process.myendpath = cms.EndPath(process.ntupler)

    return process
