# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForNtupler import *
# customizerFuncForMuonHLTNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms

def customizerFuncForMuonHLTNtupler(process, newProcessName = "MYHLT"):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput

    flag_HLTRerun = True

    from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase

    process.ntupler = ntuplerBase.clone()
    process.ntupler.offlineMuon = cms.untracked.InputTag("muons")
    process.ntupler.L3Muon = cms.untracked.InputTag("hltIterL3MuonCandidates")
    process.ntupler.L2Muon = cms.untracked.InputTag("hltL2MuonCandidates")
    process.ntupler.triggerResults = cms.untracked.InputTag("TriggerResults", "", "HLT")
    process.ntupler.triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT")
    process.ntupler.offlineLumiScaler = cms.untracked.InputTag("scalersRawToDigi")

    if flag_HLTRerun: # -- after HLT re-run -- #
        process.ntupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis", "Muon", newProcessName)
        process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults", "", newProcessName) # -- result after rerun HLT -- #
        process.ntupler.myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD", "", newProcessName) # -- result after rerun HLT -- #
        process.ntupler.lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi", "", newProcessName)
        process.ntupler.iterL3MuonNoID   = cms.untracked.InputTag("hltIterL3MuonsNoID", "", newProcessName)
    else: # -- without HLT re-run -- #
        process.ntupler.L1Muon = cms.untracked.InputTag("gmtStage2Digis", "Muon", "RECO")
        process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults")
        process.ntupler.myTriggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD")
        process.ntupler.lumiScaler = cms.untracked.InputTag("scalersRawToDigi") # -- same with OfflineLumiScaler

    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ntuple.root"),
      closeFileFast = cms.untracked.bool(False),
      )

    process.mypath = cms.EndPath(process.ntupler)


