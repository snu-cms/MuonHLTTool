import FWCore.ParameterSet.Config as cms

MHPrintBTagMuTriggerInfo = cms.EDAnalyzer('MHPrintBTagMuTriggerInfo',
    triggerResults = cms.untracked.InputTag("TriggerResults::MYHLT"),
    triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
)


