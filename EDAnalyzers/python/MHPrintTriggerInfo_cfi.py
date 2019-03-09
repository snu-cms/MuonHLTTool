import FWCore.ParameterSet.Config as cms

MHPrintTriggerInfo = cms.EDAnalyzer('MHPrintTriggerInfo',
    triggerResults = cms.untracked.InputTag("TriggerResults::HLT"),
    triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
)


