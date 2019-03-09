import FWCore.ParameterSet.Config as cms

MHFindTriggerFiredEvents = cms.EDAnalyzer('MHFindTriggerFiredEvents',
    triggerResults = cms.untracked.InputTag("TriggerResults::HLT"),
    trigName = cms.untracked.string("HLT_IsoMu27_v"),
)


