import FWCore.ParameterSet.Config as cms

MHMETAnalyzer = cms.EDAnalyzer('MHMETAnalyzer',
    onlinePFMET = cms.untracked.InputTag("hltPFMETProducer::HLT"),
    genMET = cms.untracked.InputTag("genMetTrue::HLT"),
)


