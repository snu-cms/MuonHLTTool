import FWCore.ParameterSet.Config as cms

WmuSkimmer = cms.EDFilter("WmuSkimmer",
  pu_info_src      = cms.InputTag("addPileupInfo"),
  gen_particle_src = cms.InputTag("genParticles"),
  pu_min           = cms.double(-999.),
  pu_max           = cms.double(999.)
)
