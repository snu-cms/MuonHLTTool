import FWCore.ParameterSet.Config as cms

DYmuSkimmer = cms.EDFilter("DYmuSkimmer",
  pu_info_src      = cms.InputTag("addPileupInfo"),
  gen_particle_src = cms.InputTag("genParticles"),
  pu_min           = cms.double(-999.),
  pu_max           = cms.double(999.)
)
