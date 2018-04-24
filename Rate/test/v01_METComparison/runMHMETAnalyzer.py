import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonHLT")

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/k/kplee/work/private/ROOTFile_Test/GENSIMRAW_ttbar_92XForTSG.root'),
	secondaryFileNames = cms.untracked.vstring(),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

from MuonHLTTool.Rate.MHMETAnalyzer_cfi import *

process.myMHMETAnalyzer = MHMETAnalyzer.clone()
# process.myMHMETAnalyzer.onlinePFMET = cms.untracked.InputTag("hltPFMETProducer::HLT")
# process.myMHMETAnalyzer.genMET = cms.untracked.InputTag("genMetTrue::HLT")

process.mypath = cms.Path(process.myMHMETAnalyzer)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("METPlots.root"),
	closeFileFast = cms.untracked.bool(False),
)

