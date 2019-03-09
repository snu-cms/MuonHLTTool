import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonHLT")

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/k/kplee/work/private/MuonHLT/v20180608_v01_BTagMuIssue_muonInJet/v04_PrintOut/CMSSW_10_1_7/src/WorkingDir/output_menu2p1_BTagMuFailed.root'),
	secondaryFileNames = cms.untracked.vstring(),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v9'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

from MuonHLTTool.EDAnalyzers.MHPrintBTagMuTriggerInfo_cfi import *

process.myMHPrintBTagMuTriggerInfo = MHPrintBTagMuTriggerInfo.clone()
# process.myMHPrintBTagMuTriggerInfo.triggerResults = cms.untracked.InputTag("TriggerResults::HLT")
# process.myMHPrintBTagMuTriggerInfo.triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT")

process.mypath = cms.Path(process.myMHPrintBTagMuTriggerInfo)