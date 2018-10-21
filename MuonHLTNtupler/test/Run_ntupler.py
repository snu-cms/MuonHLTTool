import FWCore.ParameterSet.Config as cms

process = cms.Process("ntupler")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2018A/SingleMuon/AOD/PromptReco-v1/000/316/187/00000/1CCE3B04-E457-E811-A80C-FA163E0178DF.root'),
    secondaryFileNames = cms.untracked.vstring(),
    # lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v9'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

# -- ntupler -- #
flag_HLTRerun = False

from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase

process.ntupler = ntuplerBase.clone()
process.ntupler.OfflineMuon = cms.untracked.InputTag("muons")
process.ntupler.L3Muon = cms.untracked.InputTag("hltIterL3MuonCandidates")
process.ntupler.L2Muon = cms.untracked.InputTag("hltL2MuonCandidates")
process.ntupler.TriggerResults = cms.untracked.InputTag("TriggerResults", "", "HLT")
process.ntupler.TriggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT")
process.ntupler.OfflineLumiScaler = cms.untracked.InputTag("scalersRawToDigi")

if flag_HLTRerun: # -- after HLT re-run -- #
  process.ntupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis", "Muon", newProcessName)
  process.ntupler.MyTriggerResults = cms.untracked.InputTag("TriggerResults", "", newProcessName) # -- result after rerun HLT -- #
  process.ntupler.MyTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD", "", newProcessName) # -- result after rerun HLT -- #
  process.ntupler.LumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi", "", newProcessName)
  process.ntupler.IterL3MuonNoID   = cms.untracked.InputTag("hltIterL3MuonsNoID", "", newProcessName)
else: # -- without HLT re-run -- #
  process.ntupler.L1Muon = cms.untracked.InputTag("gmtStage2Digis", "Muon", "RECO")
  process.ntupler.MyTriggerResults = cms.untracked.InputTag("TriggerResults")
  process.ntupler.MyTriggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD")
  process.ntupler.LumiScaler = cms.untracked.InputTag("scalersRawToDigi")

process.mypath = cms.EndPath(process.ntupler)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("ntuple.root"),
	closeFileFast = cms.untracked.bool(False),
	)

process.MessageLogger = cms.Service( "MessageLogger",
	destinations = cms.untracked.vstring("cerr"),
	cerr = cms.untracked.PSet( threshold = cms.untracked.string('ERROR'), ),
	)