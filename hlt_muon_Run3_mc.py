# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: hlt_muon --python_filename=hlt_muon_Run3_mc.py --step HLT:MuonHLT --process MYHLT --era=Run3 --mc --conditions=auto:phase1_2021_realistic --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeDoubleMuIsoFix --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForDoubletRemoval --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForCscSegment --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForGEM --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForPatatrackWithIsoAndTriplets --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeIOSeedingPatatrack --filein=/store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120003/e786c41e-21ba-489f-880c-42d0a248e59e.root -n 300 --no_output --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('MYHLT',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HLTrigger.Configuration.HLT_MuonHLT_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(300),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120003/e786c41e-21ba-489f-880c-42d0a248e59e.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('hlt_muon nevts:300'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# -- L2 seed stat recovery -- #
#process.hltIterL3MuonPixelTracksTrackingRegions.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.hltIter3IterL3MuonL2Candidates.src = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.hltL3MuonsIterL3IO.L3TrajBuilderParameters.MuonTrackingRegionBuilder.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence = cms.Sequence(
#    process.HLTIterL3OImuonTkCandidateSequence +
#    process.hltIterL3OIL3MuonsLinksCombination +
#    process.hltIterL3OIL3Muons +
#    process.hltIterL3OIL3MuonCandidates +
#    #process.hltL2SelectorForL3IO +
#    process.HLTIterL3IOmuonTkCandidateSequence +
#    process.hltIterL3MuonsFromL2LinksCombination
#)

# -- Ignoring filters -- #
process.HLT_Mu50_v13 = cms.Path(
    process.HLTBeginSequence +
    cms.ignore(process.hltL1sSingleMu22or25) +
    cms.ignore(process.hltPreMu50) +
    cms.ignore(process.hltL1fL1sMu22or25L1Filtered0) +
    process.HLTL2muonrecoSequence +
    cms.ignore(process.hltL2fL1sMu22or25L1f0L2Filtered10Q) +
    process.HLTL3muonrecoSequence +
    cms.ignore(process.hltL1fForIterL3L1fL1sMu22or25L1Filtered0) +
    cms.ignore(process.hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q) +
    process.HLTEndSequence
)

skimDY = True         # set True (False) for DY (otherwise)
isDIGI = True         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW)
MvaVersion = "v8Pre"  # set v7Fast or v8Pre
WP = 0.0
doOI = False

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", skimDY, isDIGI, MvaVersion)
process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")
from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", skimDY, isDIGI, MvaVersion)
process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered",       "", "MYHLT")
process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered", "", "MYHLT")

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v13,
     process.HLT_Mu50_v13,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myseedpath
)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3 import customizeDoubleMuIsoFix,customizeMuonHLTForDoubletRemoval,customizeMuonHLTForCscSegment,customizeMuonHLTForGEM,customizeMuonHLTForPatatrackWithIsoAndTriplets,customizeIOSeedingPatatrack
process = customizeDoubleMuIsoFix(process)
process = customizeMuonHLTForDoubletRemoval(process)
process = customizeMuonHLTForCscSegment(process)
process = customizeMuonHLTForGEM(process)
process = customizeMuonHLTForPatatrackWithIsoAndTriplets(process)
process = customizeIOSeedingPatatrack(process, mvaCutBs = (WP, WP), mvaCutEs = (WP, WP))

if doOI == True :
  from RecoMuon.TrackerSeedGenerator.customizeOIseeding import customizeOIseeding
  process = customizeOIseeding(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

if 'MessageLogger' in process.__dict__:
    process.MessageLogger.TriggerSummaryProducerAOD = cms.untracked.PSet()
    process.MessageLogger.L1GtTrigReport = cms.untracked.PSet()
    process.MessageLogger.L1TGlobalSummary = cms.untracked.PSet()
    process.MessageLogger.HLTrigReport = cms.untracked.PSet()
    process.MessageLogger.FastReport = cms.untracked.PSet()
    process.MessageLogger.ThroughputService = cms.untracked.PSet()
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
