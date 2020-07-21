# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase2_realistic_T15 -n 10 --era Phase2C9 --eventcontent RECOSIM,DQM --runUnscheduled -s RAW2DIGI,RECO:reconstruction_trackingOnly,VALIDATION:@trackingOnlyValidation,DQM:@trackingOnlyDQM --datatier GEN-SIM-RECO,DQMIO --geometry Extended2026D49 --filein file:step2.root --fileout file:step3.root


import FWCore.ParameterSet.Config as cms
# from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
# process = cms.Process('MYHLT',Phase2C9)

from Configuration.StandardSequences.Eras import eras
process = cms.Process("MYHLT", eras.Phase2C9)
# process = cms.Process("MYHLT", eras.Phase2C9_trigger)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('Configuration.StandardSequences.Validation_cff')
# process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
# process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/relval/CMSSW_11_0_0_pre11/RelValZMM_14/GEN-SIM-DIGI-RAW/110X_mcRun4_realistic_v2_2026D41noPU-v1/10000/EE56DFF6-6AEB-BD4C-97C0-CE2D32452D25.root'),
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/relval/CMSSW_11_1_0_pre4_GEANT4/RelValZMM_14/GEN-SIM-DIGI-RAW/110X_mcRun4_realistic_v3_2026D49noPU-v1/10000/FEB6F9C7-9090-5949-9877-4FA2EFD8CDD8.root'),
        #root://xrootd-cms.infn.it//store/mc/PhaseIITDRSpring19DR/Mu_FlatPt2to100-pythia8-gun/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v2/70000/FFCFF986-ED0B-B74F-B253-C511D19B8249.root'),
       # 'file:step2.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string("DQMIO.root"),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# Path and EndPath definitions
# process.raw2digi_step = cms.Path(process.RawToDigi)
# process.reconstruction_step = cms.Path(process.reconstruction_trackingOnly)
# process.prevalidation_step = cms.Path(process.globalPrevalidationTrackingOnly)
# process.validation_step = cms.EndPath(process.globalValidationTrackingOnly)
# process.dqmoffline_step = cms.EndPath(process.DQMOfflineTracking)
# process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
# process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
# process.DQMoutput_step = cms.EndPath(process.DQMoutput)


# -- L1TT -- #
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")

from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)

# L1TRKALGO == 'HYBRID'
process.load("L1Trigger.TrackFindingTracklet.Tracklet_cfi")
L1TRK_PROC  =  process.TTTracksFromTrackletEmulation
L1TRK_NAME  = "TTTracksFromTrackletEmulation"
L1TRK_LABEL = "Level1TTTracks"

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")  # defined above
process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")
process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag(L1TRK_NAME, L1TRK_LABEL) )

process.TTTracksEmulation = cms.Path(process.offlineBeamSpot*L1TRK_PROC)
process.TTTracksEmulationWithTruth = cms.Path(process.offlineBeamSpot*process.TTTracksFromTrackletEmulation*process.TrackTriggerAssociatorTracks)

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.L1simulation_step = cms.Path(process.SimL1Emulator)
# -- #


# -- Setup -- #
process.source.fileNames = cms.untracked.vstring(
    "root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRWinter20DIGI/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/PU200_pilot2_110X_mcRun4_realistic_v3-v2/270000/FB5CCEB0-775F-9D4D-B422-15DDE09B88A9.root"
    # "file:/eos/cms/store/mc/Phase2HLTTDRWinter20DIGI/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/110000/0A86D9A3-925B-1A47-963C-097E662902C1.root"
    # "file:/eos/user/m/moh/TestSamples/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/4043F7D2-1BF4-FE40-82D9-10786D005454.root"
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/B7D073FD-FB58-1647-A60C-9D14FAF5AD0E.root",
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/FF34EE7B-54D1-8547-A44C-56D95CEEAE23.root",
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/FF109638-EEC3-924F-B812-468287EA690C.root",
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/FEF22648-0D68-3B4E-BDBE-0912CFD05FAA.root",
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/FEE1D95E-8B38-F141-AE4D-B711B970D59F.root",
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/FE04F8CA-FB66-3949-8C04-36841AA00400.root",
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/FD93A4A3-A7D2-8B46-BEF2-0BFD48301349.root",
    # "file:/data/user/moh/DYToLL_M-50_TuneCP5_14TeV-pythia8__Phase2HLTTDRWinter20DIGI-PU200_pilot_110X_mcRun4_realistic_v3-v2__GEN-SIM-DIGI-RAW/FD72C652-FEC4-1240-8AB5-45975D2D17AA.root"
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( -1 )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
    numberOfThreads = cms.untracked.uint32( 4 ),
    numberOfStreams = cms.untracked.uint32( 0 ),
    sizeOfStackForThreadsInKB = cms.untracked.uint32( 10*1024 )
)

if 'MessageLogger' in process.__dict__:
    process.MessageLogger.categories.append('TriggerSummaryProducerAOD')
    process.MessageLogger.categories.append('L1GtTrigReport')
    process.MessageLogger.categories.append('L1TGlobalSummary')
    process.MessageLogger.categories.append('HLTrigReport')
    process.MessageLogger.categories.append('FastReport')
# -- #

# -- Schedule -- #
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("outfile.root"),
    closeFileFast = cms.untracked.bool(False)
)

process.GenMuAnalyzerTEST = cms.EDAnalyzer("GenMuAnalyzer",
    genParticle_src = cms.InputTag("genParticles"),
    L1TT_src = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
    L1TkMuon_src = cms.InputTag("L1TkMuons"),
    pt_min = cms.double(0.0)
)

process.mypath = cms.Path( process.GenMuAnalyzerTEST )

process.schedule = cms.Schedule( process.TTTracksEmulationWithTruth, process.L1simulation_step, process.mypath )


from HLTrigger.Configuration.Eras import modifyHLTforEras
modifyHLTforEras(process)


# print process.dumpPython()

# process.Timing = cms.Service("Timing",
#     summaryOnly = cms.untracked.bool(True),
#     useJobReport = cms.untracked.bool(True)
# )

# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#     ignoreTotal = cms.untracked.int32(1)
# )

