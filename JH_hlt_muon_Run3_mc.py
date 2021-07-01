# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: hlt_muon --python_filename=JH_hlt_muon_Run3_mc.py --step RAW2DIGI,HLT:MuonHLT --process MYHLT --era=Run3 --mc --conditions=113X_mcRun3_2021_realistic_v10 --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForDoubletRemoval --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForCscSegment --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeMuonHLTForGEM --customise=HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizerFuncForMuonHLTSeeding --filein=root://xrootd-cms.infn.it//store/mc/Run3Winter20DRPremixMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8_HCAL/GEN-SIM-DIGI-RAW/110X_mcRun3_2021_realistic_v6-v2/280000/01F2B624-56C1-8940-A735-14F363072EBA.root -n 500 --no_output --no_exec
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
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff") #JH
process.load('HLTrigger.Configuration.HLT_MuonHLT_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/mc/Run3Winter20DRPremixMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8_HCAL/GEN-SIM-DIGI-RAW/110X_mcRun3_2021_realistic_v6-v2/280000/01F2B624-56C1-8940-A735-14F363072EBA.root'),
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
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
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
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('hlt_muon nevts:500'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '113X_mcRun3_2021_realistic_v10', '')

# -- Ntupler -- #
doNtuple = True
if doNtuple:
    skimDY = True  # set True (False) for DY (otherwise)
    isDIGI = True  # set True (False) for DIGI (otherwise)
    from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
    process = customizerFuncForMuonHLTNtupler(process, "MYHLT", skimDY, isDIGI, "Run3v6")
    from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
    process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", skimDY, isDIGI, "Run3v6")
    process.HLTAnalyzerEndpath = cms.EndPath( process.hltPreHLTAnalyzerEndpath + process.hltL1TGlobalSummary + process.hltTrigReport )
    process.HLTSchedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_IsoMu24_v13, process.HLT_Mu50_v13,  process.HLTriggerFinalPath, process.mypath, process.myseedpath, process.HLTAnalyzerEndpath])
    #process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step)
process.schedule = cms.Schedule() #JH
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3 import customizeMuonHLTForDoubletRemoval,customizeMuonHLTForCscSegment,customizeMuonHLTForGEM,customizerFuncForMuonHLTSeeding 

#call to customisation function customizeMuonHLTForDoubletRemoval imported from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
process = customizeMuonHLTForDoubletRemoval(process)

#call to customisation function customizeMuonHLTForCscSegment imported from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
process = customizeMuonHLTForCscSegment(process)

#call to customisation function customizeMuonHLTForGEM imported from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
process = customizeMuonHLTForGEM(process)

#call to customisation function customizerFuncForMuonHLTSeeding imported from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
#process = customizerFuncForMuonHLTSeeding(process)
process = customizerFuncForMuonHLTSeeding(process, newProcessName='MYHLT', doSort=False, mvaCutBs = (0., 0.), mvaCutEs = (0., 0.)) #JH

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
