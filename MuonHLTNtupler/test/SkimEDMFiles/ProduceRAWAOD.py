# -- input: RAW & AOD
# -- ouptut: RAWAOD

# -- it can be used for the efficinecy study which require frequent re-running HLT
# -- as CRAB fails a lot if it require "useParant" option (file open error)

import FWCore.ParameterSet.Config as cms

process = cms.Process( "User" )

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/k/kplee/work/private/ROOTFile_Test/AOD_HLTPhysics_Run2017Cv3_Run301567.root',
    ),
    secondaryFileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/k/kplee/work/private/ROOTFile_Test/RAW_HLTPhysics_Run2017Cv1_Run301567_Parent1.root',
        'file:/afs/cern.ch/user/k/kplee/work/private/ROOTFile_Test/RAW_HLTPhysics_Run2017Cv1_Run301567_Parent2.root',
        'file:/afs/cern.ch/user/k/kplee/work/private/ROOTFile_Test/RAW_HLTPhysics_Run2017Cv1_Run301567_Parent3.root'
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.outputRAWAOD = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "RAWAOD.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string( 'RAWAOD' ),
        filterName = cms.untracked.string( '' )
    ),    
    # -- UPDATE: drop irrelvant branches
    outputCommands = cms.untracked.vstring( 'keep *' )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

process.outputPath = cms.EndPath( process.outputRAWAOD )