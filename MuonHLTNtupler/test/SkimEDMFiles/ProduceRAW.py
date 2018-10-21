# -- input; RAW
# -- output: RAW

# -- it can be used for frequent rate study
# -- submit CRAB job with a JSON & specific run range

import FWCore.ParameterSet.Config as cms

process = cms.Process( "User" )

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        'root://eoscms.cern.ch//eos/cms/store/data/Run2017F/EphemeralHLTPhysics1/RAW/v1/000/305/636/00000/3CC1B3D9-95B9-E711-89B0-02163E01A2D2.root',
    ),

    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# -- output
process.output = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "RAW.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string( 'RAW' ),
        filterName = cms.untracked.string( '' )
    ),
    outputCommands = cms.untracked.vstring( 'keep *' )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

process.outputPath = cms.EndPath( process.output )