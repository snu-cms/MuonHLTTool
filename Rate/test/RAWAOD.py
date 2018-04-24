import FWCore.ParameterSet.Config as cms

process = cms.Process( "RAWAOD" )

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/k/kplee/work/private/ROOTFile_Test/AODSIM_ttbar_92XForTSG.root',
    ),
    secondaryFileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/k/kplee/work/private/ROOTFile_Test/GENSIMRAW_ttbar_92XForTSG.root',
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
    outputCommands = cms.untracked.vstring( 'keep *', 'drop *_*_*_SIM' )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

process.outputPath = cms.EndPath( process.outputRAWAOD )