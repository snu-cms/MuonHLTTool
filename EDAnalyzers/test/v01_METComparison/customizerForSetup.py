# add specific customizations
import FWCore.ParameterSet.Config as cms

def customizeForSetup( process ):
	inputFileNames = [
	'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/50000/0043237F-0DA1-E711-AD49-346AC29F11B8.root',
	'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/50000/0078B52F-F5A0-E711-BCF6-0090FAA57360.root',
	'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/50000/008E6A3C-BEA0-E711-BB0F-FA163E86F14F.root',
	'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/50000/008F1C4B-10A1-E711-8E3C-008CFAC942DC.root',
	'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/50000/00A1E603-FEA0-E711-A3B4-FA163E12D9E7.root',
	]

	_customInfo2 = {}
	_customInfo2['menuType'  ]= "GRun"
	_customInfo2['globalTags']= {}
	_customInfo2['globalTags'][True ] = "auto:run2_hlt_GRun"
	_customInfo2['globalTags'][False] = "auto:run2_mc_GRun"
	_customInfo2['inputFiles']={}
	_customInfo2['inputFiles'][True]  = "file:RelVal_Raw_GRun_DATA.root"
	_customInfo2['inputFiles'][False] = "file:RelVal_Raw_GRun_MC.root"
	_customInfo2['maxEvents' ]= -1
	_customInfo2['globalTag' ]= "101X_mc2017_realistic_TSG_2018_04_09_20_43_53"
	_customInfo2['inputFile' ]= inputFileNames
	_customInfo2['realData'  ]= False

	from HLTrigger.Configuration.customizeHLTforALL import customizeHLTforAll
	process = customizeHLTforAll(process,"GRun",_customInfo2)