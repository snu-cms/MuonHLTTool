import sys, os

names = [
    'HLTPhysics2018C_Run319941_RAWAOD',
    'HLTPhysics2018E_Run325308_RAWAOD',
    'SingleMuon2018C_Run319941_NMu1_Pt27to1000000000_PU40to60_RAWAOD',
    'SingleMuon2018C_Run319941_NMu2_Pt27to1000000000_PU40to60_RAWAOD'
]

for name in names:
    f = open('HLTCfgData_ntuple_IsoMu24.py', 'r')
    fData = f.read()
    f.close()
    newData = fData.replace('NAMEREPLACE',name)

    fnew = 'HLTCfgData_ntuple_IsoMu24_'+name+'.py'
    fn = open(fnew,'w')
    fn.write(newData)
    fn.close()

    cmd = 'cmsRun  %s  >&log_%s.log&' % (fnew, name)
    print cmd
    os.system(cmd)


# xrdcp -r /u/user/msoh/MuonHLT/HLTTDR/CMSSW_10_2_6/src/MuonHLTTool/MuonHLTNtupler/test/Local/test_skim_20190820/ntuple/Outputs root://cluster142.knu.ac.kr//store/user/moh/RAWAODSkimV00/v20190820/Ntuples


