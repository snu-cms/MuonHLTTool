#!/usr/bin/env python
import subprocess
import os

crabDirBase = "CRABDir"
proxy = '"/tmp/x509up_u72027"'

FileList = os.listdir("./%s" % crabDirBase)
List_CRABDir = []
print "available crabDir list: "
for filename in FileList:
	if "crab_" in filename:
		print "'"+filename+"',"
		List_CRABDir.append( filename )

# CRABDirs = [
# # 'crab_DYntuple_v20170428_80XReMiniAOD_FixEvtNum_SingleMuon_Run2016G_v1_GoldenJSON_271036_to_284044',
# ]

# CRABDirs.sort()

# print "Selected file list: "
# print CRABDirs

CRABDirs = List_CRABDir # -- all crab directories in crabDirBase -- #
CRABDirs.sort()

# for File in FileList:
# 	if "crab_" in File and os.path.isdir( File ):
# 		print "Recognized crab directory: " + File
# 		CRABDirs.append( File )

ResubmtCMD = []
CompletedList = []
UnknownList = []
OthersList = []

for crabDir in CRABDirs:
	crabDirPath = "%s/%s" % (crabDirBase, crabDir)
	# outputDir = "v" + crabDir.split("_v")[1]
	
	cmd = 'crab status "'+crabDirPath+'" --proxy='+proxy
	result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(stdout, stderr) = result.communicate()
	print "#" * 100
	print cmd+'\n'
	print "[stdout]"
	print stdout
	print "[stderr]"
	print stderr
	print "#" * 100 +'\n\n'

	if "failed" in stdout:
		ResubmtCMD += ['crab resubmit '+crabDirPath+' --proxy='+proxy]

	elif "COMPLETED" in stdout:
		CompletedList.append( crabDir )

	elif "UNKNOWN" in stdout:
		UnknownList.append( crabDir )

	else:
		OthersList.append( crabDir )

print "[Completed list]"
for one in CompletedList:
	print "'"+one+"',"

print "[Unknown list]"
for one in UnknownList:
	print "'"+one+"',"

print "[Others]"
for one in OthersList:
	print "'"+one+"',"

print "\n[CRAB jobs which should be resubmitted]"
for CMD in ResubmtCMD:
	print CMD

for CMD in ResubmtCMD:
	# 	os.system( CMD )
	result = subprocess.Popen(CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(stdout, stderr) = result.communicate()
	print "#" * 100
	print CMD+'\n'
	print "[stdout]"
	print stdout
	print "[stderr]"
	print stderr
	print "#" * 100 +'\n\n'
