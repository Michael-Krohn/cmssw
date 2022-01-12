import os
import ROOT


logFiles = 'SingleMuon_Run2018D_rootFiles.txt'

with open(logFiles) as f:
    lines = f.read().splitlines()

jobs=[]
for line in lines:
    line = line.strip()
    SignalFile = ROOT.TFile.Open(str(line))
    if len(SignalFile.GetListOfKeys()) == 0:
	jobs.append(line)

outputTxt=''
for job in jobs:
    outputTxt+=job+'\n'

'''
condorSubmitFile = 'condor_submit.txt'
with open(condorSubmitFile) as f:
    lines2 = f.read().splitlines()

newCondorSubmitFileTxt=''
lineNum=0
saveNextLine=False
for line in lines2:
    if lineNum < 8: newCondorSubmitFileTxt+=line+'\n'
    elif 'Arguments' in line:
#	print line.split()
	for job in jobs:
	    if job in line.split()[5]:
		print "KEEP THIS JOB: ", job
		newCondorSubmitFileTxt+=line+'\n'
		saveNextLine=True
    elif len(line.split()) == 1 and saveNextLine:
	saveNextLine=False
	newCondorSubmitFileTxt+=line+'\n'
    else: continue
	
    lineNum=lineNum+1
'''

outputFile = open('BadRootFiles_Run2018D.txt','w')
outputFile.write(outputTxt)
outputFile.close()



