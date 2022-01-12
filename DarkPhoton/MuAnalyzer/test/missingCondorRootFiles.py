import os


#rootFiles = 'DYJetsToLL_rootFiles.txt'
rootFiles = 'SingleMuon_Run2018D_rootFiles.txt'

with open(rootFiles) as f:
    lines = f.read().splitlines()

jobs=[]
for line in lines:
#    print line.split('/')[6][41:-5]
    jobs.append(line.split('/')[6][41:-5] + "_")



condorSubmitFile = 'condor_submit.txt'
with open(condorSubmitFile) as f:
    lines2 = f.read().splitlines()

newCondorSubmitFileTxt=''
lineNum=0
saveNextLine=False
keepLine=True
for line in lines2:
    if lineNum < 8: newCondorSubmitFileTxt+=line+'\n'
    elif 'Arguments' in line:
	#print line.split()[5]
	for job in jobs:
	    if job in line.split()[5]:
	#	print "REMOVE THIS JOB: ", job
		keepLine=False
		continue
	if keepLine:
	    print "KEEP THIS JOB: ", line.split()[5]
	    newCondorSubmitFileTxt+=line+'\n'
	    saveNextLine=True
	keepLine=True
    elif len(line.split()) == 1 and saveNextLine:
	saveNextLine=False
	newCondorSubmitFileTxt+=line+'\n'
    else: continue
	
    lineNum=lineNum+1

outputFile = open('condor_submit_SingleMuon_Run2018D.txt','w')
outputFile.write(newCondorSubmitFileTxt)
outputFile.close()


