import os


logFiles = 'DYJetsToLL_missingLogFiles.txt'

with open(logFiles) as f:
    lines = f.read().splitlines()

jobs=[]
for line in lines:
    print line.split('/')[7][41:-4]
    jobs.append(line.split('/')[7][41:-4])


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


outputFile = open('condor_submit_DYJetsToLL.txt','w')
outputFile.write(newCondorSubmitFileTxt)
outputFile.close()



