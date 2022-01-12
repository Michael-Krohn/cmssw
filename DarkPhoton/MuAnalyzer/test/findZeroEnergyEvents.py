import os


logFiles = 'DYJetsToLL_logFiles.txt'

with open(logFiles) as f:
    lines = f.read().splitlines()

logs=''
for line in lines:
    print line.split('/')[7]
    os.system('cp '+line+' logFiles_DYJetsToLL/')
    os.system('gzip -d logFiles_DYJetsToLL/'+line.split('/')[7])
    logs+='logFiles_DYJetsToLL/'+line.split('/')[7][:-3]+'\n'


outputFile = open('DYJetsToLL_LogFiles.txt','w')
outputFile.write(logs)
outputFile.close()

logFilesUntar = 'DYJetsToLL_LogFiles.txt'

with open(logFilesUntar) as f:
    lines2 = f.read().splitlines()

eventNumbersTxt=''
for logFile in lines2:
    with open(logFile) as f:
	lines3 = f.read().splitlines()
    for line in lines3:
	if 'No cone energy!' in line:
	    print "THIS FILE"
	    print "logFile: ", logFile
#        if 'luminosityBlock' in line: 
#	    print line
#	    print "line.split(): ", line.split()
#	    lumi = line.split()[1]
#        elif 'event()' in line:
#	    eventNum = line.split()[1]
#        elif 'run()' in line:
#	    runNum = line.split()[1]
#	    eventNumbersTxt+=runNum+":"+lumi+":"+eventNum+"\n"
        else: continue

#outputFile = open('Run2018C_EventNumbers.txt','w')
#outputFile.write(eventNumbersTxt)
#outputFile.close()



