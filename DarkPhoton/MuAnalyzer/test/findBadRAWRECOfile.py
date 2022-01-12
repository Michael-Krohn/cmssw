import os


logFiles = 'BadRootFiles_Run2018D.txt'

with open(logFiles) as f:
    lines = f.read().splitlines()

#/local/cms/user/krohn045/ConfFile_MuPXAnalyzer_HardBrem_Run2018D_RECO_Analyzer/ConfFile_MuPXAnalyzer_HardBrem_Run2018D_RECO_Analyzer_502.root
#/local/cms/user/krohn045/ConfFile_MuPXAnalyzer_HardBrem_Run2018D_RECO_Analyzer/

logs=''
for line in lines:
#    print line[:-62]
#    print line[79:-5]
    os.system('cp '+line[:-62]+'log/'+line[79:-5]+'.err.gz'+' logFiles_Run2018D_RAWRECO/')
    os.system('gzip -d logFiles_Run2018D_RAWRECO/'+line[79:-5]+'.err.gz')
    logs+='logFiles_Run2018D_RAWRECO/'+line[79:-5]+'.err'+'\n'


outputFile = open('Run2018D_LogFiles.txt','w')
outputFile.write(logs)
outputFile.close()

logFilesUntar = 'Run2018D_LogFiles.txt'

with open(logFilesUntar) as f:
    lines2 = f.read().splitlines()

badFilesTxt=''
for logFile in lines2:
#    print logFile
    with open(logFile) as f:
	lines3 = f.read().splitlines()
    for line in lines3:
        if 'Input' in line: 
#	    print line
#	    print "line.split(): ", line.split()
	    badFilesTxt+=line.split()[3]+'\n'
	elif 'Closed' in line:
#            print line
#            print "line.split(): ", line.split()
	    badFilesTxt+=line.split()[5]+'\n'
        else: continue

outputFile = open('Run2018D_BadRAWRECOfiles.txt','w')
outputFile.write(badFilesTxt)
outputFile.close()



