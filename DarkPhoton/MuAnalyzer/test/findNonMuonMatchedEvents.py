import os

'''
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
'''
logFilesUntar = 'DYJetsToLL_LogFiles.txt'

with open(logFilesUntar) as f:
    lines2 = f.read().splitlines()

pdgidNumbersTxt=''
numberOfMuons=0.
for logFile in lines2:
    with open(logFile) as f:
	lines3 = f.read().splitlines()
    for line in lines3:
        if 'Track matched' in line: 
	    #print line
	    #print "line.split(): ", line.split()
	    pdgid = line.split()[3]
	elif 'Well isolated track not matched to loose Muon' in line:
	    print line
	    print pdgid
	    pdgidNumbersTxt+=pdgid+"\n"
	    if abs(float(pdgid))==13:
		numberOfMuons = numberOfMuons + 1
        else: continue

print "numberOfMuons: ", numberOfMuons
outputFile = open('DYJetsToLL_PDGIDNumbers.txt','w')
outputFile.write(pdgidNumbersTxt)
outputFile.close()



