import os
import ROOT

logFiles = 'SingleMuon_Run2018D_rootFiles.txt'
#logFiles = 'DYJetsToLL_rootFiles.txt'

with open(logFiles) as f:
    lines = f.read().splitlines()

jobs=[]
for line in lines:
    line = line.strip()
    SignalFile = ROOT.TFile.Open(str(line))
    if len(SignalFile.GetListOfKeys()) == 0:
#    if SignalFile == None:
        print line
        os.system('rm '+line)




