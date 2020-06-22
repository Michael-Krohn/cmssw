## To run the code locally
```
cmsRun python/ConfFile_MCSignal_cfg.py isMC=True runLocally=True isSig=True hasDpho=True
```

## To run on condor
MuAnalyzer
```
perl condor_filelist.perl python/ConfFile_FilteredAOD_cfg.py datafiles/Filtered_Files_DY_2017.txt isMC=True runRandomTrack=False runLocally=False --prodSpace /data/cmszfs1/user/revering --batch 1 --jobname DY_2017_HighPtMuon
```
SigAnalyzer
```
perl condor_filelist.perl python/ConfFile_MCSignal_cfg.py datafiles/map_1p0_unweighted.txt isMC=True runRandomTrack=False runLocally=False isSig=True hasDpho=True --prodSpace /data/cmszfs1/user/myusername --batch 1 --jobname YourJobName
```


