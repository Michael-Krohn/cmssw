## To run the code locally
```
cmsRun python/ConfFile_FilteredAOD_cfg.py runLocally=True
```

## To run on condor
```
perl condor_filelist.perl python/ConfFile_FilteredAOD_cfg.py Filtered_Files_DY_2017.txt isMC=True runRandomTrack=False runLocally=False --prodSpace /data/cmszfs1/user/revering --batch 1 --jobname DY_2017_HighPtMuon
```


## Plotting script
Makes 2D efficiency plots by dividing histograms

```
python test/plotTrackEfficiency.py  -d DY_2017/ -i DY_2017.root --numerator demo/TrackerTrackMatched_EtaPhi_posEta_dR0p35 --denominator demo/TrackerTrack_EtaPhi_posEta
```
