## Current Selection Criteria:

Preselection:

	-Tracks:
		pt>20 GeV
        	highPurity
        	|dxy|<0.05
        	Endcap tracks: 1.3<|eta|<2.4
	-Muons:
        	reco CutBasedIdTight
        	PFIsoTight
        	TkIsoTight        !!Need to remove one of Iso reqs, keep forgetting to do so!!
        	pt>26 GeV
        	|eta|<2.4
        	GlobalMuon

	- Tag/Probe Pairing:
        	First try using highest pt muon as tag, try 2nd highest if first has no valid probes.
		Probe track has opposite charge from tag.
	        Fitted vertex of tag+probe has chi^2/ndof <3.
		Invariant mass of tag+probe pair between 50 and 150 GeV
        	Track must belong to a vertex
        	|iEta| of HE cells along track trajectory <19 and >28
		Nearest HE cells along trajectory not in dead region (iEta<-16, 52<iPhi<64)
        	Track Isolation Requirements:
        		Sum of pt from all tracks within dR of 0.3 originating from the same vertex/pt of probe < 0.2.
			Sum of E of all ECAL hits >0.3 GeV within dR of 0.4 < 30 GeV.
        	No other tracks with dR>0.3 to probe allowed to pair with tag muon (same pairing requirements, other tracks allowed at any eta and otherwise identical reqs).

Current Signal Selections (Still undergoing significant development)

        -Tight Category:
		>2 hits along probe trajectory in HE >0.1 GeV.
                <10 GeV in HE hits along probe trajectory
                No calo jets within dR of 0.5 to probe.
                Min dR to nearest CSC hit >0.05
                Min dR to nearest standalone muon >0.5

	-Medium Category:
		>2 hits along probe trajectory in HE.
                No calo jets withing dR of 0.5 to probe.
		Standalone muon within dR of 0.5
	        Nearby standalone muon has <20% of initial probe momentum.
                <2 GeV in HE hits along probe trajectory.
        

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


