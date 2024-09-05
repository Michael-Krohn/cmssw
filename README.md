
To install code:

```
scram project CMSSW_10_6_17
cd CMSSW_10_6_17/src
cmsenv
git cms-init
git cms-merge-topic revering:DarkPhoton_10_6_17
git checkout revering/DarkPhoton_10_6_17
scram b
```

Git adds a bunch of unnecessary packages to the sparse checkout for some reason. If they show up when you run cms-merge-topic (Validation/ SimRomanPlot/ and L1Trigger/ like to show up for some reason), remove go into src/.git/info/sparse-checkout and delete every folder except the /DarkPhoton/ and /SimG4Core/CustomPhysics/ ones. Then do
```
git read-tree -mu HEAD
```
to get rid of the extra folders and reduce the compile time.
