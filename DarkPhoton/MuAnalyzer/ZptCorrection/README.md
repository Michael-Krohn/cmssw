# DY Zpt correction

## How to
Below example used ROOT ACLiC. But you can compile it with any method you prefer.
```
root [0] .L ZptCorrection.C+
root [1] ZptCorrection a("ZptWeight_DYJets.root") ### initialize with a proper data file matching with your DY sample
[ZptCorrection::Setup] using file ZptWeight_DYJets.root
(ZptCorrection &) @0x7faae4e63010
root [2] a.GetZptWeight(20.)
(double) 1.0900299
```

There are three correction methods. The first one is the simplest one only considering pt dependence. The others consider additional dependence on rapidity and mass. The first one is recommended as it gives enough correction in general. All input parameters are the value of generator level Z boson before QED FSR.
* `double GetZptWeight(double pt) const`
  * Zpt weight as a function of Z pt
* `double GetZptWeight(double pt,double rapidity) const`
  * Zpt weight as a function of Z pt and rapidity
* `double GetZptWeight(double pt,double rapidity) const`
  * Zpt weight as a function of Z pt, rapidity and mass

## Data files
* `ZptWeight_DYJets.root`: aMC@NLO; DYJetsToLL
* `ZptWeight_DYJets_MG.root`: MadGraph; DYJetsToLL
* `ZptWeight_MiNNLO.root`: POWHEG MiNNLO; DYJetsToEE and DYJetsToMuMu
