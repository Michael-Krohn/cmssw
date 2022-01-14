# Information

Performs skim on RAW dataset. 

First by RECO'ing dataset and then applying Hard Brem control region cuts. Final output is the AOD event content with the HCAL RECO hits added.

## Config file: 

python/ConfFile_MuPlusX_RAW_RECO_cfg.py - Runs the RECO setup, the skim plugin (plugins/MuPlusXFilter_HardBrem_RAW_RECO.cc), and modifies the output to contain AOD + HCAL RECO hits.

