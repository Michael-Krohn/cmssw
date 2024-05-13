# Injecting leptons

## Brief overview

Creating injecting leptons occurs in 5 steps:
1. Select semi-leptonic ttbar events (in data or MC). And save the lead AK8 jet eta,phi directions and the PV coordinates
2. Create a skim of the RAW events selected
3. For each event selected, create an LHE event of the lepton to be injected using the jet direction and PV coordinates from 1.
4. Process the LHE file and RAW skim through RECO
5. Merge the RECO collections of the LHE lepton and the original RAW event and run them through MINIAOD

## Detailed instructions of each step

1. This step is done outside of this code framework. It can be done in any analysis structure. An output txt file with the following information per line should be created:
   runNum lumiSec evtNum jet0_pt jet0_eta jet0_phi pvx pvy pvz

2. Skim the RAW events that were selected and also run the RECO over them. Create the config with the proper cmsDriver command for the era of data/MC you are running over and add:
   --customise TauAnalysis/MCEmbeddingTools/customisers.customiseSelecting_Reselect

   An example config that was used on 2016 data can be found here: TauAnalysis/MCEmbeddingTools/python/selection.py
   
3.	The RECO step has to be run again, but this time add the following option to the cmsDriver command:
   
   --customise TauAnalysis/MCEmbeddingTools/customisers.customiseLHEandCleaning_Reselect

   This uses the plugin EmbeddingLHEProducer.cc, which is very hacky. It reads in the local file created in step 1, but the file is hardcoded here (feel free to change this!): https://github.com/Michael-Krohn/cmssw/blob/lepton-injection/TauAnalysis/MCEmbeddingTools/plugins/EmbeddingLHEProducer.cc#L184

   An example config that was used on 2016 data can be found here: TauAnalysis/MCEmbeddingTools/python/lheprodandcleaning.py

4. To process the LHE file, use the same cmsDriver RECO command but add:
   
   cmsDriver TauAnalysis/MCEmbeddingTools/python/EmbeddingPythia8Hadronizer_cfi.py \
   ...\
   --customise TauAnalysis/MCEmbeddingTools/customisers.customiseGenerator_Reselect

   An example config that was used on 2016 data can be found here: TauAnalysis/MCEmbeddingTools/python/generator.py

5. Use the standard cmsDriver MiniAOD comand for the era and MC/data that you are processing, but add the following option:

   --customise TauAnalysis/MCEmbeddingTools/customisers.customiseMerging_Reselect

   An example config that was used on 2016 data can be found here: TauAnalysis/MCEmbeddingTools/python/merging.py


