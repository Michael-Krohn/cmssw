# Injecting leptons

## Brief overview

Creating injecting leptons occurs in 5 steps:
1. Select semi-leptonic ttbar events (in data or MC). And save the lead AK8 jet eta,phi directions and the PV coordinates
2. Create a skim of the RAW events selected
3. For each event selected, create an LHE event of the lepton to be injected using the jet direction and PV coordinates from 1.
4. Process the LHE file and RAW skim through RECO
5. Merge the collections of the LHE lepton and the original RAW event

## Detailed instructions

1. This step was done outside of this code framework. It can be done in any analysis structure. An output txt file with the following information per line should be created:
   runNum lumiSec evtNum jet0_pt jet0_eta jet0_phi pvx pvy pvz
2.
3. Run over ttbar/data events and select events with a semi-leptonic decay. You want to create a list of the events (runNum, lumiSec, evtNum) along with the eta, phi direction of the lead AK8 jet (hadronically decaying jet) and the PV coordinates2.	You need to create a skim of the RAW events that were selected.
4.	Create LHE events for the individual lepton, use the jet direction and the PV information (txt file hardcoded in EmbeddingLHEProducer.cc)
5.	Run the LHE files through RECO
6.	Merge the individual lepton events with the original ttbar/data events

