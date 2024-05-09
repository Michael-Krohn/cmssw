#Injecting leptons

##Brief overview

Creating injecting leptons occurs in 5 steps:
1. Select semi-leptonic ttbar events (in data or MC). And save the lead AK8 jet eta,phi directions and the PV coordinates
2. Create a skim of the RAW events selected
3. For each event selected, create an LHE event of the lepton to be injected using the jet direction and PV coordinates from 1.
4. Process the LHE file and RAW skim through RECO
5. Merge the collections of the LHE lepton and the original RAW event

##Details

First step of the lepton injection workflow are done outside of this code framework. We had some run over ntuples and print out the information:
1.	Run over ttbar/data events and select events with a semi-leptonic decay. You want to create a list of the events (runNum, lumiSec, evtNum) along with the eta, phi direction of the lead AK8 jet (hadronically decaying jet) and the PV coordinates
2.	You need to create a skim of the RAW events that were selected.
3.	Create LHE events for the individual lepton, use the jet direction and the PV information (txt file hardcoded in EmbeddingLHEProducer.cc)
4.	Run the LHE files through RECO
5.	Merge the individual lepton events with the original ttbar/data events

