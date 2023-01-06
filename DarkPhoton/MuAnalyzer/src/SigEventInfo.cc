#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DarkPhoton/MuAnalyzer/interface/SigEventInfo.h"
#include <iostream>

SigEventInfo::SigEventInfo() {
  eventWeight = 1;
  genEventWeight = 1;
  muRocWeight = 1;
  filledWeights = false;
  dbremZ = -1;
  dbremPerp = -1;
  sigCat = -1;
  bremDepth = -1;
  bremR = 525;
  dbremInitialE = -1;
  dbremFinalE = -1;
  deflectingAngle = 20;
  nTagMuons = -1;
  nPairedProbes = -1;
  probeTrackIso = -1;
  probeEcalIso = -1;
  probeEcalIso0p2 = -1;
  probeEcalIso0p1 = -1;
  probeEcalIso0p05 = -1;
  isLoose = false;
  isTight = false;
  isHighPurity = false;
  minGenMuDr = -1;
  minGenDE = -10;
  probeReducedChi = -1;
  tagMuonEta = 0;
  tagMuonPt = -1;
  minSimCSCDr = 100;
  minSimCSCELoss = -1;
  probeTrackEta = 0;
  probeTrackPhi = 0;
  probeTrackPt = -1.;
  pairVtxChi = -1;
  diMuonMass = -1;
  minCSCImpactParameter = -1;
  cscMinDR = -1;
  cscRecHitMinDR = -1;
  cscChiSq = -1;
  cscHitAngleChange = 100;
  closestApproach = -1;
  standaloneE = -1;
  staMinDr = -1;
  nCellsFound = -1;
  nNeighborCellsFound = -1;
  iEta = 0;
  pfMet_pt=0;
  muMetMT=0;
  muAndProbeMetMT=0;
  cellCenterDEta = -1;
  cellCenterDPhi = -1;
  coneEnergy = -1;
  for (int i = 0; i < 7; i++) {
    hitEnergies[i] = -1;
    neighborHitEnergies[i] = -1;
  }
  hitsOverThresh = -1;
  vtxZ = -1;
  bremDepth = -1;
  nearestJetDr = -1;
  nearJetMaxE = -1;
  nearestJetE = -1;
  minCaloJetDr = -1;
  caloJetHcalE = -1;
  caloJetEcalE = -1;
  caloJetTotalE = -1;
  caloConeE = -1;
  caloConeHits = -1;
  hitsOverThresh = -1;
  genProbeMatchdE = -1;
  genProbeMatchdR = -1;
  probeMuondE = -1;
  globalMuonE = -1;
  nEle = 0;
  nExtraEle = 0;
  nExtraMuon = 0;
  probeCharge = 0;
  CSCdEta = -10;
  CSCdPhi = -10;
  stadEta = -10;
  stadPhi = -10;
  staNHits = 0;
  staChi2 = -1;
  staTransDCA = -1;
  staLongDCA = -1;
  finalMuE = -1.;
  disStaE = -1.;
  disStaDEoverE = 10.;
  nPUmean = -1.;
  pileupWeight = -1.;
  usedBin = -1;
  energyBinSF = 1;
}

bool SigEventInfo::passTriggers(const edm::Event& iEvent, 
		                edm::EDGetToken m_trigResultsToken,
				edm::EDGetToken m_trigEventToken,
				std::vector<std::string> m_muonPathsToPass) 
{ 
   bool passTriggers;
   edm::Handle<edm::TriggerResults> triggerResults;
   edm::Handle<trigger::TriggerEvent> trigEvent;
   iEvent.getByToken(m_trigEventToken,trigEvent);
   iEvent.getByToken(m_trigResultsToken, triggerResults);
   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);
   for (size_t i = 0; i < trigNames.size(); ++i) {
      const std::string& name = trigNames.triggerName(i);
      for (auto& pathName : m_muonPathsToPass) {
         if ((name.find(pathName) != std::string::npos)) {
            if (triggerResults->accept(i)) { 
               passTriggers = true;
               trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","","HLT"));
	       if(filterIndex<trigEvent->sizeFilters())
	       {
	          const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex);
		  const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
		  for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt)
		  {
		     const trigger::TriggerObject& obj = trigObjColl[*keyIt];
                     passingMuons.push_back(obj);
		  }
	       }
            }
         }
      }
   }

   return passTriggers;
}   

bool SigEventInfo::matchTagAndTrigger(const reco::Muon* selectedMuon)
{
   for(std::vector<trigger::TriggerObject>::const_iterator trigObject = passingMuons.begin(); trigObject!=passingMuons.end();++trigObject)
   {
      double dR=deltaR(trigObject->eta(),trigObject->phi(),selectedMuon->eta(),selectedMuon->phi());
      if(dR<0.1) return true;
   }
   return false;
}
