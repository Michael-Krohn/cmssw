#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include <iostream>

EventInfo::EventInfo() {
  eventWeight = 1;
  diMuonMass = -1;
  heConeEnergy = -1;
  probeHcalIso = -1;
  tagHcalIso = -1;
  minDR = -1;
  ProbeEta = 10;
  ProbePhi = 10;
  probeiEta = 0;
  centerInCone = -1;
  centerOutOfCone = -1;
  minDR_muonHCAL = -1;
  minDR_randomHCAL = -1;
  expectedHits = 0;
  coneHits = -1;
  coneEnergy = -1;
  randomConeHits = -1;
  randomConeEnergy = -1;
  hitsOverThresh = -1;
  randomHitsOverThresh = -1;
  cellsFound = -1;
  for(int depth=0;depth<7;depth++)
  {
    xShiftEnergies[depth] = -1;
    yShiftEnergies[depth] = -1;
    zShiftEnergies[depth] = -1;
    xShifts[depth] = 10;
    yShifts[depth] = 10;
    zShifts[depth] = 10;
  }
}

bool EventInfo::goodPrimaryVertex(const edm::Event& iEvent,
                                  edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label) {
  using namespace reco;

  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(primaryVertices_Label, vertices);

  std::vector<const reco::Vertex*> PVertices;
  std::vector<reco::Vertex>::const_iterator firstGoodVertex = vertices->end();

  for (std::vector<reco::Vertex>::const_iterator it = vertices->begin(); it != firstGoodVertex; ++it) {
    if (!it->isFake() && it->ndof() > 4 && it->position().Rho() < 2. && std::abs(it->position().Z()) < 24.) {
      if (firstGoodVertex == vertices->end()) {
        firstGoodVertex = it;
        PVertices.push_back(&(*it));
      }
      break;
    }
  }

  if (firstGoodVertex == vertices->end()) {
    std::cout << "NO GOOD VERTEX" << std::endl;
    return false;
  } else {
    bestVtx = *(firstGoodVertex);
    return true;
  }
}

bool EventInfo::passTriggers(const edm::Event& iEvent,
                             edm::EDGetToken m_trigResultsToken,
                             std::vector<std::string> m_muonPathsToPass) {
  bool passTriggers;

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(m_trigResultsToken, triggerResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);
  for (size_t i = 0; i < trigNames.size(); ++i) {
    const std::string& name = trigNames.triggerName(i);
    for (auto& pathName : m_muonPathsToPass) {
      if ((name.find(pathName) != std::string::npos)) {
        if (triggerResults->accept(i)) {
          passTriggers = true;
        }
      }
    }
  }

  return passTriggers;
}
