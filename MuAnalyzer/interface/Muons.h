#ifndef Muons_h
#define Muons_h

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


class Muons{

  public:
    Muons();

    void SelectMuons(const edm::Event& iEvent, edm::EDGetToken m_recoMuonToken);

    std::vector<const reco::Muon*> selectedMuons;
    const reco::Muon* highPtSelectedMuon;


};

#endif
