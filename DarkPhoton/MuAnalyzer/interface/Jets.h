#ifndef Jets_h
#define Jets_h

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/JetReco/interface/PFJet.h"


class Jets{

  public:
    Jets();

    void SelectJets(const edm::Event& iEvent, edm::EDGetToken m_jetsToken);

    std::vector<const reco::PFJet*> selectedJets;
};

#endif
