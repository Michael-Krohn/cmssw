#ifndef EVENTINFO_H
#define EVENTINFO_H

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


class EventInfo{

  public:
    EventInfo();

    bool goodPrimaryVertex(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label);
    bool passTriggers(const edm::Event& iEvent, edm::EDGetToken m_trigResultsToken, std::vector<std::string> m_muonPathsToPass);

    reco::Vertex bestVtx;


};

#endif

