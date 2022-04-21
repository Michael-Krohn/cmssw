#ifndef EventInfo_h
#define EventInfo_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"  
#include "TH1F.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class EventInfo{
  public:
    EventInfo();
    bool goodPrimaryVertex(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label);
    bool passTriggers(const edm::Event& iEvent, edm::EDGetToken m_trigResultsToken, std::vector<std::string> m_muonPathsToPass);
    reco::Vertex bestVtx;
    double eventWeight;
    double diMuonMass;
    double minDR;
    double ProbeEta;
    double ProbePhi;
    double heConeEnergy;
    double probeHcalIso;
    double tagHcalIso;
    int probeiEta;
    int centerInCone;
    int centerOutOfCone;
    double minDR_muonHCAL;
    double minDR_randomHCAL;
    int hitDepth_muonHCAL;
    int hitDepth_randomHCAL;
    int coneHits;
    double coneEnergy;
    int randomConeHits;
    double randomConeEnergy;
    int hitsOverThresh;
    int randomHitsOverThresh;
    int cellsFound;
    double hitEnergies[7];
    double randHitEnergies[7];
    double neighborHitEnergies[7];
    double neighborHitEnergiesPhi[7];
    double cellCenterDEta;
    double cellCenterDPhi;
};
#endif
