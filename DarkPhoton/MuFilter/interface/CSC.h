#ifndef CSC_H
#define CSC_H

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class CSC{
  public:
    CSC();

    void ExtrapolateTrackToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, std::vector<const reco::Track*>::const_iterator& iTrack, GlobalVector one_momentum, std::vector<reco::TransientTrack> tracksToVertex, GlobalPoint VertexPosition);
    void ExtrapolateTrackToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, const reco::Track* iTrack, reco::TransientTrack tracksToVertex, GlobalPoint VertexPosition);
   void ExtrapolateMuonToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, const reco::Muon* iMuon, GlobalVector two_momentum, std::vector<reco::TransientTrack> tracksToVertex);

    GlobalPoint MuonGlobalPoint;
    GlobalPoint TrackGlobalPoint;
    GlobalPoint TrackVertex;
    GlobalVector CSCTraj;
    double CSCChiSq;
    double MuonEta;
    double MuonPhi;
    double MuonP;
    double MuonEta_dR;
    double MuonPhi_dR;
    double MuonP_dR;
    double minDR;
    double minDR_Muon;
    double minTotalImpactParameter;
    double minTotalImpactParameter_Muon;
    double TrackP;
    double TrackEta;
    double TrackPhi;
    double TrackP_dR;
    double TrackEta_dR;
    double TrackPhi_dR;
};

#endif
