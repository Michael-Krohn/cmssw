#ifndef TRACKS_H
#define TRACKS_H

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

class Tracks{
  public:
    Tracks();

    void SelectTracks(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label);
    bool PairTracks(std::vector<const reco::Track*>::const_iterator& Track, const reco::TrackRef MuonTrack, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder);
    bool PairTrackerTracks(std::vector<const reco::Track*>::const_iterator& Track, std::vector<const reco::Track*>::const_iterator& Track_2nd, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder);
    double GetIsolation(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, double eta, double phi, double conesize, double vetopt);
    std::vector<reco::TransientTrack> tracksToVertex;
    std::vector<const reco::Track*> selectedEndcapTracks;
    std::vector<const reco::Track*> selectedTracks;    
    GlobalVector one_momentum;
    GlobalVector two_momentum;
    float MuonTrackMass;
    CachingVertex<5> fittedVertex;

};

#endif
