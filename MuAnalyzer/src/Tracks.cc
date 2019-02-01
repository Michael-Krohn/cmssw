#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "FWCore/Framework/interface/Event.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" 

Tracks::Tracks(){}

void Tracks::SelectTracks(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label){

  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);

  for(std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) {
      selectedTracks.push_back(&(*iTrack));

      if (fabs(iTrack->eta()) > 2.4 || fabs(iTrack->eta()) < 1.653) continue;

      selectedEndcapTracks.push_back(&(*iTrack));
  }
}

bool Tracks::PairTracks(std::vector<const reco::Track*>::const_iterator& Track, const reco::TrackRef MuonTrack, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder){

  tracksToVertex.push_back(transientTrackBuilder->build(*Track));
  tracksToVertex.push_back(transientTrackBuilder->build(*MuonTrack));

  // try to fit these two tracks to a common vertex
  KalmanVertexFitter vertexFitter;
  fittedVertex = vertexFitter.vertex(tracksToVertex);

  // some poor fits will simply fail to find a common vertex
  if (fittedVertex.isValid()  &&  fittedVertex.totalChiSquared() >= 0.  &&  fittedVertex.degreesOfFreedom() > 0) {
    // others we can exclude by their poor fit
    if (fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom() < 3.) {
      // important! evaluate momentum vectors AT THE VERTEX
      TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
      TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
      one_momentum = one_TSCP.momentum();
      two_momentum = two_TSCP.momentum();

      double total_energy = sqrt(one_momentum.mag2() + 0.106*0.106) + sqrt(two_momentum.mag2() + 0.106*0.106);
      double total_px = one_momentum.x() + two_momentum.x();
      double total_py = one_momentum.y() + two_momentum.y();
      double total_pz = one_momentum.z() + two_momentum.z();
      MuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
    }else{
      MuonTrackMass = -10000;
    }

    if (MuonTrackMass < 80 || MuonTrackMass > 100) return false;

  }else{
    return false;
  }

  return true;

}


bool Tracks::PairTrackerTracks(std::vector<const reco::Track*>::const_iterator& Track, std::vector<const reco::Track*>::const_iterator& Track_2nd, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder){

  tracksToVertex.push_back(transientTrackBuilder->build(*Track));
  tracksToVertex.push_back(transientTrackBuilder->build(*Track_2nd));

  // try to fit these two tracks to a common vertex
  KalmanVertexFitter vertexFitter;
  fittedVertex = vertexFitter.vertex(tracksToVertex);

  // some poor fits will simply fail to find a common vertex
  if (fittedVertex.isValid()  &&  fittedVertex.totalChiSquared() >= 0.  &&  fittedVertex.degreesOfFreedom() > 0) {
    // others we can exclude by their poor fit
    if (fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom() < 3.) {
      TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
      TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
      one_momentum = one_TSCP.momentum();
      two_momentum = two_TSCP.momentum();
      return true;
    }else{
      return false;
    }
  
  }else{
    return false;
  }
  

}
