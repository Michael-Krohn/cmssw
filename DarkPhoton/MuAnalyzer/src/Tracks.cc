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
  highendcaptrackpt = 0;
  for(std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) {
      selectedTracks.push_back(&(*iTrack));

      if (fabs(iTrack->eta()) > 2.4 || fabs(iTrack->eta()) < 1.3) continue;
      if (iTrack->pt()<15) continue;
      if (iTrack->pt()>highendcaptrackpt)
      {
         highptendcaptrack = &(*iTrack);
         highendcaptrackpt = iTrack->pt();
      }
      selectedEndcapTracks.push_back(&(*iTrack));
  }
}

double Tracks::GetIsolation(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,double eta, double phi, double conesize, edm::Handle<reco::VertexCollection> vtxHandle, const reco::Track* MainTrack)
{
   bool foundtrack = false;
   unsigned int vtxindex = 0;
   unsigned int trackindex = 0;
   double Isolation = 0;
   //double Isolation = 0;a
   for(unsigned int i=0;i<vtxHandle->size();i++)
   {
      reco::VertexRef vtx(vtxHandle, i);
      if(!vtx->isValid()){continue;}
      for(unsigned int j=0; j<vtx->tracksSize();j++)
      {
         if(vtx->trackRefAt(j)->pt()==MainTrack->pt()&&vtx->trackRefAt(j)->eta()==MainTrack->eta()&&vtx->trackRefAt(j)->phi()==MainTrack->phi())
         {
            vtxindex = i;
            trackindex = j;
            foundtrack = true;
         }
      }
   }
   if(!foundtrack)
   {
      return -1;
   }
 
   reco::VertexRef primaryVtx(vtxHandle,vtxindex);
    
   for(unsigned int i=0;i<primaryVtx->tracksSize();i++)
   {
      if(i==trackindex){continue;}
      reco::TrackBaseRef secondarytrack = primaryVtx->trackRefAt(i);
      double dphi = fabs(phi-secondarytrack->phi());
      if(dphi>ROOT::Math::Pi()) dphi -= 2*ROOT::Math::Pi();
      double Dr = sqrt( pow(eta-secondarytrack->eta(),2.0) + pow(dphi,2.0));
      if(Dr>conesize){continue;}
       Isolation += secondarytrack->pt();
   }
   return Isolation;
}

double Tracks::GetIsolation(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,double eta, double phi, double conesize, edm::Handle<reco::VertexCollection> vtxHandle, const reco::TrackRef MainTrack)
{
   bool foundtrack = false;
   unsigned int vtxindex = 0;
   unsigned int trackindex = 0;
   double Isolation = 0;
   //double Isolation = 0;a
   for(unsigned int i=0;i<vtxHandle->size();i++)
   {
      reco::VertexRef vtx(vtxHandle, i);
      if(!vtx->isValid()){continue;}
      for(unsigned int j=0; j<vtx->tracksSize();j++)
      {
         if(fabs(vtx->trackRefAt(j)->pt()-MainTrack->pt())<0.01)
        // if((*vtx->trackRefAt(j))==(*MainTrack))
         {
            vtxindex = i;
            trackindex = j;
            foundtrack = true;
         }
      }
   }
   if(!foundtrack)
   {
      return 10000;
   }
 
   reco::VertexRef primaryVtx(vtxHandle,vtxindex);
    
   for(unsigned int i=0;i<primaryVtx->tracksSize();i++)
   {
      if(i==trackindex){continue;}
      reco::TrackBaseRef secondarytrack = primaryVtx->trackRefAt(i);
      double dphi = fabs(phi-secondarytrack->phi());
      if(dphi>ROOT::Math::Pi()) dphi -= 2*ROOT::Math::Pi();
      double Dr = sqrt( pow(eta-secondarytrack->eta(),2.0) + pow(dphi,2.0));
      if(Dr>conesize){continue;}
       Isolation += secondarytrack->pt();
   }
   return Isolation;
}

bool Tracks::PairTracks(std::vector<const reco::Track*>::const_iterator& Track, const reco::TrackRef MuonTrack, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder){
  tracksToVertex.clear();
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
      //if(sqrt(pow(one_momentum.x(),2)+pow(one_momentum.y(),2))<15){return false;}
      two_momentum = two_TSCP.momentum();

      double total_energy = sqrt(one_momentum.mag2() + 0.106*0.106) + sqrt(two_momentum.mag2() + 0.106*0.106);
      double total_px = one_momentum.x() + two_momentum.x();
      double total_py = one_momentum.y() + two_momentum.y();
      double total_pz = one_momentum.z() + two_momentum.z();
      MuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
    }else{
      MuonTrackMass = -10000;
    }

    if (MuonTrackMass < 50 || MuonTrackMass > 150) return false;
  }else{
    return false;
  }
  pairvertexchi=fittedVertex.totalChiSquared()/fittedVertex.degreesOfFreedom();
  return true;

}

bool Tracks::PairTracks(const reco::Track* Track, const reco::TrackRef  MuonTrack, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder){
  tracksToVertex.clear();
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
      //if(sqrt(pow(one_momentum.x(),2)+pow(one_momentum.y(),2))<15){return false;}
      two_momentum = two_TSCP.momentum();

      double total_energy = sqrt(one_momentum.mag2() + 0.106*0.106) + sqrt(two_momentum.mag2() + 0.106*0.106);
      double total_px = one_momentum.x() + two_momentum.x();
      double total_py = one_momentum.y() + two_momentum.y();
      double total_pz = one_momentum.z() + two_momentum.z();
      MuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
    }else{
      MuonTrackMass = -10000;
    }

    if (MuonTrackMass < 50 || MuonTrackMass > 150) return false;
  }else{
    return false;
  }

  return true;

}

void Tracks::PairSimTracks(const reco::Track* Track, const reco::TrackRef  MuonTrack, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder){

  tracksToVertex.push_back(transientTrackBuilder->build(*Track));
  tracksToVertex.push_back(transientTrackBuilder->build(*MuonTrack));
  // try to fit these two tracks to a common vertex
  KalmanVertexFitter vertexFitter;
  fittedVertex = vertexFitter.vertex(tracksToVertex);

  // some poor fits will simply fail to find a common vertex
  if (fittedVertex.isValid()  &&  fittedVertex.totalChiSquared() >= 0.  &&  fittedVertex.degreesOfFreedom() > 0) 
  {
    // others we can exclude by their poor fit
    simVtxChi = fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom(); 
  }
  else{simVtxChi = -1;}
      // important! evaluate momentum vectors AT THE VERTEX
  TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
  TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
  one_momentum = one_TSCP.momentum();
  two_momentum = two_TSCP.momentum();
  
  double total_energy = sqrt(one_momentum.mag2() + 0.106*0.106) + sqrt(two_momentum.mag2() + 0.106*0.106);
  double total_px = one_momentum.x() + two_momentum.x();
  double total_py = one_momentum.y() + two_momentum.y();
  double total_pz = one_momentum.z() + two_momentum.z();
  SimMuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
  SimProbeTrackPt = one_momentum.perp();
  SimProbeTrackEta = one_momentum.eta();
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
