// -*- C++ -*-
//
// Package:    MuMu/MuFilter_AOD
// Class:      MuFilter_AOD
// 
/**\class MuFilter_AOD MuFilter_AOD_AOD.cc MuMu/MuFilter_AOD/plugins/MuFilter_AOD_AOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Krohn
//         Created:  Mon, 18 Jun 2018 21:22:23 GMT
//
//
// system include files
#include <memory>
#include <iomanip>
#include <iostream>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/TrackReco/interface/Track.h"

// for vertexing
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuFilter_AOD : public edm::stream::EDFilter<>  {
   public:
      explicit MuFilter_AOD(const edm::ParameterSet&);
      ~MuFilter_AOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<reco::Track>::const_iterator&);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trakCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;

    bool m_isMC;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuFilter_AOD::MuFilter_AOD(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trakCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",true))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
  }

}


MuFilter_AOD::~MuFilter_AOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
bool
MuFilter_AOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;
   using namespace pat;

  reco::Vertex bestVtx;
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(primaryVertices_Label, vertices);

  std::vector<const reco::Vertex*> PVertices;
  std::vector<reco::Vertex>::const_iterator firstGoodVertex = vertices->end();

  //Looping over vertices and verifying there is at least 1 good one. Might be unnecessary since vertex fitting is performed later.
  for (std::vector<reco::Vertex>::const_iterator it=vertices->begin(); it!=firstGoodVertex; ++it) {
    if (!it->isFake() && it->ndof()>4 && it->position().Rho()<2. && std::abs(it->position().Z())<24.) {
      if(firstGoodVertex == vertices->end()){
        firstGoodVertex = it;
        PVertices.push_back(&(*it));
      }
      break;
    }
  }

  if(firstGoodVertex == vertices->end()){
    std::cout<<"NO GOOD VERTEX" << std::endl;
    return false;
  }

  bestVtx = *(firstGoodVertex);

  int anyMuonPass = 0;
  int nMuonTrackCand = 0;
  int nTracksNoMuon = 0;

  float MuonTrackMass = 0.;

  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

   // this wraps tracks with additional methods that are used in vertex-calculation
   edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  //Looping over the reconstructed Muons
  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {
    std::cout << "CHECKING IF WE HAVE A LOOSE MUON" << std::endl;
    if (!(iMuon->isPFMuon() && (iMuon->isGlobalMuon() || iMuon->isTrackerMuon()))) continue; //Loose ID Muon requirements

    std::cout << "CHECKING IF MUON PASSES pT, eta" << std::endl;
    if (iMuon->pt() < 10 || fabs(iMuon->eta()) > 2.8) continue;
    anyMuonPass++;

    //Looping over tracks
    for(std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) {

        if (iMuon->charge()==iTrack->charge()) continue;

        if (!iTrack->quality(Track::highPurity)) continue;

        if (fabs(iTrack->eta()) > 2.8 || fabs(iTrack->eta()) < 1.653) continue;//Track required to be in the endcap

        // make a pair of TransientTracks to feed the vertexer
        std::vector<reco::TransientTrack> tracksToVertex;
        tracksToVertex.push_back(transientTrackBuilder->build(*iTrack));
	if(iMuon->isGlobalMuon()){
          tracksToVertex.push_back(transientTrackBuilder->build(iMuon->globalTrack()));
	}else if(iMuon->isTrackerMuon()){
          tracksToVertex.push_back(transientTrackBuilder->build(iMuon->innerTrack()));
	}else{
	  continue;
	}
        // try to fit these two tracks to a common vertex
        KalmanVertexFitter vertexFitter;
        CachingVertex<5> fittedVertex = vertexFitter.vertex(tracksToVertex);

        // some poor fits will simply fail to find a common vertex
        if (fittedVertex.isValid()  &&  fittedVertex.totalChiSquared() >= 0.  &&  fittedVertex.degreesOfFreedom() > 0) {
           // others we can exclude by their poor fit
           if (fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom() < 3.) {
              // important! evaluate momentum vectors AT THE VERTEX
              TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
              TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
              GlobalVector one_momentum = one_TSCP.momentum();
              GlobalVector two_momentum = two_TSCP.momentum();


              double total_energy = sqrt(one_momentum.mag2() + 0.106*0.106) + sqrt(two_momentum.mag2() + 0.106*0.106);
	      double total_px = one_momentum.x() + two_momentum.x();
	      double total_py = one_momentum.y() + two_momentum.y();
	      double total_pz = one_momentum.z() + two_momentum.z();
              MuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
           }else{
	      continue;
	   }
        }else{
	   continue;
	}

/*	TLorentzVector* TrackVector;
        TLorentzVector* MuonVector;
	TrackVector->SetPtEtaPhiM(iTrack->pt(),iTrack->eta(),iTrack->phi(),iTrack->m());*/
//        if ((iMuon->p4() + iTrack->p4()).mass() < 85 || (iMuon->p4() + iTrack->p4()).mass() > 105) continue;
	if (MuonTrackMass < 80 || MuonTrackMass > 100) continue;

        if (m_isMC){
	  if (!isTrackMatchedToMuon(iEvent, iTrack)) continue;//Matching the track to a GEN muon
	}

        nMuonTrackCand++;
//	MuonTrackMass = (iMuon->p4() + iTrack->p4()).mass();

	nTracksNoMuon=0.;
        for(std::vector<reco::Muon>::const_iterator iMuon2 = recoMuons->begin(); iMuon2 != recoMuons->end(); iMuon2++) {

           if (!(iMuon2->isPFMuon() && (iMuon2->isGlobalMuon() || iMuon2->isTrackerMuon()))) continue;
           if (fabs(iMuon2->eta()) > 2.8 || fabs(iMuon2->eta()) < 1.653) continue;
	   if(deltaR(iTrack->eta(), iTrack->phi(), iMuon2->eta(), iMuon2->phi()) < 0.1){
             nTracksNoMuon++;

	   }
	}

    }
  }

  if (nMuonTrackCand > 0){
    std::cout <<"PASSES"<<std::endl;
    return true;
  }else{
    std::cout << "FAILS"<<std::endl;
    return false;
  }		
}

bool MuFilter_AOD::isTrackMatchedToMuon(const edm::Event& iEvent, std::vector<reco::Track>::const_iterator& Track){

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);

  bool matched = false;
  int nMatchedParticles = 0;

  for(std::vector<reco::GenParticle>::const_iterator iGENparticle = genParticles->begin(); iGENparticle != genParticles->end(); iGENparticle++) {
     if(std::abs(iGENparticle->pdgId()) != 13) continue;

     if(deltaR(Track->eta(), Track->phi(), iGENparticle->eta(), iGENparticle->phi()) > 0.1) continue;
     nMatchedParticles++;
     matched = true;
  }
  return matched;

}
// ------------ method called once each job just before starting event loop  ------------
void 
MuFilter_AOD::beginStream(edm::StreamID)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuFilter_AOD::endStream() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuFilter_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuFilter_AOD);
