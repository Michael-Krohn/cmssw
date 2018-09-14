// -*- C++ -*-
//
// Package:    MuMu/MuMuAnalyzer_FilteredAOD
// Class:      MuMuAnalyzer_FilteredAOD
// 
/**\class MuMuAnalyzer_FilteredAOD MuMuAnalyzer_FilteredAOD_AOD.cc MuMu/MuMuAnalyzer_FilteredAOD/plugins/MuMuAnalyzer_FilteredAOD_AOD.cc

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
#include <TH2.h>

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

// for Global Coordinates
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuMuAnalyzer_FilteredAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuAnalyzer_FilteredAOD(const edm::ParameterSet&);
      ~MuMuAnalyzer_FilteredAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<reco::Track>::const_iterator&);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trakCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label;
 
    bool m_isMC;

    TH1F* h_MuonTrackMass;
    TH1F* h_GENMuonTrackdR;
    TH1F* h_NmatchedGenMuonTracks;
    TH1F* h_nMuonTrackCand;
    TH1F* h_CandTrackRecoMuonDR;
    TH1F* h_nTracksNoMuon;
    TH2F* m_histogram_fittedVertex_xy;
    TH1F* m_histogram_fittedVertex_z;
    TH2F* m_histogram_primaryVertex_xy;
    TH1F* m_histogram_primaryVertex_z;
    TH2F* m_histogram_differenceVertex_xy;
    TH1F* m_histogram_differenceVertex_z;
    TH1F* m_TransverseImpactParameter;
    TH1F* m_LongiudinalImpactParameter;
    TH1F* m_MinTransverseImpactParameter;
    TH1F* m_MinLongiudinalImpactParameter;
    TH2F* m_histogram_TrackerTrack_EtaPhi_negEta;
    TH2F* m_histogram_MuonTrack_EtaPhi_negEta;
    TH2F* m_histogram_TrackerTrack_EtaPhi_posEta;
    TH2F* m_histogram_MuonTrack_EtaPhi_posEta;

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
MuMuAnalyzer_FilteredAOD::MuMuAnalyzer_FilteredAOD(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trakCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",true))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
  }

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  h_GENMuonTrackdR = fs->make<TH1F>("DeltaR_GENMuonTrack", ";DeltaR(GENMuon,Track) ;Events", 30,  0.  , 0.2  );
  h_NmatchedGenMuonTracks = fs->make<TH1F>("nMatchedGenMuonsTracks", "; # of Matched Gen Muons To Tracks;Events", 10, -0.5  , 9.5  );
  h_MuonTrackMass = fs->make<TH1F>("MuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50  , 150  );
  h_nMuonTrackCand = fs->make<TH1F>("nMuonTrackCand", "; # of MuonTrack Candidates ;Events", 8, -0.5  , 7.5  );
  h_nTracksNoMuon = fs->make<TH1F>("nTrackNoMuon", "; # of MuonTrack Candidates ;Events", 8, -0.5  , 7.5  );
  h_CandTrackRecoMuonDR = fs->make<TH1F>("DeltaR_CandTrackRecoMuon", ";DeltaR(CandTrack,RecoMuon) ;Events", 50,  0.  , 0.01  );
  m_histogram_fittedVertex_xy = fs->make<TH2F>("fittedVertex_xy", "", 300, -3., 3., 300, -3., 3.);
  m_histogram_fittedVertex_z = fs->make<TH1F>("fittedVertex_z", "", 30, -30., 30.);
  m_histogram_primaryVertex_xy = fs->make<TH2F>("primaryVertex_xy", "", 300, -3., 3., 300, -3., 3.);
  m_histogram_primaryVertex_z = fs->make<TH1F>("primaryVertex_z", "", 30, -30., 30.);
  m_histogram_differenceVertex_xy = fs->make<TH2F>("differenceVertex_xy", "", 300, -3., 3., 300, -3., 3.);
  m_histogram_differenceVertex_z = fs->make<TH1F>("differenceVertex_z", "", 60, -10., 10.);
  m_TransverseImpactParameter = fs->make<TH1F>("TransverseImpactParameter", "", 300, -50, 50);
  m_LongiudinalImpactParameter = fs->make<TH1F>("LongiudinalImpactParameter", "", 300, -50, 50);
  m_MinTransverseImpactParameter = fs->make<TH1F>("MinTransverseImpactParameter", "", 300, -50, 50);
  m_MinLongiudinalImpactParameter = fs->make<TH1F>("MinLongiudinalImpactParameter", "", 300, -50, 50);
  m_histogram_TrackerTrack_EtaPhi_negEta = fs->make<TH2F>("TrackerTrack_EtaPhi_negEta", "", 100, -2.5, -1.63, 100, -3.2, 3.2); 
  m_histogram_MuonTrack_EtaPhi_negEta = fs->make<TH2F>("MuonTrack_EtaPhi_negEta", "", 100, -2.5, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta = fs->make<TH2F>("TrackerTrack_EtaPhi_posEta", "", 100, 1.63, 2.5, 100, -3.2, 3.2);
  m_histogram_MuonTrack_EtaPhi_posEta = fs->make<TH2F>("MuonTrack_EtaPhi_posEta", "", 100, 1.63, 2.5, 100, -3.2, 3.2);

}


MuMuAnalyzer_FilteredAOD::~MuMuAnalyzer_FilteredAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuAnalyzer_FilteredAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    return;
  }

  bestVtx = *(firstGoodVertex);

  int anyMuonPass = 0;
  int nMuonTrackCand = 0;
  int nTracksNoMuon = 0;

  float MuonTrackMass = 0.;

  bool TrackDissapears = true;

  double TrackEta = -999;
  double TrackPhi = -999;

  GlobalVector one_momentum;
  GlobalVector two_momentum; 

  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

   // this wraps tracks with additional methods that are used in vertex-calculation
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {
    std::cout << "CHECKING IF WE HAVE A LOOSE MUON" << std::endl;
    if (!(iMuon->isPFMuon() && (iMuon->isGlobalMuon() || iMuon->isTrackerMuon()))) continue;

    std::cout << "CHECKING IF MUON PASSES pT, eta" << std::endl;
    if (iMuon->pt() < 10 || fabs(iMuon->eta()) > 2.8) continue;
    anyMuonPass++;

    for(std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) {

        if (iMuon->charge()==iTrack->charge()) continue;

        if (!iTrack->quality(Track::highPurity)) continue;

        if (fabs(iTrack->eta()) > 2.8 || fabs(iTrack->eta()) < 1.653) continue;

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

              // this is a good vertex: fill histograms
//	      m_histogram_fittedVertex_xy->Fill(fittedVertex.position().x(), fittedVertex.position().y());
//              m_histogram_fittedVertex_z->Fill(fittedVertex.position().z());
                                     
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
           }
	   else{
	      MuonTrackMass = -10000;
	   }
        }

	if (MuonTrackMass < 80 || MuonTrackMass > 100) continue;

        if (m_isMC){
	  if (!isTrackMatchedToMuon(iEvent, iTrack)) continue;
	}

	TrackEta = one_momentum.eta();
	TrackPhi = one_momentum.phi();

        m_histogram_fittedVertex_xy->Fill(fittedVertex.position().x(), fittedVertex.position().y());
        m_histogram_fittedVertex_z->Fill(fittedVertex.position().z());
        m_histogram_differenceVertex_xy->Fill(bestVtx.position().x() - fittedVertex.position().x(), bestVtx.position().y() - fittedVertex.position().y());
	m_histogram_differenceVertex_z->Fill((bestVtx.position().z() - fittedVertex.position().z()));

        nMuonTrackCand++;
//	MuonTrackMass = (iMuon->p4() + iTrack->p4()).mass();
        h_MuonTrackMass->Fill(MuonTrackMass);

	double minTransverseImpactParameter = 1000;
	double minLongiudinalImpactParameter = 1000;
	if( TheCSCSegments.isValid() ){
	  std::cout << "Looping over CSC segments for a single Track" << std::endl;
          for(CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end(); iSegment++){
	     CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
	     if(iTrack->eta() < 0 && iDetId.endcap() == 1) continue;
	     if(iTrack->eta() > 0 && iDetId.endcap() == 2) continue;

/*	     std::cout << "iDetId.station: " << iDetId.station() << std::endl;
             std::cout << "iDetId.ring: " << iDetId.ring() << std::endl;
             std::cout << "iDetId.chamber: " << iDetId.chamber() << std::endl;*/

	     DetId TheDetUnitId(iSegment->cscDetId());
	     const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
//	     std::cout << "TheUnit->position().eta(): " << TheUnit->position().eta() << " TheUnit->position().phi(): " << TheUnit->position().phi() << std::endl;
	     LocalPoint TheLocalPosition = iSegment->localPosition();
	     const BoundPlane& TheSurface = TheUnit->surface();
	     GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);

	     TrajectoryStateClosestToPoint  traj = tracksToVertex[0].trajectoryStateClosestToPoint(TheGlobalPosition);

	     std::cout << "traj.perigeeParameters().transverseImpactParameter(): " << traj.perigeeParameters().transverseImpactParameter() << std::endl;
	     std::cout << "traj.perigeeError().transverseImpactParameterError(): " << traj.perigeeError().transverseImpactParameterError() << std::endl;
	     std::cout << "traj.perigeeParameters().longitudinalImpactParameter(): " << traj.perigeeParameters().longitudinalImpactParameter() << std::endl;
	     std::cout << "traj.perigeeError().longitudinalImpactParameterError(): " << traj.perigeeError().longitudinalImpactParameterError() << std::endl;
//	     std::cout << "TheLocalPosition.eta(): " << TheLocalPosition.eta() << " TheLocalPosition.phi(): " << TheLocalPosition.phi() << std::endl;
	     std::cout << "TheGlobalPosition.eta(): " << TheGlobalPosition.eta() << " TheGlobalPosition.phi(): " << TheGlobalPosition.phi() << std::endl;
	     std::cout << "iTrack->eta(): " << iTrack->eta() << " iTrack->phi(): " << iTrack->phi() << std::endl;
	     m_TransverseImpactParameter->Fill(traj.perigeeParameters().transverseImpactParameter());
	     m_LongiudinalImpactParameter->Fill(traj.perigeeParameters().longitudinalImpactParameter());
	     if (std::abs(minTransverseImpactParameter) > std::abs(traj.perigeeParameters().transverseImpactParameter()) && std::abs(minLongiudinalImpactParameter) > std::abs(traj.perigeeParameters().longitudinalImpactParameter())){
		minTransverseImpactParameter = traj.perigeeParameters().transverseImpactParameter();
		minLongiudinalImpactParameter = traj.perigeeParameters().longitudinalImpactParameter();
	     }
	  }
	}
        m_MinTransverseImpactParameter->Fill(minTransverseImpactParameter);
        m_MinLongiudinalImpactParameter->Fill(minLongiudinalImpactParameter);

	m_histogram_TrackerTrack_EtaPhi_negEta->Fill(TrackEta, TrackPhi);
	m_histogram_TrackerTrack_EtaPhi_posEta->Fill(TrackEta, TrackPhi);

	if(minTransverseImpactParameter > 1){// && minLongiudinalImpactParameter > 1){
	  TrackDissapears = true;
	}else{
	  TrackDissapears = false;
	}

	if(TrackDissapears== false){
	  m_histogram_MuonTrack_EtaPhi_negEta->Fill(TrackEta, TrackPhi);
	  m_histogram_MuonTrack_EtaPhi_posEta->Fill(TrackEta, TrackPhi);

	}
	
/*	nTracksNoMuon=0.;
        for(std::vector<reco::Muon>::const_iterator iMuon2 = recoMuons->begin(); iMuon2 != recoMuons->end(); iMuon2++) {

           if (!(iMuon2->isPFMuon() && (iMuon2->isGlobalMuon() || iMuon2->isTrackerMuon()))) continue;
           if (fabs(iMuon2->eta()) > 2.8 || fabs(iMuon2->eta()) < 1.653) continue;
//           h_CandTrackRecoMuonDR->Fill(deltaR(iTrack->eta(), iTrack->phi(), iMuon2->eta(), iMuon2->phi()));
	   if(deltaR(iTrack->eta(), iTrack->phi(), iMuon2->eta(), iMuon2->phi()) < 0.1){
	     h_CandTrackRecoMuonDR->Fill(deltaR(iTrack->eta(), iTrack->phi(), iMuon2->eta(), iMuon2->phi()));
             nTracksNoMuon++;

	   }
	}*/

    }
  }

  if (nMuonTrackCand > 0){
    h_nMuonTrackCand->Fill(nMuonTrackCand);
    h_nTracksNoMuon->Fill(nTracksNoMuon);
    m_histogram_primaryVertex_xy->Fill(bestVtx.position().x(), bestVtx.position().y());
    m_histogram_primaryVertex_z->Fill(bestVtx.position().z());
    std::cout <<"PASSES"<<std::endl;
  }else{
    std::cout << "FAILS"<<std::endl;
  }		
}

bool MuMuAnalyzer_FilteredAOD::isTrackMatchedToMuon(const edm::Event& iEvent, std::vector<reco::Track>::const_iterator& Track){

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);

  bool matched = false;
  int nMatchedParticles = 0;

  for(std::vector<reco::GenParticle>::const_iterator iGENparticle = genParticles->begin(); iGENparticle != genParticles->end(); iGENparticle++) {
     if(std::abs(iGENparticle->pdgId()) != 13) continue;

     if(deltaR(Track->eta(), Track->phi(), iGENparticle->eta(), iGENparticle->phi()) > 0.1) continue;
     nMatchedParticles++;
     h_GENMuonTrackdR->Fill(deltaR(Track->eta(), Track->phi(), iGENparticle->eta(), iGENparticle->phi()));
     matched = true;
  }
  h_NmatchedGenMuonTracks->Fill(nMatchedParticles);
  return matched;

}
// ------------ method called once each job just before starting event loop  ------------
void 
MuMuAnalyzer_FilteredAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuAnalyzer_FilteredAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuAnalyzer_FilteredAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuAnalyzer_FilteredAOD);
