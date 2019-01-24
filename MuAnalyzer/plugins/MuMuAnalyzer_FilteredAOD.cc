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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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

//Triggers
#include "FWCore/Common/interface/TriggerNames.h"

#include "MuMu/MuMuAnalyzer/interface/Histograms.h"
#include "MuMu/MuMuAnalyzer/interface/CSC.h"
#include "MuMu/MuMuAnalyzer/interface/EventInfo.h"
#include "MuMu/MuMuAnalyzer/interface/Muons.h"
#include "MuMu/MuMuAnalyzer/interface/Tracks.h"


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
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, Histograms);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label;
    edm::EDGetToken m_trigResultsToken;
    std::vector<std::string> m_muonPathsToPass;

    bool m_isMC;
    bool m_runRandomTrackEfficiency;

    Histograms myHistograms;

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
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",true)),
  m_runRandomTrackEfficiency (iConfig.getUntrackedParameter<bool>("runRandomTrackEfficiency",false))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
  }
  m_trigResultsToken = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("trigResults"));
  m_muonPathsToPass   = iConfig.getParameter<std::vector<std::string> >("muonPathsToPass");

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  myHistograms.book(fs);

/*  Muons myMuons;
  Histograms myHistograms;
  EventInfo myEventInfo;
  Tracks myTracks;
  CSC myCSCs;*/


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

  Muons myMuons;
  EventInfo myEventInfo;
  Tracks myTracks;
  CSC myCSCs;

  myHistograms.h_eventCount->Fill(0.5);
  int cutProgress = 0;

  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
  cutProgress++;

  if(!myEventInfo.passTriggers(iEvent, m_trigResultsToken, m_muonPathsToPass)) return;
  cutProgress++;

  myMuons.SelectMuons(iEvent, m_recoMuonToken);

  myTracks.SelectTracks(iEvent, trackCollection_label);

  int anyMuonPass = 0;
  int anyTrackPass = 0;
  int anyGoodFittedVertices = 0;
  int anyPairPassMass = 0;
  int nMuonTrackCand = 0;
  int nTracksPairedPerMuon;
  int nTracksNoMuon = 0;
  int nMuonsPairedPerEvent = 0;

   // this wraps tracks with additional methods that are used in vertex-calculation
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  myCSCs.minTotalImpactParameter = 1000;
  myCSCs.minDR = 1000;
  myCSCs.minTotalImpactParameter_Muon = 1000;
  myCSCs.minDR_Muon = 1000;

//  bool foundMuonTrackPair = false;
  if(!m_runRandomTrackEfficiency){

    if(myMuons.selectedMuons.size() > 0){

      nTracksPairedPerMuon = 0;

      for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack ) {

          myTracks.tracksToVertex.clear();

	  //Track and Muon have opposite charge
          if (myMuons.highPtSelectedMuon->charge()==(*iTrack)->charge()) continue;
          if (!((*iTrack)->quality(Track::highPurity))) continue;


	  //Using the highest pT track that pairs with a muon
//	  if (foundMuonTrackPair) continue;

	  // make a pair of TransientTracks to feed the vertexer
	  if(myMuons.highPtSelectedMuon->isGlobalMuon()){
	    if(!myTracks.PairTracks(iTrack, myMuons.highPtSelectedMuon->globalTrack(), transientTrackBuilder)) continue;
	  }else if(myMuons.highPtSelectedMuon->isTrackerMuon()){
	    if(!myTracks.PairTracks(iTrack, myMuons.highPtSelectedMuon->innerTrack(), transientTrackBuilder)) continue;
          }else{
	    continue;
	  }

	  if (m_isMC){
	    if (!isTrackMatchedToMuon(iEvent, iTrack, myHistograms)) continue;
	  }

	  nTracksPairedPerMuon++;


          nMuonTrackCand++;

	  myHistograms.h_MuonTrackMass->Fill(myTracks.MuonTrackMass);

//	  foundMuonTrackPair = true;

	  myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, myTracks.one_momentum, myTracks.two_momentum, myTracks.tracksToVertex);

	  myHistograms.m_MinTotalImpactParameter->Fill(myCSCs.minTotalImpactParameter);
	  myHistograms.m_MinDR->Fill(myCSCs.minDR);

      }
      myHistograms.h_nTracksPairedPerMuon->Fill(nTracksPairedPerMuon);
      if (nTracksPairedPerMuon > 0)nMuonsPairedPerEvent++;
    }

    myHistograms.h_nMuonsPairedPerEvent->Fill(nMuonsPairedPerEvent);

    if(nMuonTrackCand > 0){

      myHistograms.PlotTrackDisappearance(myCSCs.TrackP, myCSCs.TrackEta, myCSCs.TrackPhi, myCSCs.minDR, myCSCs.minTotalImpactParameter, myCSCs.TrackP_dR, myCSCs.TrackEta_dR, myCSCs.TrackPhi_dR);
  
      myHistograms.m_histogram_TrackerTack_P->Fill(myCSCs.TrackP);

    }

    if (nMuonTrackCand > 0){
      myHistograms.h_nMuonTrackCand->Fill(nMuonTrackCand);
      myHistograms.h_nTracksNoMuon->Fill(nTracksNoMuon);
      std::cout <<"PASSES"<<std::endl;
    }else{
      std::cout << "FAILS"<<std::endl;
    }	


    if(anyMuonPass > 0){
      cutProgress++;
      if(anyTrackPass > 0){
        cutProgress++;
        if(anyGoodFittedVertices > 0){
	  cutProgress++;
	  if(anyPairPassMass > 0){
	    cutProgress++;
	  }
        }
      }
    }


    while(cutProgress > 0){
      myHistograms.h_cutProgress->Fill(cutProgress);
      cutProgress--;
    }

    //foundMuonTrackPair = false;

    //Looping over Muons and all Tracks to study extrapolating reco Muons to the CSCs.
    if(myMuons.selectedMuons.size() > 0){

      nTracksPairedPerMuon = 0;

    for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedTracks.begin(); iTrack != myTracks.selectedTracks.end(); ++iTrack ) {

         myTracks.tracksToVertex.clear();

         //Track and Muon have opposite charge
         if (myMuons.highPtSelectedMuon->charge()==(*iTrack)->charge()) continue;
         if (!((*iTrack)->quality(Track::highPurity))) continue;

         //Using the highest pT track that pairs with a muon
//         if (foundMuonTrackPair) continue;

         // make a pair of TransientTracks to feed the vertexer
         if(myMuons.highPtSelectedMuon->isGlobalMuon()){
           if(!myTracks.PairTracks(iTrack, myMuons.highPtSelectedMuon->globalTrack(), transientTrackBuilder)) continue;
         }else if(myMuons.highPtSelectedMuon->isTrackerMuon()){
           if(!myTracks.PairTracks(iTrack, myMuons.highPtSelectedMuon->innerTrack(), transientTrackBuilder)) continue;
         }else{
           continue;
         }

         myCSCs.ExtrapolateMuonToCSC(iEvent, iSetup, CSCSegment_Label, myMuons.highPtSelectedMuon, myTracks.one_momentum, myTracks.two_momentum, myTracks.tracksToVertex);

      }
    }

    if(fabs(myCSCs.MuonEta) > 1.653 && fabs(myCSCs.MuonEta) < 2.4){
      std::cout << "Plotting myCSCs.minDR_Muon: " << myCSCs.minDR_Muon << " myCSCs.minTotalImpactParameter_Muon: " << myCSCs.minTotalImpactParameter_Muon << std::endl;
      myHistograms.m_histogram_MuonTrack_P->Fill(myCSCs.MuonP);
      myHistograms.m_MinDR_Muon->Fill(myCSCs.minDR_Muon);
      myHistograms.m_MinTotalImpactParameterMuon->Fill(myCSCs.minTotalImpactParameter_Muon);
    }
  }else{
    //Making efficiencies for random tracks
    for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack ) {
       if (!(*iTrack)->quality(Track::highPurity)) continue;

       for(std::vector<const reco::Track*>::const_iterator iTrack_2nd = myTracks.selectedTracks.begin(); iTrack_2nd != myTracks.selectedTracks.end(); ++iTrack_2nd ) {

          myTracks.tracksToVertex.clear();

	  if (!(*iTrack_2nd)->quality(Track::highPurity)) continue;

	  if(!myTracks.PairTrackerTracks(iTrack, iTrack_2nd, transientTrackBuilder)) continue;

	  myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, myTracks.one_momentum, myTracks.two_momentum, myTracks.tracksToVertex);

	  myHistograms.m_MinTotalImpactParameter->Fill(myCSCs.minTotalImpactParameter);
          myHistograms.m_MinDR->Fill(myCSCs.minDR);
       }
    }
    myHistograms.PlotTrackDisappearance(myCSCs.TrackP, myCSCs.TrackEta, myCSCs.TrackPhi, myCSCs.minDR, myCSCs.minTotalImpactParameter, myCSCs.TrackP_dR, myCSCs.TrackEta_dR, myCSCs.TrackPhi_dR);
  }

}

bool MuMuAnalyzer_FilteredAOD::isTrackMatchedToMuon(const edm::Event& iEvent, std::vector<const reco::Track*>::const_iterator& Track, Histograms myHistograms){

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);

  bool matched = false;

  for(std::vector<reco::GenParticle>::const_iterator iGENparticle = genParticles->begin(); iGENparticle != genParticles->end(); iGENparticle++) {
     if(std::abs(iGENparticle->pdgId()) != 13) continue;

     if(deltaR((*Track)->eta(), (*Track)->phi(), iGENparticle->eta(), iGENparticle->phi()) > 0.1) continue;
     myHistograms.h_GENMuonTrackdR->Fill(deltaR((*Track)->eta(), (*Track)->phi(), iGENparticle->eta(), iGENparticle->phi()));
     matched = true;
  }
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
