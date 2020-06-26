// -*- C++ -*-
//
// Package:    DarkPhoton/MuAnalyzer
// Class:      MuAnalyzer
// 
/**\class MuAnalyzer MuAnalyzer_AOD.cc DarkPhoton/MuAnalyzer/plugins/MuAnalyzer_AOD.cc

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
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

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

//for HCAL info
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

//Triggers
#include "FWCore/Common/interface/TriggerNames.h"

#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"


// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuAnalyzer(const edm::ParameterSet&);
      ~MuAnalyzer();

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
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;
    edm::EDGetToken m_simTracksToken;


    bool m_isMC;
    bool m_isSig;
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
MuAnalyzer::MuAnalyzer(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",false)),
  m_isSig (iConfig.getUntrackedParameter<bool>("isSig",false)),
  m_runRandomTrackEfficiency (iConfig.getUntrackedParameter<bool>("runRandomTrackEfficiency",false))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
    if(m_isSig){m_simTracksToken=consumes<edm::SimTrackContainer> (iConfig.getParameter<edm::InputTag>("g4SimHits"));}
  }
  else
  {
     m_trigResultsToken = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("trigResults"));
     m_muonPathsToPass   = iConfig.getParameter<std::vector<std::string> >("muonPathsToPass");
  }
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  myHistograms.book(fs);

}


MuAnalyzer::~MuAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
MuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  Muons myMuons;
  EventInfo myEventInfo;
  Tracks myTracks;
  CSC myCSCs;
  HCAL myHCAL;

  myHistograms.m_eventCount->Fill(0.5);
  int cutProgress = 0;

  //if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
  cutProgress++;
  if(!m_isMC)
  {
     if(!myEventInfo.passTriggers(iEvent, m_trigResultsToken, m_muonPathsToPass)) return;
  }
  math::XYZTLorentzVectorD dpho;
  if(m_isSig)
  {
    edm::Handle<SimTrackContainer> simtracks;
    iEvent.getByToken(m_simTracksToken, simtracks);
    bool dphofound=false;
    for(SimTrackContainer::const_iterator isimtrk = simtracks->begin(); isimtrk!=simtracks->end(); ++isimtrk)
    {
       if(!isimtrk->noVertex()&&isimtrk->type()==9994){dphofound=true;}
       dpho=isimtrk->momentum();
    }
    if(!dphofound){return;}
  }

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
	  myHistograms.m_MuonEtaDist->Fill(myMuons.highPtSelectedMuon->eta());

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

	  myHistograms.m_MuonTrackMass->Fill(myTracks.MuonTrackMass);

	  myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, myTracks.one_momentum, myTracks.tracksToVertex, myTracks.fittedVertex.position());
          if(myCSCs.minTotalImpactParameter>10){continue;}
          nTracksPairedPerMuon++;

          nMuonTrackCand++;

	  myHistograms.m_MinTotalImpactParameter->Fill(myCSCs.minTotalImpactParameter);
	  myHistograms.m_MinDR->Fill(myCSCs.minDR);

      }
      myHistograms.m_nTracksPairedPerMuon->Fill(nTracksPairedPerMuon);
      if (nTracksPairedPerMuon > 0)nMuonsPairedPerEvent++;
    }

    myHistograms.m_nMuonsPairedPerEvent->Fill(nMuonsPairedPerEvent);

    if(nMuonTrackCand > 0){
      cutProgress++;
/*      if(fabs(myCSCs.TrackEta) > 2.3 && myCSCs.minTotalImpactParameter < 5){
	std::cout << "PASSING TRACK ON EDGE OF DETECTOR" << std::endl;
	std::cout << "myCSCs.TrackEta: " << myCSCs.TrackEta << " myCSCs.TrackPhi: " << myCSCs.TrackPhi << std::endl;
	std::cout << "myCSCs.TrackEta_dR: " << myCSCs.TrackEta_dR << " myCSCs.TrackPhi_dR: " << myCSCs.TrackPhi_dR << std::endl;
	std::cout << "myCSCs.minTotalImpactParameter: " << myCSCs.minTotalImpactParameter << " myCSCs.minDR: " << myCSCs.minDR << std::endl;
      }*/
      myHistograms.PlotTrackDisappearance(myCSCs.TrackP, myCSCs.TrackEta, myCSCs.TrackPhi, myCSCs.minDR, myCSCs.minTotalImpactParameter, myCSCs.TrackP_dR, myCSCs.TrackEta_dR, myCSCs.TrackPhi_dR);
      if(myCSCs.minDR < 0.1){myHCAL.CheckHCAL(iEvent, iSetup, HBHERecHit_Label);}
      myHistograms.m_histogram_TrackerTack_P->Fill(myCSCs.TrackP);
      myHistograms.m_nMuonTrackCand->Fill(nMuonTrackCand);
      myHistograms.m_nTracksNoMuon->Fill(nTracksNoMuon);
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

    //Looping over Muons and all Tracks to study extrapolating reco Muons to the CSCs.
    if(myMuons.selectedEndcapMuons.size() > 0){

      nTracksPairedPerMuon = 0;
      for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedTracks.begin(); iTrack != myTracks.selectedTracks.end(); ++iTrack ) {

         myTracks.tracksToVertex.clear();

         //Track and Muon have opposite charge
         if (myMuons.highPtSelectedEndcapMuon->charge()==(*iTrack)->charge()) continue;
         if (!((*iTrack)->quality(Track::highPurity))) continue;
         //Using the highest pT track that pairs with a muon
         // make a pair of TransientTracks to feed the vertexer
	 continue;
         if(myMuons.highPtSelectedEndcapMuon->isGlobalMuon()){
           if(!myTracks.PairTracks(iTrack, myMuons.highPtSelectedEndcapMuon->innerTrack(), transientTrackBuilder)) continue;
         }else if(myMuons.highPtSelectedEndcapMuon->isTrackerMuon()){
           if(!myTracks.PairTracks(iTrack, myMuons.highPtSelectedEndcapMuon->innerTrack(), transientTrackBuilder)) continue;
         }else{
           continue;
         }
         if(myTracks.GetIsolation(iEvent, trackCollection_label,myTracks.two_momentum.eta(),myTracks.two_momentum.phi(),0.3,myTracks.two_momentum.perp())>3.0){continue;}
         myCSCs.ExtrapolateMuonToCSC(iEvent, iSetup, CSCSegment_Label, myMuons.highPtSelectedEndcapMuon, myTracks.two_momentum, myTracks.tracksToVertex);
      }
    }
    if(fabs(myCSCs.MuonEta_dR) > 1.653 && fabs(myCSCs.MuonEta_dR) < 2.4) {
      cutProgress++;
      myHistograms.m_histogram_MuonTrack_P->Fill(myCSCs.MuonP_dR);
      myHistograms.m_MinDR_Muon->Fill(myCSCs.minDR_Muon);
      myHistograms.m_MinTotalImpactParameterMuon->Fill(myCSCs.minTotalImpactParameter_Muon);
      }

    double MatchedP, MatchedMinDr;
    MatchedMinDr = 10;
    GlobalPoint MatchedGlobalPoint, VertexPosition;
    bool TrackHCAL = true;
    if(TrackHCAL)
    {
       MatchedP = myCSCs.TrackP_dR;
       MatchedMinDr = myCSCs.minDR;
       MatchedGlobalPoint = myCSCs.TrackGlobalPoint;
       VertexPosition = myCSCs.TrackVertex;
    }
    else
    {
       MatchedP = myCSCs.MuonP_dR;
       MatchedMinDr = myCSCs.minDR_Muon;
       MatchedGlobalPoint = myCSCs.MuonGlobalPoint;
    }
     
    if(MatchedMinDr < 0.1)
    {
      cutProgress++;
      double randphi = MatchedGlobalPoint.phi()+2.0;
      if(randphi>ROOT::Math::Pi()){randphi-=2*ROOT::Math::Pi();}
      if(MatchedGlobalPoint.eta()<0&&randphi<-0.9&&randphi>-1.6){randphi=randphi-4.0+2*ROOT::Math::Pi();}
      GlobalPoint RandGlobalPoint(GlobalPoint::Polar(MatchedGlobalPoint.theta(),randphi,MatchedGlobalPoint.mag()));
      bool GoodRand=true;
      if(myTracks.GetIsolation(iEvent, trackCollection_label,RandGlobalPoint.eta(),RandGlobalPoint.phi(),0.3,0)>3.0){GoodRand=false;}
      double minDR_MuonHCAL = myHCAL.MuonMindR(iEvent, iSetup, HBHERecHit_Label, MatchedGlobalPoint);
      if((pow(MatchedGlobalPoint.eta()-dpho.eta(),2)+pow(MatchedGlobalPoint.phi()-dpho.phi(),2))>0.01){return;}
      myHCAL.HitsPlots(iEvent, iSetup, HBHERecHit_Label, MatchedGlobalPoint, RandGlobalPoint, GoodRand, myHistograms, MatchedP, VertexPosition);
      myHistograms.m_MinDR_MuonHCAL->Fill(minDR_MuonHCAL);
      myHistograms.m_HitEnergy_MinDR_MuonHCAL->Fill(myHCAL.MuonHitEnergy);
      myHistograms.PlotCSCHits(iEvent,iSetup,CSCSegment_Label);
      //myHistograms.PlotHCALHits(iEvent,iSetup,HBHERecHit_Label);
    }

    while(cutProgress > 0)
    {
      myHistograms.m_cutProgress->Fill(cutProgress);
      cutProgress--;
    }

  }else{
    bool useFirstTrackToPair = false;
    //Making efficiencies for random tracks
    cout << "Making random tracks!!\n";
    for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack ) {

       if (useFirstTrackToPair) continue;

       if (!(*iTrack)->quality(Track::highPurity)) continue;

       for(std::vector<const reco::Track*>::const_iterator iTrack_2nd = myTracks.selectedTracks.begin(); iTrack_2nd != myTracks.selectedTracks.end(); ++iTrack_2nd ) {
	  myCSCs.minDR = 1000;
	  myCSCs.minTotalImpactParameter = 1000;

          myTracks.tracksToVertex.clear();

	  if (!(*iTrack_2nd)->quality(Track::highPurity)) continue;

	  if(!myTracks.PairTrackerTracks(iTrack_2nd, iTrack, transientTrackBuilder)) continue;

          if(myTracks.GetIsolation(iEvent, trackCollection_label,myTracks.one_momentum.eta(),myTracks.one_momentum.phi(),0.3,myTracks.one_momentum.perp())>3.0){continue;}

	  myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack_2nd, myTracks.one_momentum, myTracks.tracksToVertex, myTracks.fittedVertex.position());
          
	  myHistograms.m_MinTotalImpactParameter->Fill(myCSCs.minTotalImpactParameter);
          myHistograms.m_MinDR->Fill(myCSCs.minDR);

	  myHistograms.PlotTrackDisappearance(myCSCs.TrackP, myCSCs.TrackEta, myCSCs.TrackPhi, myCSCs.minDR, myCSCs.minTotalImpactParameter, myCSCs.TrackP_dR, myCSCs.TrackEta_dR, myCSCs.TrackPhi_dR);
       }
       useFirstTrackToPair = true;
    }
  }
}

bool MuAnalyzer::isTrackMatchedToMuon(const edm::Event& iEvent, std::vector<const reco::Track*>::const_iterator& Track, Histograms myHistograms)
{
  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);

  bool matched = false;

  for(std::vector<reco::GenParticle>::const_iterator iGENparticle = genParticles->begin(); iGENparticle != genParticles->end(); iGENparticle++) {
     if(std::abs(iGENparticle->pdgId()) != 13) continue;

     if(deltaR((*Track)->eta(), (*Track)->phi(), iGENparticle->eta(), iGENparticle->phi()) > 0.1) continue;
     myHistograms.m_GENMuonTrackdR->Fill(deltaR((*Track)->eta(), (*Track)->phi(), iGENparticle->eta(), iGENparticle->phi()));
     matched = true;
  }
  return matched;
}
// ------------ method called once each job just before starting event loop  ------------
void 
MuAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuAnalyzer::endJob() 
{
   myHistograms.Normalize();  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuAnalyzer);
