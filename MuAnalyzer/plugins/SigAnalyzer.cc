// -*- C++ -*-
//
// Package:    DarkPhoton/SigAnalyzer
// Class:      SigAnalyzer
// 
/**\class SigAnalyzer SigAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/SigAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Revering (adapted from analyzer developed by Michael Krohn)
//         Created:  Friday, 2 Jan 2020 13:34:23 GMT
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

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
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

#include "DarkPhoton/MuAnalyzer/interface/MCHistograms.h"
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

class SigAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SigAnalyzer(const edm::ParameterSet&);
      ~SigAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, MCHistograms);
    bool MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, MCHistograms myHistograms);
    bool MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, MCHistograms myHistograms);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetToken m_simTracksToken;
    edm::EDGetToken m_simVerticesToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label;
    edm::EDGetToken m_trigResultsToken;
    std::vector<std::string> m_muonPathsToPass;
    const reco::Track* selectedTrack;
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;

    bool m_isMC;
    bool m_isSig;
    bool m_runRandomTrackEfficiency;
    double bremDepth;

    MCHistograms myHistograms;

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
SigAnalyzer::SigAnalyzer(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",true)),
  m_isSig (iConfig.getUntrackedParameter<bool>("isSig",true)),
  m_runRandomTrackEfficiency (iConfig.getUntrackedParameter<bool>("runRandomTrackEfficiency",false))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
    m_simTracksToken = consumes<edm::SimTrackContainer> (iConfig.getParameter<edm::InputTag>("g4SimHits"));
    m_simVerticesToken = consumes<edm::SimVertexContainer> (iConfig.getParameter<edm::InputTag>("g4SimHits"));
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


SigAnalyzer::~SigAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void SigAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<SimTrackContainer> simtracks;
  iEvent.getByToken(m_simTracksToken, simtracks);
  edm::Handle<SimVertexContainer> simvertices;
  iEvent.getByToken(m_simVerticesToken, simvertices);
  myHistograms.ResetCutFlow();
  int dphovtx = -1;
  bool founddpho=false;
  double dphoE = 0;
  math::XYZTLorentzVectorD dpho, mother;
  //find dark brem vertex
  for(SimTrackContainer::const_iterator isimtrk = simtracks->begin(); isimtrk!=simtracks->end(); ++isimtrk)
  {
     if(!isimtrk->noVertex()&&isimtrk->type()==9994)
     {
        const SimVertex &vtx = (*simvertices)[isimtrk->vertIndex()];
	if(double(vtx.position().z())>388&&double(vtx.position().rho())<268&&double(vtx.position().rho())<(1.2857*vtx.position().z()-302.8))
	{
	   bremDepth = vtx.position().r();
           myHistograms.m_DBremLocation->Fill(vtx.position().z(),vtx.position().rho());
	   dpho=isimtrk->momentum();
	   dphovtx=isimtrk->vertIndex();
	   dphoE=isimtrk->momentum().E();
	   myHistograms.m_dphoEnergy->Fill(dphoE);
	   founddpho=true;
	}
     }
  }

  //Find deflected muon
  double FmuE = 0;
  bool foundsib=false;
  if(m_isSig){printf("isSig is true.\n");}
  else{printf("isSig is false.\n");}
  math::XYZTLorentzVectorD FinalMu;
  for(SimTrackContainer::const_iterator isimtrk = simtracks->begin(); isimtrk!=simtracks->end(); ++isimtrk)
  {
     if(!isimtrk->noVertex())
     {
        if(isimtrk->vertIndex()==dphovtx&&isimtrk->type()!=9994)
	{
	   FmuE=isimtrk->momentum().E();
           FinalMu = isimtrk->momentum();
	   double iE = sqrt(pow(FinalMu.x()+dpho.x(),2)+pow(FinalMu.y()+dpho.y(),2)+pow(FinalMu.z()+dpho.z(),2)+pow(0.1056,2));
           if(m_isSig){mother.SetPxPyPzE(FinalMu.x()+dpho.x(),FinalMu.y()+dpho.y(),FinalMu.z()+dpho.z(),iE);}
	   else{mother.SetPxPyPzE(FinalMu.x(),FinalMu.y(),FinalMu.z(),FmuE);}
	   foundsib=true;
	}
     }
  }

  myHistograms.m_eventCount->Fill(0.5);
  myHistograms.IncCutFlow();
  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
  if(!founddpho){return;}
  myHistograms.IncCutFlow();
  if(!foundsib){return;}
  myHistograms.IncCutFlow();
  if(!m_isMC)
  {
     if(!myEventInfo.passTriggers(iEvent, m_trigResultsToken, m_muonPathsToPass)) return;
  }
  myHistograms.IncCutFlow();

  myMuons.SelectMuons(iEvent, m_recoMuonToken);  
  
  //Match the dark brem to a track
  bool Matched = MatchMotherToTrack(iEvent, mother, trackCollection_label, myHistograms);
  if(Matched)
  {
     myHistograms.IncCutFlow();
   // this wraps tracks with additional methods that are used in vertex-calculation
     edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
     GlobalPoint FillerPoint; 
     myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, selectedTrack,transientTrackBuilder->build(selectedTrack),  FillerPoint);
     GlobalPoint TrackGlobalPoint = GlobalPoint(GlobalPoint::Polar(selectedTrack->theta(),selectedTrack->phi(),bremDepth));
     myHCAL.FindMuonHits(iEvent, iSetup, HBHERecHit_Label, TrackGlobalPoint, myHistograms);
     double dotpro = mother.x()*FinalMu.x()+mother.y()*FinalMu.y()+mother.z()*FinalMu.z();
     myHistograms.m_DeflectingAngle->Fill(std::acos(dotpro/mother.P()/FinalMu.P()));
     myHistograms.m_FractionalELost->Fill(FmuE/(FmuE+dphoE));
     myHistograms.m_finalMuE->Fill(FmuE);
     if(m_isSig){myHistograms.m_initialMuE->Fill(FmuE+dphoE);}
     else{myHistograms.m_initialMuE->Fill(FmuE);}
     myHistograms.m_finalMuPt->Fill(FinalMu.Pt());
     myHistograms.m_initialMuPt->Fill(mother.Pt());
     myHistograms.m_CSC_dR->Fill(myCSCs.minDR);
     bool MuonMatched = MatchTrackToMuon(iEvent,mother,myMuons,myHistograms);
     if(MuonMatched){myHistograms.IncCutFlow();}    
  }

//Pair track with a muon to make a Z

}

bool SigAnalyzer::MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, MCHistograms myHistograms)
{
  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);

  bool matched = false;
  double mindR;
  double mindE = 0.2;
  for(std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) 
  {
     if (fabs(iTrack->eta()) > 2.4 || fabs(iTrack->eta()) < 1.653) continue;
     if (iTrack->pt()<26) continue;
     double dR = deltaR(iTrack->eta(),iTrack->phi(), mother.eta(), mother.phi());
     double dEOverE = std::abs(std::sqrt(pow(iTrack->outerP(),2)+pow(0.1056,2))-mother.E())/mother.E();
     if(dR > 0.1) continue;
     if(dEOverE<mindE)
     {
        mindR=dR;
	mindE=dEOverE;
        selectedTrack = (&(*iTrack));
	matched = true;
     }
  }
  if(matched)
  {
    myHistograms.m_TrackMotherdE->Fill(mindE);
    myHistograms.m_TrackMotherdR->Fill(mindR);
  }

  return matched;
}

bool SigAnalyzer::MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, MCHistograms myHistograms)
{
  bool matched = false;
  double mindE = 1.;
  double mindR = 0.1;
  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
  {
     if(!((*iMuon)->isGlobalMuon())) 
     {
        myHistograms.m_TrackerVGlobal->Fill(1);
        continue;
     }
     if(!((*iMuon)->isTrackerMuon())){myHistograms.m_TrackerVGlobal->Fill(2);}
     else{myHistograms.m_TrackerVGlobal->Fill(3);}
     double dR = deltaR((*iMuon)->eta(),(*iMuon)->phi(),mother.eta(),mother.phi());
     double dEOverE =  std::abs(std::sqrt(pow((*iMuon)->p(),2)+pow(0.1056,2))-mother.E())/mother.E();
     if(dR > 0.1) continue;
     if(dEOverE<mindE)
     {
        mindR=dR;
        mindE=dEOverE;
	matched=true;
     }
  }
  if(matched)
  {
    myHistograms.m_MuonMotherdE->Fill(mindE);
    myHistograms.m_MuonMotherdR->Fill(mindR);
  }
  return matched;
}
// ------------ method called once each job just before starting event loop  ------------
void 
SigAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SigAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SigAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SigAnalyzer);
