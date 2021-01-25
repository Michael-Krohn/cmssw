// -*- C++ -*-
//
// Package:    DarkPhoton/MuPXAnalyzer
// Class:      MuPXAnalyzer
// 
/**\class MuPXAnalyzer MuPXAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/MuPXAnalyzer.cc

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

//for Standalone Muon Tracking
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

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

//for ECAL info
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

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

//Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DarkPhoton/MuAnalyzer/interface/MuPXHistograms.h"
#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"
#include "DarkPhoton/MuAnalyzer/interface/ECAL.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuPXAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuPXAnalyzer(const edm::ParameterSet&);
      ~MuPXAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, MuPXHistograms);
    bool MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, MuPXHistograms myHistograms, double FmuE, double minCSCdr);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetToken m_simTracksToken;
    edm::EDGetToken m_simVerticesToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label;
    edm::EDGetToken m_trigResultsToken;
    edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
    std::vector<std::string> m_muonPathsToPass;
    const reco::Track* selectedTrack;
    const reco::Muon* selectedMuon;
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label;
    edm::EDGetToken m_theSTAMuonLabel;
    bool m_isMC;
    double weight_;
    double standaloneE;
    MuPXHistograms myHistograms;

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
MuPXAnalyzer::MuPXAnalyzer(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  reducedEndcapRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
  reducedBarrelRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",true))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
    m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
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


MuPXAnalyzer::~MuPXAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void MuPXAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  Muons myMuons;
  EventInfo myEventInfo;
  Tracks myTracks;
  HCAL myHCAL;
  ECAL myECAL;

  myHistograms.ResetCutFlow();

  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  myHistograms.m_eventCount->Fill(0.5);
  //if(!founddpho){return;}
  myHistograms.IncCutFlow();
 
  myTracks.SelectTracks(iEvent, trackCollection_label);
  myMuons.SelectMuons(iEvent, m_recoMuonToken);
  myHistograms.m_NPassingTag->Fill(myMuons.selectedMuons.size());
 
  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);  

  bool Paired = false;
  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  int nPairedTracks =0;
  double pairVtxChi;
  double highpt=0;
  double muTrackMass=0;
  //Pair our candidate track with a tagging muon to reconstruct a Z
  if(myMuons.selectedMuons.size()>0)
  {
     for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
     {
        double hightrackpt=0;
        int nPairs=0;
        if((*iMuon)->pt()<highpt) continue;
        for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
        {
           if((*iMuon)->charge()==(*iTrack)->charge()) continue;
           //if(!(*iMuon)->isGlobalMuon()) continue;
           if((*iMuon)->isGlobalMuon())
           {
              if(myTracks.PairTracks(iTrack, (*iMuon)->globalTrack(), transientTrackBuilder))
              {
                Paired=true;
                nPairs++;
                if((*iTrack)->pt()>hightrackpt)
                {
                   highpt=(*iMuon)->pt();
                   pairVtxChi=myTracks.pairvertexchi;
                   muTrackMass=myTracks.MuonTrackMass;
                   selectedMuon = (*iMuon);
                   selectedTrack = (*iTrack); 
                }
              }
           } 
           else if((*iMuon)->isTrackerMuon())
           {
              if(myTracks.PairTracks(iTrack, (*iMuon)->innerTrack(), transientTrackBuilder))
              {
                Paired=true;
                nPairs++;
                if((*iTrack)->pt()>hightrackpt)
                {
                   highpt=(*iMuon)->pt();
                   pairVtxChi=myTracks.pairvertexchi;
                   muTrackMass=myTracks.MuonTrackMass;
                   selectedMuon = (*iMuon);
                   selectedTrack = (*iTrack); 
                }
              }
           }
        }
        if(nPairs>0){nPairedTracks=nPairs;} 
     }
  }
  
  myHistograms.IncCutFlow();
  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
 // myHistograms.IncCutFlow();
  if(myMuons.highmuonpt<26.){return;}
  myHistograms.IncCutFlow();
  myHistograms.m_NPassingProbe->Fill(nPairedTracks);
 
  edm::Handle <reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);
 
  if(Paired)
  {
     myHistograms.m_ProbeTrackIso->Fill(myTracks.GetIsolation(iEvent,trackCollection_label,selectedTrack->momentum().eta(),selectedTrack->momentum().phi(),0.3,vtxHandle,selectedTrack)/selectedTrack->pt());
     myHistograms.m_ProbeEcalIso->Fill(myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*selectedTrack)));
     myHistograms.m_ProbeHcalIso->Fill(myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*selectedTrack)));
     myHistograms.m_TagEta->Fill(selectedMuon->eta());
     myHistograms.m_TagPhi->Fill(selectedMuon->phi());
     myHistograms.m_TagEtaPhi->Fill(selectedMuon->eta(),selectedMuon->phi());
     myHistograms.m_TagProbeVtxChi->Fill(pairVtxChi);
     myHistograms.m_TagPt->Fill(selectedMuon->pt());
     myHistograms.m_ProbePt->Fill(selectedTrack->pt());
     myHistograms.m_ProbeEta->Fill(selectedTrack->eta());
     myHistograms.m_ProbePhi->Fill(selectedTrack->phi());
     myHistograms.m_ProbeEtaPhi->Fill(selectedTrack->eta(),selectedTrack->phi());
     myHistograms.m_MuonTrackMass->Fill(muTrackMass);
     myHistograms.IncCutFlow();
  }
}

bool MuPXAnalyzer::MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, MuPXHistograms myHistograms, double FmuE, double minCSCdr)
{
  bool matched = false;
  double mindE = 1.;
  if(myMuons.selectedMuons.size()==0){return matched;}
  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
  {
     //if(!(*iMuon)->isGlobalMuon()) {continue;}
     double dR = deltaR((*iMuon)->eta(),(*iMuon)->phi(),mother.eta(),mother.phi());
     double dEOverE =  std::abs(std::sqrt(pow((*iMuon)->p(),2)+pow(0.1056,2))-mother.E())/mother.E();
     if(dR > 0.2) continue;
     if(dEOverE<mindE)
     {
        mindE=dEOverE;
	matched=true;
     }
  }
  return matched;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuPXAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuPXAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuPXAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuPXAnalyzer);
