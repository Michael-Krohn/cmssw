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

//Jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DarkPhoton/MuAnalyzer/interface/MuPXHistograms.h"
#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "DarkPhoton/MuAnalyzer/interface/Jets.h"
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
    bool MatchTrackToMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, Muons myMuons);

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
    std::vector<const reco::Track*> pairedTracks;
    const reco::Track* selectedTrack;
    const reco::Muon* selectedMuon;
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label;
    edm::EDGetToken m_theSTAMuonLabel;
    edm::EDGetToken m_pfJetCollection_label;
    bool m_isMC;
    double weight_;
    double standaloneE;
    MuPXHistograms myHistograms;
    MuPXHistograms muProbe;
    MuPXHistograms nonMuonProbe;
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
  m_pfJetCollection_label(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PFJets"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",false))
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

  //myHistograms.book(fs->mkdir("allEvents"));
  //muProbe.book(fs->mkdir("muProbe"));
  nonMuonProbe.book(fs->mkdir("nonMuonProbe"));
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
  MuPXEventInfo info;
  Tracks myTracks;
  Jets myJets;
  HCAL myHCAL;
  ECAL myECAL;

  //myHistograms.ResetCutFlow();

  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  if(m_isMC)
  {
     edm::Handle<GenEventInfoProduct> eventInfo;
     iEvent.getByToken(m_genInfoToken, eventInfo);
     weight_  = eventInfo->weight();
  }
  else
  {
     weight_=1;
  }
  //if(!founddpho){return;}
  info.cutProgress++;
 
  myTracks.SelectTracks(iEvent, trackCollection_label);
  myMuons.SelectMuons(iEvent, m_recoMuonToken);
  myJets.SelectJets(iEvent, m_pfJetCollection_label);
//  myHistograms.m_NPassingTag->Fill(myMuons.selectedMuons.size(),weight_);
 
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
  double MuonTrackMass=0;


  //Looping over vertices and verifying there is at least 1 good one. Might be unnecessary since vertex fitting is performed later.

  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
     {
        double hightrackpt=0;
        int nPairs=0;
        //if((*iMuon)->pt()<highpt) continue;
        for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
        {
           if((*iMuon)->charge()==(*iTrack)->charge()) continue;
           //if(!(*iMuon)->isGlobalMuon()) continue;
           if((*iMuon)->pt()<highpt){continue;}
           if((*iMuon)->isGlobalMuon())
           {
              if(myTracks.PairTracks(iTrack, (*iMuon)->globalTrack(), transientTrackBuilder))
              {
                Paired=true;
                if(nPairs==0){pairedTracks.clear();}
                nPairs++;
                pairedTracks.push_back((*iTrack));
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
                if(nPairs==0){pairedTracks.clear();}
                nPairs++;
                pairedTracks.push_back((*iTrack));
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
  
  info.cutProgress++;
  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
  if(myMuons.highmuonpt<26.){return;}
  info.cutProgress++;
  info.eventWeight=weight_;
  info.nPassingProbe=nPairedTracks;
  info.nPassingTag=myMuons.selectedMuons.size();
  info.paired=Paired;
  edm::Handle <reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);
 
  if(Paired)
  {
     info.cutProgress++;
     info.muonTrackMass=muTrackMass;
     info.probeTrack=selectedTrack;
     info.tagMuon=selectedMuon;
     info.probeTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,selectedTrack->momentum().eta(),selectedTrack->momentum().phi(),0.3,vtxHandle,selectedTrack);
     info.probeEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*selectedTrack));
     info.probeHcalIso=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*selectedTrack));
     info.tagProbeVtxChi=pairVtxChi;
     //Study multiple paired tracks
     double largestDR = 0;
     double drSum = 0;
     int drcount = 0;
     for(std::vector<const reco::Track*>::const_iterator iTrack = pairedTracks.begin(); iTrack != pairedTracks.end(); ++iTrack)
     {
        for(std::vector<const reco::Track*>::const_iterator iTrack2 = pairedTracks.begin(); iTrack2 != pairedTracks.end(); ++iTrack2)
        {
           if(iTrack!=iTrack2)
           {
             double DR = deltaR((*iTrack)->eta(),(*iTrack)->phi(),(*iTrack2)->eta(),(*iTrack2)->phi());
             drcount++;
             drSum+=DR;
             if(DR>largestDR){largestDR=DR;} 
           }
        }
     }
     if(drcount>0)
     {
        info.smallestCone=largestDR/2.;
        info.averageDr=drSum/drcount;
     }
     info.tagTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,selectedMuon->momentum().eta(),selectedMuon->momentum().phi(),0.3,vtxHandle,selectedMuon->innerTrack());
     info.tagEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build((*selectedMuon).globalTrack()));
     info.tagHcalIso=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build((*selectedMuon).globalTrack()));

    //Jet analysis
    int njets=0;
    double nearestJetDr=-1;
    double nearestJetE;
    for(std::vector<const reco::PFJet*>::const_iterator it = myJets.selectedJets.begin(); it!=myJets.selectedJets.end();++it)
    {
      njets++;
      double jetDr=deltaR(selectedTrack->eta(),selectedTrack->phi(),(*it)->eta(),(*it)->phi());
      if(jetDr<nearestJetDr||nearestJetDr<0)
      {
        nearestJetDr=jetDr;
        nearestJetE=(*it)->energy();
      }
      //myHistograms.m_JetPt->Fill((*it)->pt(),weight_);
    }
    info.nearestJetE=nearestJetE;
    info.nearestJetDr=nearestJetDr;
    info.nJets=njets;
  }
  //myHistograms.FillHists(info);
  bool matched = MatchTrackToMuon(iEvent, selectedTrack, myMuons);
  //if(matched){muProbe.FillHists(info);}
  if(!matched){nonMuonProbe.FillHists(info);}

}

bool MuPXAnalyzer::MatchTrackToMuon(const edm::Event& iEvent,const reco::Track* selectedTrack, Muons myMuons)
{
  bool matched = false;
  double mindPt = 1.;
  if(myMuons.selectedMuons.size()==0){return matched;}
  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
  {
     if(!(*iMuon)->isGlobalMuon()) {continue;}
     double dR = deltaR((*iMuon)->eta(),(*iMuon)->phi(),selectedTrack->eta(),selectedTrack->phi());
     double dPtOverPt =  std::abs(((*iMuon)->pt()-selectedTrack->pt())/selectedTrack->pt());
     if(dR > 0.2) continue;
     if(dPtOverPt<mindPt)
     {
        mindPt=dPtOverPt;
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
