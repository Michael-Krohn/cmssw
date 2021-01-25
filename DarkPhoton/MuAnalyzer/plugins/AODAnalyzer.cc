// -*- C++ -*-
//
// Package:    DarkPhoton/AODAnalyzer
// Class:      AODAnalyzer
// 
/**\class AODAnalyzer AODAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/AODAnalyzer.cc

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

class AODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AODAnalyzer(const edm::ParameterSet&);
      ~AODAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, MCHistograms);
    bool MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, MCHistograms myHistograms);
    bool MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, MCHistograms myHistograms, double FmuE, double minCSCdr);
    void MatchSimTrks(math::XYZTLorentzVectorD mother, math::XYZTLorentzVectorD tagmu, Muons myMuons, Tracks myTracks, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, MCHistograms myHistograms);

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
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;
    edm::EDGetToken m_theSTAMuonLabel;
    bool m_isMC;
    double weight_;
    double standaloneE;
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
AODAnalyzer::AODAnalyzer(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
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


AODAnalyzer::~AODAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void AODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  Muons myMuons;
  EventInfo myEventInfo;
  Tracks myTracks;
  HCAL myHCAL;

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
  math::XYZTLorentzVectorD FinalMu;
  
  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);  

  bool SiblingMatched = false;
  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  int nMatches =0;
  bool backupmatch = false;
  double tageta,tagpt, pairVtxChi;
  //Pair out dark brem candidate track with a tagging muon to reconstruct a Z
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
  {
     if(myMuons.selectedMuons.size()<1){continue;}
     if(myMuons.highPtSelectedMuon->charge()==(*iTrack)->charge()) continue;
     if(!myMuons.highPtSelectedMuon->isGlobalMuon()) continue;
     if(myTracks.PairTracks(iTrack, myMuons.highPtSelectedMuon->globalTrack(), transientTrackBuilder))
     {
       SiblingMatched=true;
       nMatches++;
       pairVtxChi=myTracks.pairvertexchi;
       tageta=myMuons.highPtSelectedMuon->eta();
       tagpt=myMuons.highPtSelectedMuon->pt();
       selectedTrack = (*iTrack); 
     } 
  }
/*  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
  {
     for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
     {
        if((*iMuon)->charge()==(*iTrack)->charge()) continue;
        if(!(*iMuon)->isGlobalMuon()) continue;
        if(myTracks.PairTracks(iTrack, (*iMuon)->globalTrack(), transientTrackBuilder))
        {
          SiblingMatched=true;
          nMatches++;
          tageta=(*iMuon)->eta();
          tagpt=(*iMuon)->pt();
          selectedTrack = (*iTrack); 
        } 
     }
     nselectedmuons++;
  }*/

  myHistograms.m_nSelectedMuons->Fill(myMuons.selectedMuons.size());
  if(!SiblingMatched)
  {  
     for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
     {
        if(myMuons.selectedMuons.size()<1){continue;}
        if(myMuons.secondhighPtSelectedMuon->charge()==(*iTrack)->charge()) continue;
        if(!myMuons.secondhighPtSelectedMuon->isGlobalMuon()) continue;
        if(myTracks.PairTracks(iTrack, myMuons.secondhighPtSelectedMuon->globalTrack(), transientTrackBuilder))
        {
          SiblingMatched=true;
          nMatches++;
          tageta=myMuons.secondhighPtSelectedMuon->eta();
          tagpt=myMuons.secondhighPtSelectedMuon->pt();
          selectedTrack = (*iTrack); 
          backupmatch = true;
        } 
     }
  }
  myHistograms.IncCutFlow();
  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
 // myHistograms.IncCutFlow();
  if(myMuons.highmuonpt<26.){return;}
  myHistograms.IncCutFlow();
  myHistograms.m_NMatchedMuons->Fill(nMatches);
 
  edm::Handle <reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);
 
  if(SiblingMatched)
  {
     myHistograms.m_EventsWithDpho->Fill(0.5);
     myHistograms.m_TrackIsolation->Fill(myTracks.GetIsolation(iEvent,trackCollection_label,selectedTrack->momentum().eta(),selectedTrack->momentum().phi(),0.3,vtxHandle,selectedTrack)/selectedTrack->pt());
     myHistograms.m_TagMuonEta->Fill(tageta);
     myHistograms.m_PairVtxChi->Fill(pairVtxChi);
     if(fabs(tageta)>1.6)
     {
        myHistograms.m_TagMuonEta->Fill(selectedTrack->eta());
     }
     myHistograms.m_TagMuonPt->Fill(tagpt);
     if(!backupmatch)
     {
        myHistograms.m_NoBackupTagMuonEta->Fill(tageta);
        myHistograms.m_NoBackupTagMuonPt->Fill(tagpt);
     }
     myHistograms.m_DiMuonMass->Fill(myTracks.MuonTrackMass);
     myHistograms.IncCutFlow();
     GlobalPoint TrackGlobalPoint = GlobalPoint(GlobalPoint::Polar(selectedTrack->theta(),selectedTrack->phi(),50));
   // this wraps tracks with additional methods that are used in vertex-calculation
     if(!myHCAL.FindMuonHits(iEvent, iSetup, HBHERecHit_Label, TrackGlobalPoint, myHistograms, standaloneE, weight_, selectedTrack->dsz(), selectedTrack->charge(), transientTrackBuilder->build(selectedTrack))){return;}
     myHistograms.m_HECellsFound->Fill(myHCAL.CellsFound);
  }
}

bool AODAnalyzer::MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, MCHistograms myHistograms)
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

void AODAnalyzer::MatchSimTrks(math::XYZTLorentzVectorD mother, math::XYZTLorentzVectorD tagmu, Muons myMuons, Tracks myTracks, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, MCHistograms myHistograms)
{
  const reco::Track* selectedTrack;
  const reco::Muon* selectedMuon;
  double mindE = 10;
  bool mumatched = false;
  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
  {
     double dR = deltaR((*iMuon)->eta(),(*iMuon)->phi(),tagmu.eta(),tagmu.phi());
     double dEOverE =  std::abs(std::sqrt(pow((*iMuon)->p(),2)+pow(0.1056,2))-tagmu.E())/tagmu.E();
     if(dR > 0.2) continue;
     if(dEOverE<mindE)
     {
        mindE=dEOverE;
	mumatched=true;
        selectedMuon=(*iMuon);
     }
  }
  mindE = 10;
  bool trackmatched = false;
  if(mother.pt()>0&&mumatched) //I never actually check if a mother track was found, so do this to throw out cases where it wasn't.
  {
     myHistograms.m_SimProbeTrackPt->Fill(mother.pt());
     myHistograms.m_SimProbeTrackEta->Fill(mother.eta());
  }
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
  {
     double dR = deltaR((*iTrack)->eta(),(*iTrack)->phi(),mother.eta(),mother.phi());
     double dEOverE =  std::abs(std::sqrt(pow((*iTrack)->p(),2)+pow(0.1056,2))-mother.E())/mother.E();
     if(dR > 0.2) continue;
     if(dEOverE<mindE)
     {
        mindE=dEOverE;
	trackmatched=true;
        selectedTrack=(*iTrack);
     }
  }
  if(trackmatched&&mumatched)
  {
    myTracks.PairSimTracks(selectedTrack, selectedMuon->globalTrack(), transientTrackBuilder);
    myHistograms.m_SimVtxChi->Fill(myTracks.simVtxChi);
    myHistograms.m_SimInvMass->Fill(myTracks.SimMuonTrackMass);
  }
}

bool AODAnalyzer::MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, MCHistograms myHistograms, double FmuE, double minCSCdr)
{
  bool matched = false;
  double mindE = 1.;
  double globalmuonE;
  double mindR = 0.1;
  double standalonedR = 0.;
  if(myMuons.selectedMuons.size()==0){return matched;}
  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
  {
     if(((*iMuon)->isTrackerMuon())&&((*iMuon)->isGlobalMuon())){myHistograms.m_TrackerVGlobal->Fill(3);}
     else if((*iMuon)->isGlobalMuon())
     {
        myHistograms.m_TrackerVGlobal->Fill(1);
     }
     else if((*iMuon)->isTrackerMuon())
     {
        myHistograms.m_TrackerVGlobal->Fill(2);
     }
     else{myHistograms.m_TrackerVGlobal->Fill(4);}
     //if(!(*iMuon)->isGlobalMuon()) {continue;}
     double dR = deltaR((*iMuon)->eta(),(*iMuon)->phi(),mother.eta(),mother.phi());
     double dEOverE =  std::abs(std::sqrt(pow((*iMuon)->p(),2)+pow(0.1056,2))-mother.E())/mother.E();
     if(dR > 0.2) continue;
     if(dEOverE<mindE)
     {
        mindR=dR;
        mindE=dEOverE;
	matched=true;
        if((*iMuon)->isAValidMuonTrack(reco::Muon::OuterTrack))
        {
          standaloneE = (*iMuon)->outerTrack()->p();
          standalonedR = deltaR((*iMuon)->outerTrack()->eta(),(*iMuon)->outerTrack()->phi(),(*iMuon)->eta(),(*iMuon)->phi());
        }
        else
        {
          standaloneE = -1.;
          standalonedR = 10.;
        }
        globalmuonE = (*iMuon)->p();
     }
  }
  if(matched)
  {
    myHistograms.m_MuonMotherdE->Fill(mindE);
    myHistograms.m_MuonMotherdR->Fill(mindR);
    myHistograms.m_StandaloneMuonE->Fill(standaloneE);
    myHistograms.m_GlobalMuonE->Fill(globalmuonE);
    myHistograms.m_StandalonedE->Fill((standaloneE-FmuE)/FmuE);
    if(minCSCdr>0.2)
    {
      myHistograms.m_NHitStandaloneMuonE->Fill(standaloneE);
      myHistograms.m_NHitStandalonedE->Fill((standaloneE-FmuE)/FmuE);      
      myHistograms.m_NHitStandalonedR->Fill(standalonedR);
      myHistograms.m_BigCSCLocations->Fill(mother.eta(),mother.phi());
      if(standalonedR<0.2){myHistograms.m_BigCSCdrEta->Fill(mother.eta());}
    }
  }
  else
  {
    standaloneE=0;
    myHistograms.m_NonMatchedE->Fill(FmuE);     
  }
  return matched;
}
// ------------ method called once each job just before starting event loop  ------------
void 
AODAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AODAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AODAnalyzer);
