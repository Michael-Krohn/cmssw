// -*- C++ -*-
//
// Package:    DarkPhoton/MuPlusXFilter_HardBrem_RAW_RECO
// Class:      MuPlusXFilter_HardBrem_RAW_RECO
// 
/**\class MuPlusXFilter_HardBrem_RAW_RECO MuPlusXFilter_HardBrem_RAW_RECO.cc DarkPhoton/MuFilter/plugins/MuPlusXFilter_HardBrem_RAW_RECO.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Revering (adapted from filterr developed by Michael Krohn)
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
#include "Math/VectorUtil.h"

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
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1Trigger/interface/CaloSpare.h"
#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
//Triggers
#include "FWCore/Common/interface/TriggerNames.h"


//Custom Classes
#include "DarkPhoton/MuFilter/interface/MuPXHistograms.h"
#include "DarkPhoton/MuFilter/interface/CSC.h"
#include "DarkPhoton/MuFilter/interface/EventInfo.h"
#include "DarkPhoton/MuFilter/interface/Muons.h"
#include "DarkPhoton/MuFilter/interface/Tracks.h"
#include "DarkPhoton/MuFilter/interface/Jets.h"
#include "DarkPhoton/MuFilter/interface/HCAL.h"
#include "DarkPhoton/MuFilter/interface/ECAL.h"
#include "DarkPhoton/MuFilter/interface/ParticleListDrawer.h"

// class declaration
//

// If the filterr does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuPlusXFilter_HardBrem_RAW_RECO : public edm::stream::EDFilter<>  {
   public:
      explicit MuPlusXFilter_HardBrem_RAW_RECO(const edm::ParameterSet&);
      ~MuPlusXFilter_HardBrem_RAW_RECO();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, MuPXHistograms);
    bool MatchTrackToMuon(const edm::Event& iEvent, const reco::Track* selTrack, Muons myMuons, edm::EDGetToken m_recoMuonToken, MuPXEventInfo& info);
    bool MatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double *genInfo);
    bool LooseMatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double *genInfo, std::vector<int> &nearbyIDs, std::vector<int> &highPtIds);
    int PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, HCAL myHCAL,  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetToken m_recoPhotonToken;
    edm::EDGetToken m_simTracksToken;
    edm::EDGetToken m_simVerticesToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label;
    edm::EDGetToken m_trigResultsToken;
    std::vector<std::string> m_muonPathsToPass;
    std::vector<const reco::Track*> pairedTracks;
    const reco::Track* selectedTrack;
    const reco::Muon* selectedMuon;
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label;
    edm::EDGetToken m_theSTAMuonLabel;
    edm::EDGetToken m_pfJetCollection_label;
    edm::EDGetToken m_caloTower_label;
    edm::EDGetToken m_caloJet_label;
    edm::EDGetTokenT<reco::GsfElectronCollection> m_electronsToken;
    edm::EDGetToken m_PUInfoToken;
    edm::EDGetTokenT<reco::CandidateView> candToken_;
    bool m_isMC;
    bool m_signalRegion;
    bool m_makeTuples;
    bool m_makeSplits;
    double weight_;
    double nPUmean_;
    double standaloneE;
    double pairVtxChi;
    double muTrackMass;
    double probeLIP;
    double probeTIP;
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
MuPlusXFilter_HardBrem_RAW_RECO::MuPlusXFilter_HardBrem_RAW_RECO(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  reducedEndcapRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
  reducedBarrelRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",false)),
  m_signalRegion(iConfig.getUntrackedParameter<bool>("signalRegion",false)),
  m_makeTuples(iConfig.getUntrackedParameter<bool>("makeTuples",false)),
  m_makeSplits (iConfig.getUntrackedParameter<bool>("makeSplits",true))
{
   //now do what ever initialization is needed


}


MuPlusXFilter_HardBrem_RAW_RECO::~MuPlusXFilter_HardBrem_RAW_RECO()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
bool 
MuPlusXFilter_HardBrem_RAW_RECO::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  Muons myMuons;
  EventInfo myEventInfo;
  MuPXEventInfo info;
  Tracks myTracks;
  HCAL myHCAL;
  ECAL myECAL;
  CSC myCSCs;

  myCSCs.minTotalImpactParameter = 1000;
  myCSCs.minDR = 1000;
  myCSCs.minTotalImpactParameter_Muon = 1000;
  myCSCs.minDR_Muon = 1000;

  info.cutProgress=0;

  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  weight_=1;
  nPUmean_=0;


  info.eventWeight = weight_;
  info.cutProgress++;
  edm::Handle <reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);

  edm::Handle<std::vector<reco::Track> > thePATTrackHandle; 
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);  
  
  std::vector<double> HCALisolationInfo;
  int muProgress = myMuons.SelectMuons(iEvent, m_recoMuonToken);
  if(muProgress<4){return false;}
  bool Paired = false;
  int maxTrackProgress=0;
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  for(std::vector<const reco::Muon*>::const_iterator muon = myMuons.selectedMuons.begin(); muon!=myMuons.selectedMuons.end();muon++)
  {  
     if(maxTrackProgress==11){continue;}
     myTracks.selectedTracks.clear(); 
     myTracks.selectedEndcapTracks.clear();     
     int trackProgress = myTracks.SelectTracks(iEvent, trackCollection_label, vtxHandle, *muon, transientTrackBuilder);
  //Pair our candidate track with a tagging muon to reconstruct a Z
     if(trackProgress==6){trackProgress = trackProgress+PairTagAndProbe(iEvent,iSetup,myTracks,*muon,myECAL, myHCAL, transientTrackBuilder,trackCollection_label, vtxHandle,reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label);}
     if(trackProgress>maxTrackProgress){maxTrackProgress=trackProgress;}
  }
  if(maxTrackProgress<11){return false;}
  Paired = true;
 
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedTracks.begin(); iTrack != myTracks.selectedTracks.end(); ++iTrack)
  {  
     if(deltaR((*iTrack)->eta(),(*iTrack)->phi(),selectedTrack->eta(),selectedTrack->phi())<0.3) continue;
     if(selectedMuon->charge()==(*iTrack)->charge()) continue;
     //if(fabs((*iTrack)->dxy())>0.05) continue;
     if(myTracks.PairTracks(iTrack, selectedMuon->globalTrack(), transientTrackBuilder))
     {        
//       double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label,0.3,vtxHandle,*iTrack);
//       double probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
       //if((probeTrackIso/(*iTrack)->pt())>0.2||probeEcalIso>30.) continue;
       Paired=false;
     }
  }
  
  
  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return false;
  info.nPassingTag=myMuons.selectedMuons.size();
  info.paired=Paired;
 
  if(Paired)
  {
    std::cout<< "Passing event" << std::endl;
    return true;
  }else return false;

}

bool MuPlusXFilter_HardBrem_RAW_RECO::MatchTrackToMuon(const edm::Event& iEvent,const reco::Track* selTrack, Muons myMuons, edm::EDGetToken m_recoMuonToken, MuPXEventInfo& info)
{
  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  bool matched = false;
  double mindPt = 1.;
  if(myMuons.selectedMuons.size()==0){return matched;}
  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); ++iMuon)
  {
     if(!iMuon->isGlobalMuon()) {continue;}
     double dR = deltaR(iMuon->eta(),iMuon->phi(),selTrack->eta(),selTrack->phi());
     double dPtOverPt =  std::abs((iMuon->pt()-selTrack->pt())/selTrack->pt());
     if(dR > 0.2) continue;
     if(dPtOverPt<mindPt)
     {
        mindPt=dPtOverPt;
	matched=true;
        info.matchedProbeGlobalMuonE = iMuon->energy();
    }
  }
  return matched;
}

bool MuPlusXFilter_HardBrem_RAW_RECO::MatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double* genInfo)
{
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(m_genParticleToken, genParticles);
   int nmuons=0;
   int ntotalmuons=0;
   double farMuE=0;
   double tagMuDr=1000;
   double farMuDr=0;
   double farMuEta=0; 
   double farMuProbeDr=-1;
   for (const auto& particle : *(genParticles.product()))
   {
      if (std::abs(particle.pdgId()) != 13)
         continue;
      ntotalmuons++;
      if (particle.mother()->pdgId() != 23)
         continue;
      nmuons++;
      double dR = deltaR(particle.eta(), particle.phi(), selTrack->eta(), selTrack->phi());
      double tagDr = deltaR(particle.eta(), particle.phi(), selectedMuon->eta(), selectedMuon->phi());
      if(tagDr<tagMuDr)
      {
         tagMuDr=tagDr;
      }
      if(tagDr>farMuDr)
      {
         farMuE=particle.pt();
         farMuDr=tagDr;
         farMuProbeDr=dR;
         farMuEta=particle.eta();
      }
      if(dR>0.2) continue;
      return true;
   }
   if(nmuons==2)
   {
      genInfo[0]=farMuE;
      genInfo[1]=farMuProbeDr;
      genInfo[2]=farMuEta;
      genInfo[3]=tagMuDr;
   }
   return false;
}

bool MuPlusXFilter_HardBrem_RAW_RECO::LooseMatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double* genInfo, std::vector<int> &nearbyIDs, std::vector<int> &highPtIds)
{
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(m_genParticleToken, genParticles);
   bool match=false;
   int nmuons=0;
   int ntotalmuons=0;
   int motherPDGId=0;
   int tagPDGId=0;
   double farMuE=0;
   double tagMuDr=1000;
   double farMuDr=-1;
   double farMuEta=0;
   double farMuProbeDr=-1;
   genInfo[7]=0;
   for (const auto& particle : *(genParticles.product()))
   {
      if (std::abs(particle.pdgId()) != 13)
      {
         double dR = deltaR(particle.eta(), particle.phi(), selTrack->eta(), selTrack->phi());
         if(dR<0.2)
         {
            if(particle.pt()>15){highPtIds.push_back(particle.pdgId());}
         }
         continue;
      }
   //   if (fabs(particle.mother()->pdgId()) != 13)
  //      continue;
      if (fabs(particle.mother()->pdgId())==23)
         nmuons++;
      ntotalmuons++;
      double dR = deltaR(particle.eta(), particle.phi(), selTrack->eta(), selTrack->phi());
      double tagDr = deltaR(particle.eta(), particle.phi(), selectedMuon->eta(), selectedMuon->phi());
      if(tagDr<tagMuDr)
      {
         tagMuDr=tagDr;
         const reco::Candidate* mother = particle.mother();
         while(fabs(mother->pdgId())==13)
         {
            mother=mother->mother();
         }
         tagPDGId=mother->pdgId();
      }
      if(dR<farMuDr||farMuDr<0)
      {
         farMuE=particle.pt();
         farMuDr=dR;
         farMuProbeDr=dR;
         farMuEta=particle.eta();
         const reco::Candidate* mother = particle.mother();
         while(fabs(mother->pdgId())==13)
         {
            mother=mother->mother();
         }
         motherPDGId=mother->pdgId();  
      }
      if(dR>0.2) continue;
      match=true;
   }
   //printf("Tag mo id is %i. Probe mo id is %i.\n",motherPDGId,tagPDGId);
   genInfo[0]=farMuE;
   genInfo[1]=farMuProbeDr;
   genInfo[2]=farMuEta;
   genInfo[3]=tagMuDr;
   genInfo[4]=motherPDGId;
   genInfo[5]=ntotalmuons;
   genInfo[6]=tagPDGId;
   return match;
}

int MuPlusXFilter_HardBrem_RAW_RECO::PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, HCAL myHCAL, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label)
{
  int cutpro = 0;
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
  {
     if(Tag->charge()==(*iTrack)->charge()) continue;
     if(cutpro<1){cutpro=1;}
     int pairingResult = myTracks.PairTracks(iTrack, Tag->globalTrack(), transientTrackBuilder);
     if(cutpro<(1+pairingResult)){cutpro=pairingResult+1;}
     if(pairingResult==3)
     {
       double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label,0.3,vtxHandle,*iTrack);
       double probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
       std::vector<double> HCALisolationInfo=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*iTrack));
       if(!m_signalRegion)
       {
          if((probeTrackIso/(*iTrack)->pt())<0.1&&probeEcalIso<30.) continue;
       }
       else
       {
          if((probeTrackIso/(*iTrack)->pt())>0.05) continue;
          if(probeEcalIso>10.) continue;
          if(HCALisolationInfo[0] < 30) continue;
       }
       if(cutpro<5){cutpro=5;}
       selectedTrack = (*iTrack);
       selectedMuon = Tag;
       pairVtxChi = myTracks.pairvertexchi;
       muTrackMass = myTracks.MuonTrackMass;
       reco::TransientTrack tk = transientTrackBuilder->build(*iTrack);
       TrajectoryStateClosestToPoint  traj = tk.trajectoryStateClosestToPoint(myTracks.probeVertPoint);
       probeTIP = traj.perigeeParameters().transverseImpactParameter();
       probeLIP = traj.perigeeParameters().longitudinalImpactParameter();
     }
  }
  return cutpro;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuPlusXFilter_HardBrem_RAW_RECO::beginStream(edm::StreamID)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuPlusXFilter_HardBrem_RAW_RECO::endStream() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuPlusXFilter_HardBrem_RAW_RECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuPlusXFilter_HardBrem_RAW_RECO);
