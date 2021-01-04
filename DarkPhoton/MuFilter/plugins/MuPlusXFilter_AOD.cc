// -*- C++ -*-
//
// Package:    MuMu/MuPlusXFilter_AOD
// Class:      MuPlusXFilter_AOD
// 
/**\class MuPlusXFilter_AOD MuPlusXFilter_AOD_AOD.cc MuMu/MuPlusXFilter_AOD/plugins/MuPlusXFilter_AOD_AOD.cc

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
#include "DarkPhoton/MuFilter/interface/eventHistos.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"

// HCAL Stuff
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"

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

class MuPlusXFilter_AOD : public edm::stream::EDFilter<>  {
   public:
      explicit MuPlusXFilter_AOD(const edm::ParameterSet&);
      ~MuPlusXFilter_AOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<reco::Track>::const_iterator&);

    double TrackIsolation(const edm::Event&, edm::EDGetTokenT<std::vector<reco::Track>>,double, edm::Handle<reco::VertexCollection>, std::vector<reco::Track>::const_iterator&);
    double HCALIsolation(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>, const reco::TransientTrack);
    void GetTransientProjectedCells(const CaloSubdetectorGeometry*, HcalDetId*, reco::TransientTrack);

    eventHistos m_allEvents;
    eventHistos m_passingEvents;
    edm::EDGetToken m_recoMuonToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;
    edm::EDGetTokenT<edm::TriggerResults> trigResults_Label;
    edm::EDGetToken m_genParticleToken;
    const reco::Track* selTrack;
    const reco::Muon* selMuon;
    double selVtxChi;
    double selMuonTrackMass;
    double selTrackIso;
    double selHcalIso;
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
MuPlusXFilter_AOD::MuPlusXFilter_AOD(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  trigResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResultsTag"))),
  selVtxChi(0),
  selMuonTrackMass(0),
  selTrackIso(1000),
  selHcalIso(1000),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",true))
{
  edm::Service<TFileService> fs;
  m_allEvents.book(fs->mkdir("allEvents"));
  m_passingEvents.book(fs->mkdir("passingEvents"));
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
  }

}


MuPlusXFilter_AOD::~MuPlusXFilter_AOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
bool
MuPlusXFilter_AOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  m_allEvents.ResetCutFlow();
  m_passingEvents.ResetCutFlow();
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trigResults_Label, trigResults);
  if(!m_isMC)
  {
     const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
     std::string pathName="HLT_IsoMu24_v";
     unsigned int trigIndex=0;
     bool foundtrig = false;
     for(unsigned int itrig=0; itrig!=trigResults->size(); ++itrig)
     {
        if(trigNames.triggerName(itrig).find(pathName)!=std::string::npos)
        {
           trigIndex=itrig;
           foundtrig = true;
        }
     }
     if(!foundtrig){return false;}
     bool passTrig=trigResults->accept(trigIndex);
     if(!passTrig) return false;
  }
  m_allEvents.IncCutFlow();
  m_passingEvents.IncCutFlow();
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
  m_allEvents.IncCutFlow();
  m_passingEvents.IncCutFlow();

  int anyMuonPass = 0;
  int nMuonTrackCand = 0;
  int nTracksNoMuon = 0;
  int nTotalMuonTrackCand = 0;
  float MuonTrackMass = 0.;

  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);

  // this wraps tracks with additional methods that are used in vertex-calculation
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  //Looping over the reconstructed Muons
  
  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {
    if(!(iMuon->passed(reco::Muon::CutBasedIdTight))) continue;
    if(!(iMuon->passed(reco::Muon::PFIsoTight))) continue;
    if(!(iMuon->passed(reco::Muon::TkIsoTight))) continue;
    if (iMuon->pt() < 26 || fabs(iMuon->eta()) > 2.4) continue;
    anyMuonPass++;

    //Looping over tracks
    for(std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) {

        if (iMuon->charge()==iTrack->charge()) continue;

        if (!iTrack->quality(Track::highPurity)) continue;

        if (fabs(iTrack->eta()) > 2.5 || fabs(iTrack->eta()) < 1.3) continue;//Track required to be in the endcap

        if (iTrack->pt()<15) continue;
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
        double vtxChi = 0;
        // some poor fits will simply fail to find a common vertex
        if (fittedVertex.isValid()  &&  fittedVertex.totalChiSquared() >= 0.  &&  fittedVertex.degreesOfFreedom() > 0)        {
           // others we can exclude by their poor fit
           vtxChi = fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom(); 

           if (vtxChi < 3.) 
           {
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
           }
           else{continue;}
        }
        else{continue;}

/*	TLorentzVector* TrackVector;
        TLorentzVector* MuonVector;
	TrackVector->SetPtEtaPhiM(iTrack->pt(),iTrack->eta(),iTrack->phi(),iTrack->m());*/
//        if ((iMuon->p4() + iTrack->p4()).mass() < 85 || (iMuon->p4() + iTrack->p4()).mass() > 105) continue;
	if (MuonTrackMass < 50 || MuonTrackMass > 150) continue;
        nTotalMuonTrackCand++;
        double trackIso = TrackIsolation(iEvent,trackCollection_label, 0.3, vertices, iTrack);
        double hcalIso = HCALIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*iTrack));
       
        if(nMuonTrackCand==0)
        {
           selTrack = &(*iTrack);
           selMuon = &(*iMuon);
           selVtxChi = vtxChi;
           selMuonTrackMass = MuonTrackMass;
           selTrackIso = trackIso;
           selHcalIso = hcalIso;
        }
        if(trackIso<0.15&&hcalIso<3) continue;

        if (m_isMC){
	  if (!isTrackMatchedToMuon(iEvent, iTrack)) continue;//Matching the track to a GEN muon
	}
        selTrack = &(*iTrack);
        selMuon = &(*iMuon);
        selVtxChi = vtxChi;
        selMuonTrackMass = MuonTrackMass;
        selTrackIso = trackIso;
        selHcalIso = hcalIso;
        nMuonTrackCand++;
//	MuonTrackMass = (iMuon->p4() + iTrack->p4()).mass();

	nTracksNoMuon=0.;
        for(std::vector<reco::Muon>::const_iterator iMuon2 = recoMuons->begin(); iMuon2 != recoMuons->end(); iMuon2++)
        {

           if (!(iMuon2->isPFMuon() && (iMuon2->isGlobalMuon() || iMuon2->isTrackerMuon()))) continue;
           if (fabs(iMuon2->eta()) > 2.8 || fabs(iMuon2->eta()) < 1.653) continue;
	   if(deltaR(iTrack->eta(), iTrack->phi(), iMuon2->eta(), iMuon2->phi()) < 0.1){
             nTracksNoMuon++;

	   }
	}

    }
  }
  if(anyMuonPass>0)
  {
    m_allEvents.IncCutFlow();
    m_passingEvents.IncCutFlow();
    if(nTotalMuonTrackCand>0)
    {
       m_allEvents.IncCutFlow();
       m_passingEvents.IncCutFlow();
       if(nMuonTrackCand>0)
       {
         m_allEvents.IncCutFlow();          
         m_passingEvents.IncCutFlow();      
       }
    }
  }
  m_allEvents.m_eventCount->Fill(1);
  m_allEvents.m_NPassingTag->Fill(anyMuonPass);

  if(nTotalMuonTrackCand==0) return false;  
  m_allEvents.m_MuonTrackMass->Fill(selMuonTrackMass);
  m_allEvents.m_ProbeEta->Fill(selTrack->eta());
  m_allEvents.m_ProbePt->Fill(selTrack->pt());
  m_allEvents.m_ProbePhi->Fill(selTrack->phi());
  m_allEvents.m_NPassingProbe->Fill(nMuonTrackCand);
  m_allEvents.m_TagEta->Fill(selMuon->eta());
  m_allEvents.m_TagPt->Fill(selMuon->pt());
  m_allEvents.m_TagPhi->Fill(selMuon->phi());
  m_allEvents.m_TagProbeVtxChi->Fill(selVtxChi);
  m_allEvents.m_ProbeTrackIso->Fill(selTrackIso);
  m_allEvents.m_ProbeHcalIso->Fill(selHcalIso);
  m_allEvents.m_ProbeCombinedIso->Fill(selTrackIso,selHcalIso);
 
  if (nMuonTrackCand > 0)
  {
    std::cout <<"PASSES"<<std::endl;
    m_passingEvents.m_eventCount->Fill(1);
    m_passingEvents.m_MuonTrackMass->Fill(selMuonTrackMass);
    m_passingEvents.m_ProbeEta->Fill(selTrack->eta());
    m_passingEvents.m_ProbePt->Fill(selTrack->pt());
    m_passingEvents.m_ProbePhi->Fill(selTrack->phi());
    m_passingEvents.m_NPassingProbe->Fill(nMuonTrackCand);
    m_passingEvents.m_TagEta->Fill(selMuon->eta());
    m_passingEvents.m_TagPt->Fill(selMuon->pt());
    m_passingEvents.m_TagPhi->Fill(selMuon->phi());
    m_passingEvents.m_NPassingTag->Fill(anyMuonPass);
    m_passingEvents.m_TagProbeVtxChi->Fill(selVtxChi);
    m_passingEvents.m_ProbeTrackIso->Fill(selTrackIso);
    m_passingEvents.m_ProbeHcalIso->Fill(selHcalIso);
    m_passingEvents.m_ProbeCombinedIso->Fill(selTrackIso,selHcalIso);
    return true;
  }else
  {
    std::cout << "FAILS"<<std::endl;
    return false;
  }		
}

bool MuPlusXFilter_AOD::isTrackMatchedToMuon(const edm::Event& iEvent, std::vector<reco::Track>::const_iterator& Track){

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

double MuPlusXFilter_AOD::TrackIsolation(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, double conesize, edm::Handle<reco::VertexCollection> vtxHandle, std::vector<reco::Track>::const_iterator& iTrack)
{
   bool foundtrack = false;
   unsigned int vtxindex = 0;
   unsigned int trackindex = 0;
   double Isolation  =0;
   for(unsigned int i=0;i<vtxHandle->size();i++)
   {
      reco::VertexRef vtx(vtxHandle, i);
      if(!vtx->isValid()){continue;}
      for(unsigned int j=0; j<vtx->tracksSize();j++)
      {
         if(fabs(vtx->trackRefAt(j)->pt()-iTrack->pt())<0.01)
         {
             vtxindex = i;
             trackindex = j;
             foundtrack = true;
         }    
      }
   }
   
   if(!foundtrack){return 100;}

   reco::VertexRef primaryVtx(vtxHandle,vtxindex);
  
   for(unsigned int i=0;i<primaryVtx->tracksSize();i++)
   {
      if(i==trackindex) continue;
      reco::TrackBaseRef secondarytrack = primaryVtx->trackRefAt(i);
      if(deltaR(iTrack->eta(), iTrack->phi(), secondarytrack->eta(),secondarytrack->phi())>conesize) continue;
      Isolation += secondarytrack->pt();
   }

   return Isolation/iTrack->pt();
}

double MuPlusXFilter_AOD::HCALIsolation(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, const reco::TransientTrack track)
{
   edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
   iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
   double MatchedEnergy = 0;

   edm::ESHandle<CaloGeometry> TheCALOGeometry;
   iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);

   const CaloGeometry* caloGeom = TheCALOGeometry.product();
   const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);

   HcalDetId MatchedCells[20];

   GetTransientProjectedCells(HEGeom,MatchedCells,track);

   if(!hcalRecHits.isValid()) return 1000;
   const HBHERecHitCollection *hbhe = hcalRecHits.product();
   for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
   {
      HcalDetId id(hbherechit->detid());
      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());

      HcalDetId *match = std::find(std::begin(MatchedCells), std::end(MatchedCells), id);
      if(match!=std::end(MatchedCells)) MatchedEnergy+=hbherechit->energy();
   }
 
   return MatchedEnergy;
}

void MuPlusXFilter_AOD::GetTransientProjectedCells(const CaloSubdetectorGeometry* HEGeom, HcalDetId *MatchedCells, reco::TransientTrack muTrack)
{
   double start, step, end;
   start = 320;
   step = 5;
   end = 530;
   int j=0;
   HcalDetId lastClosestCell;
   for(int i=start; i<end; i+=step)
   {
      double testPointPerp = fabs(i*tan(muTrack.track().theta()));
      double testPointX = testPointPerp*cos(muTrack.track().phi());
      double testPointY = testPointPerp*sin(muTrack.track().phi());
      double testPointZ;
      if(muTrack.track().eta()>9){testPointZ = i;}
      else{testPointZ=-i;}

      GlobalPoint testGlobalPoint = GlobalPoint(testPointX,testPointY,testPointZ);
      TrajectoryStateClosestToPoint traj = muTrack.trajectoryStateClosestToPoint(testGlobalPoint);
      HcalDetId testClosestCell = (HcalDetId)HEGeom->getClosestCell(traj.position());
      if(testClosestCell!=lastClosestCell)
      {
         lastClosestCell=testClosestCell;
         MatchedCells[j] = testClosestCell;
         j = j+1;
      }
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuPlusXFilter_AOD::beginStream(edm::StreamID)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuPlusXFilter_AOD::endStream() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuPlusXFilter_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuPlusXFilter_AOD);
