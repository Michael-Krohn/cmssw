// -*- C++ -*-
//
// Package:    MuMu/baseEventFilter
// Class:      baseEventFilter
// 
/**\class baseEventFilter baseEventFilter_AOD.cc MuMu/baseEventFilter/plugins/baseEventFilter_AOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Revering
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
#include "Math/VectorUtil.h"
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
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
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

//Jets
#include "DataFormats/JetReco/interface/PFJetCollection.h"

// Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class baseEventFilter : public edm::stream::EDFilter<>  {
   public:
      explicit baseEventFilter(const edm::ParameterSet&);
      ~baseEventFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<reco::Track>::const_iterator&);

    double TrackIsolation(const edm::Event&, edm::EDGetTokenT<std::vector<reco::Track>>,double, edm::Handle<reco::VertexCollection>, std::vector<reco::Track>::const_iterator&);
    double ECALIsolation(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<EcalRecHitCollection>, edm::EDGetTokenT<EcalRecHitCollection>, const reco::TransientTrack);

    eventHistos m_allEvents;
    eventHistos m_passingEvents;
    edm::EDGetToken m_recoMuonToken;
    edm::EDGetToken m_theSTAMuonLabel;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label;
    edm::EDGetTokenT<edm::TriggerResults> trigResults_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetToken m_genInfoToken;
    edm::EDGetToken m_pfJetCollection_label;
    const reco::Track* selTrack;
    const reco::Muon* selMuon;
    double selVtxChi;
    double selMuonTrackMass;
    double selTrackIso;
    double selMinDr;
    double selMinDEoverE;
    double selEcalIso;
    double weight_;
    bool sameSign;
    bool totalRegion;
    int eventNum;
    int lumi;
    int run;
    int category;
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
baseEventFilter::baseEventFilter(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  m_theSTAMuonLabel(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("StandaloneTracks"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  reducedEndcapRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
  reducedBarrelRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
  trigResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResultsTag"))),
  m_pfJetCollection_label(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PFJets"))),
  selVtxChi(0),
  selMuonTrackMass(0),
  selTrackIso(-1),
  selEcalIso(-1),
  eventNum(-1),
  lumi(-1),
  run(-1),
  category(-1)
{
  edm::Service<TFileService> fs;
  m_allEvents.book(fs->mkdir("allEvents"));
  m_passingEvents.book(fs->mkdir("passingEvents"));
}


baseEventFilter::~baseEventFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
bool
baseEventFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
 
  totalRegion=false;
  sameSign=false; 
  weight_ = 1; 
  int njets=0;
  run=-1;
  eventNum=-1;
  lumi=-1;
  category=-1;
  m_allEvents.m_eventCount->Fill(0.5, weight_);
  m_allEvents.ResetCutFlow();
  m_passingEvents.ResetCutFlow();
  edm::Handle<reco::TrackCollection> staTracks;
  iEvent.getByToken(m_theSTAMuonLabel, staTracks);
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trigResults_Label, trigResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
  std::string pathName="HLT_IsoMu24_v";
  std::string pathName2="HLT_IsoMu27_v";
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
  m_allEvents.IncCutFlow();
  if(!foundtrig){return false;}
  m_allEvents.IncCutFlow();
  bool passTrig=trigResults->accept(trigIndex);
  if(!passTrig)
  {
     for(unsigned int itrig=0; itrig!=trigResults->size(); ++itrig)
     {
        if(trigNames.triggerName(itrig).find(pathName2)!=std::string::npos)
        {
           trigIndex=itrig;
        }
     }
  }
  passTrig=trigResults->accept(trigIndex);
  if(!passTrig){return false;}
  m_allEvents.IncCutFlow();
  m_passingEvents.IncCutFlow();
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(primaryVertices_Label, vertices);

  m_allEvents.IncCutFlow();
  m_passingEvents.IncCutFlow();

  int anyMuonPass = 0;
  int nMuonTrackCand = 0;
  int nTracksNoMuon = 0;
  int nTotalMuonTrackCand = 0;
  selMinDEoverE=10.;
  float MuonTrackMass = 0.;

  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);

  // this wraps tracks with additional methods that are used in vertex-calculation
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  //Looping over the reconstructed Muons
  int maxProg=0;
  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {
    if(!(iMuon->isPFMuon() && iMuon->isGlobalMuon())) continue;
    if(!(iMuon->passed(reco::Muon::CutBasedIdTight))) continue;
    if(maxProg==0){maxProg=1;}
    if(!(iMuon->passed(reco::Muon::PFIsoTight))) continue;
    if(maxProg==1){maxProg=2;}
    if (iMuon->pt() < 26 || fabs(iMuon->eta()) > 2.4) continue;
    if(maxProg==2){maxProg=3;}
    anyMuonPass++;

    //Looping over tracks
    for(std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) {
        bool sign=false;
        if (iMuon->charge()==iTrack->charge()){sign=true;}

        if (!iTrack->quality(reco::Track::qualityByName("highPurity"))) continue;
        if(maxProg==3){maxProg=4;}
        if (fabs(iTrack->eta()) > 2.4 || fabs(iTrack->eta()) < 1.4) continue;//Track required to be in the endcap
        if(maxProg==4){maxProg=5;}
        //if (deltaR(iTrack->eta(),iTrack->phi(),iMuon->eta(),iMuon->phi())<0.2) continue;
        if(maxProg==5){maxProg=6;}
        if (iTrack->pt()<20) continue;
        if(maxProg==6){maxProg=7;}
        bool foundtrack=false; 
        GlobalPoint tkVtx;
        for(unsigned int i=0;i<vertices->size();i++)
        {
           reco::VertexRef vtx(vertices, i);
           if(!vtx->isValid()){continue;}
           for(unsigned int j=0; j<vtx->tracksSize();j++)
           {
              if(vtx->trackRefAt(j)->pt()==iTrack->pt()&&vtx->trackRefAt(j)->eta()==iTrack->eta()&&vtx->trackRefAt(j)->phi()==iTrack->phi())
              {
                 foundtrack=true;
                 GlobalPoint vert(vtx->x(),vtx->y(),vtx->z());
                 tkVtx=vert;
              }
           } 
        }
        if(!foundtrack) continue;
        if(maxProg==7){maxProg=8;}
        reco::TransientTrack tk = transientTrackBuilder->build(*iTrack);
        TrajectoryStateClosestToPoint traj = tk.trajectoryStateClosestToPoint(tkVtx);
        double transDCA = traj.perigeeParameters().transverseImpactParameter();
        double longDCA = traj.perigeeParameters().longitudinalImpactParameter();
        if(longDCA>0.05) continue; 
        if(transDCA>0.005) continue;
        if(maxProg==8){maxProg=9;}
        // make a pair of TransientTracks to feed the vertexer
        std::vector<reco::TransientTrack> tracksToVertex;
        tracksToVertex.push_back(transientTrackBuilder->build(*iTrack));
	if(iMuon->isGlobalMuon()){
          tracksToVertex.push_back(transientTrackBuilder->build(iMuon->globalTrack()));
	}
	else{continue;}
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
        if(maxProg==9){maxProg=10;}
        if(sign)
        {
          if(MuonTrackMass<50 || MuonTrackMass>150) continue;
        }
        else{
	  if (MuonTrackMass < 50 || MuonTrackMass > 150) continue;
          //if (MuonTrackMass > 76 && MuonTrackMass <106) continue;
        }
        if(maxProg==10){maxProg=11;}
        double trackIso = TrackIsolation(iEvent,trackCollection_label, 0.3, vertices, iTrack);
        double ecalIso = ECALIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));      
        if(nMuonTrackCand==0)
        {
           selTrack = &(*iTrack);
           selMuon = &(*iMuon);
           selVtxChi = vtxChi;
           selMuonTrackMass = MuonTrackMass;
           selTrackIso = trackIso;
           selEcalIso = ecalIso;
           nTotalMuonTrackCand++;
           sameSign=sign;
        }
        //if(trackIso<0.2&&ecalIso<30) continue;
        if(trackIso>0.05||ecalIso>10) continue;
        if(maxProg==11){maxProg=12;}

        double staMinDr=10;
        double staMinDEoverE=-10;
        if(staTracks->size()>0)
        {
           for(reco::TrackCollection::const_iterator staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack)
           {
              reco::TransientTrack track = transientTrackBuilder->build(*staTrack);
              double dR = deltaR(track.impactPointTSCP().momentum().eta(),track.impactPointTSCP().momentum().phi(),(*iTrack).eta(),(*iTrack).phi());
              double staDE=(std::sqrt(track.impactPointTSCP().momentum().mag2())-(*iTrack).p())/(*iTrack).p();
              if(dR<staMinDr)
              {
                 staMinDr=dR;
              }
              if(dR<1&&staDE>staMinDEoverE){staMinDEoverE=staDE;}
           }
        }
        if(staMinDr>1.0){totalRegion=true;}
	if(staMinDr<1.0)
	{
	   if(staMinDEoverE<selMinDEoverE)
	   {
              selMinDr=staMinDr;
              selMinDEoverE=staMinDEoverE;
              selTrack = &(*iTrack);
              selMuon = &(*iMuon);
              selVtxChi = vtxChi;
              selMuonTrackMass = MuonTrackMass;
              selTrackIso = trackIso;
              selEcalIso = ecalIso;
	   }
	}
        sameSign=sign;
        nMuonTrackCand++;
//	MuonTrackMass = (iMuon->p4() + iTrack->p4()).mass();
    }
  }
  for(int prog=0; prog<maxProg; prog++)
  {
    m_allEvents.IncCutFlow();
  }
  edm::Handle<reco::PFJetCollection> pfjetHandle;
  iEvent.getByToken(m_pfJetCollection_label, pfjetHandle);

  if(pfjetHandle.isValid())
  {
     const reco::PFJetCollection pfjets = *(pfjetHandle.product());
     for(reco::PFJetCollection::const_iterator it=pfjets.begin();it!=pfjets.end();++it)
     {
        double NHF, NEMF, CHF, MUF, CEMF;  
        int CHM, NumConst;
        NHF = it->neutralHadronEnergyFraction();
        NEMF = it->neutralEmEnergyFraction();
        CHF = it->chargedHadronEnergyFraction();
        MUF = it->muonEnergyFraction();
        CEMF = it->chargedEmEnergyFraction();
        NumConst = it->chargedMultiplicity() + it->neutralMultiplicity();
        CHM = it->chargedMultiplicity();
        bool JetId;
        if(fabs(it->eta())>2.4){continue;}
        JetId = (CEMF < 0.8 && CHM > 0 && CHF > 0 && NumConst > 1 && NEMF < 0.9 && MUF < 0.8 && NHF < 0.9);   
        if(!JetId){continue;}
        if(it->pt()<30){continue;}
        njets++;
     }
  }

/*  if(anyMuonPass>0)
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
  }*/
  m_allEvents.m_NPassingTag->Fill(anyMuonPass,weight_);
  if(nTotalMuonTrackCand==0) return false;  
  m_allEvents.m_MuonTrackMass->Fill(selMuonTrackMass,weight_);
  m_allEvents.m_ProbeEta->Fill(selTrack->eta(),weight_);
  m_allEvents.m_ProbePt->Fill(selTrack->pt(),weight_);
  m_allEvents.m_staMinDr->Fill(selMinDr,weight_);
  m_allEvents.m_staMinDEoverE->Fill(selMinDEoverE,weight_);
  m_allEvents.m_ProbePhi->Fill(selTrack->phi(),weight_);
  m_allEvents.m_ProbeEtaPhi->Fill(selTrack->eta(),selTrack->phi(),weight_);
  m_allEvents.m_NPassingProbe->Fill(nMuonTrackCand,weight_);
  m_allEvents.m_TagEta->Fill(selMuon->eta(),weight_);
  m_allEvents.m_TagPt->Fill(selMuon->pt(),weight_);
  m_allEvents.m_TagPhi->Fill(selMuon->phi(),weight_);
  m_allEvents.m_TagProbeVtxChi->Fill(selVtxChi,weight_);
  m_allEvents.m_ProbeTrackIso->Fill(selTrackIso,weight_);
  m_allEvents.m_ProbeEcalIso->Fill(selEcalIso,weight_);
  m_allEvents.m_NJets->Fill(njets);
  bool passes = false; 
  if (nMuonTrackCand > 0)
  {
    eventNum=iEvent.id().event();
    lumi=iEvent.id().luminosityBlock();
    run=iEvent.id().run();
    if(!sameSign)
    {
      m_passingEvents.m_NJets->Fill(njets);
      if (selMuonTrackMass > 76 && selMuonTrackMass <106)
      {
	 if(selMinDEoverE<-0.5){category=7;}
      }
      else
      {
         if(njets<=1)
         {
            if(totalRegion){category=1;}
            else{category=2;}
         }
         else if(njets>1)
         {
            if(totalRegion){category=3;}
            else{category=4;}
         }
      }
    }
    else
    {
      if(totalRegion){category=5;}
      else{category=6;}
    }
    if(category>0){passes=true;}
  }
  if(passes)
  {
    std::cout <<"PASSES"<<std::endl;
    m_passingEvents.m_eventCount->Fill(0.5,weight_);
    m_passingEvents.m_staMinDr->Fill(selMinDr,weight_);
    m_passingEvents.m_staMinDEoverE->Fill(selMinDEoverE,weight_);
    m_passingEvents.m_eventWeight->Fill(weight_);
    m_passingEvents.m_MuonTrackMass->Fill(selMuonTrackMass,weight_);
    m_passingEvents.m_ProbeEta->Fill(selTrack->eta(),weight_);
    m_passingEvents.m_ProbePt->Fill(selTrack->pt(),weight_);
    m_passingEvents.m_ProbePhi->Fill(selTrack->phi(),weight_);
    m_passingEvents.m_ProbeEtaPhi->Fill(selTrack->eta(),selTrack->phi(),weight_);
    m_passingEvents.m_NPassingProbe->Fill(nMuonTrackCand,weight_);
    m_passingEvents.m_TagEta->Fill(selMuon->eta(),weight_);
    m_passingEvents.m_TagPt->Fill(selMuon->pt(),weight_);
    m_passingEvents.m_TagPhi->Fill(selMuon->phi(),weight_);
    m_passingEvents.m_NPassingTag->Fill(anyMuonPass,weight_);
    m_passingEvents.m_TagProbeVtxChi->Fill(selVtxChi,weight_);
    m_passingEvents.m_ProbeTrackIso->Fill(selTrackIso,weight_);
    m_passingEvents.m_ProbeEcalIso->Fill(selEcalIso,weight_);
    m_passingEvents.FillTree(run,lumi,eventNum,category);
    eventNum=iEvent.id().event();
    lumi=iEvent.id().luminosityBlock();
    run=iEvent.id().run();
    return true;
  }else
  {
    std::cout << "FAILS"<<std::endl;
    return false;
  }		
}

bool baseEventFilter::isTrackMatchedToMuon(const edm::Event& iEvent, std::vector<reco::Track>::const_iterator& Track){

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

double baseEventFilter::TrackIsolation(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, double conesize, edm::Handle<reco::VertexCollection> vtxHandle, std::vector<reco::Track>::const_iterator& iTrack)
{
   bool foundtrack = false;
   unsigned int vtxindex = -1;
   unsigned int trackindex = -1;
   double Isolation  =0;
   for(unsigned int i=0;i<vtxHandle->size();i++)
   {
      reco::VertexRef vtx(vtxHandle, i);
      if(!vtx->isValid()){continue;}
      for(unsigned int j=0; j<vtx->tracksSize();j++)
      {
         if(fabs(vtx->trackRefAt(j)->pt()==iTrack->pt()))
         {
             vtxindex = i;
             trackindex = j;
             foundtrack = true;
             GlobalPoint vert(vtx->x(), vtx->y(), vtx->z());
         }    
      }
   }
   
   if(!foundtrack){return -1;}

   reco::VertexRef primaryVtx(vtxHandle,vtxindex);
  
   for(unsigned int i=0;i<primaryVtx->tracksSize();i++)
   {
      if(i==trackindex) continue;
      reco::TrackBaseRef secondarytrack = primaryVtx->trackRefAt(i);
      if(deltaR(iTrack->eta(), iTrack->phi(), secondarytrack->eta(),secondarytrack->phi())>conesize || deltaR(iTrack->eta(), iTrack->phi(), secondarytrack->eta(),secondarytrack->phi())<0.01) continue;
      Isolation += secondarytrack->pt();
   }

   return Isolation/iTrack->pt();
}

double baseEventFilter::ECALIsolation(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label, const reco::TransientTrack track)
{
   edm::Handle<EcalRecHitCollection> rechitsEE;
   iEvent.getByToken(reducedEndcapRecHitCollection_Label, rechitsEE);

   edm::Handle<EcalRecHitCollection> rechitsEB;
   iEvent.getByToken(reducedBarrelRecHitCollection_Label, rechitsEB);
 
   edm::ESHandle<CaloGeometry> TheCALOGeometry;
   iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);

   const CaloGeometry* caloGeom = TheCALOGeometry.product();
   TrajectoryStateClosestToPoint t0 = track.impactPointTSCP();
   double eDR = 0;

   for(EcalRecHitCollection::const_iterator hit = rechitsEE->begin(); hit!= rechitsEE->end(); hit++)
   {
      const DetId id = (*hit).detid();
      const GlobalPoint hitPos = caloGeom->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
      //Check if hit and track trajectory ar in the same endcap (transient track projects both ways)
      if ((hitPos.eta() * t0.momentum().eta()) < 0) {continue;}
      TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
      math::XYZVector idPositionRoot(hitPos.x(),hitPos.y(),hitPos.z());
      math::XYZVector trajRoot(traj.position().x(),traj.position().y(),traj.position().z());
      if(ROOT::Math::VectorUtil::DeltaR(idPositionRoot,trajRoot)<0.4&&(*hit).energy()>0.3)
      {
         eDR+= (*hit).energy();
      }            
   } 
   for (EcalRecHitCollection::const_iterator hit = rechitsEB->begin(); hit != rechitsEB->end(); hit++) {
      const DetId id = (*hit).detid();
      const GlobalPoint hitPos = caloGeom->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
      if ((hitPos.eta() * t0.momentum().eta()) < 0) {continue;}
      TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
      math::XYZVector idPositionRoot(hitPos.x(), hitPos.y(), hitPos.z());
      math::XYZVector trajRoot(traj.position().x(), traj.position().y(), traj.position().z());      if (ROOT::Math::VectorUtil::DeltaR(idPositionRoot, trajRoot) < 0.4 && (*hit).energy() > 0.3)
      {
         eDR += (*hit).energy();
      }
   }

   return eDR;
}

// ------------ method called once each job just before starting event loop  ------------
void 
baseEventFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
baseEventFilter::endStream() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
baseEventFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(baseEventFilter);
