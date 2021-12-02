// -*- C++ -*-
//
// Package:    DarkPhoton/PileupAnalyzer
// Class:      PileupAnalyzer
// 
/**\class PileupAnalyzer PileupAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/PileupAnalyzer.cc

 Description: Analyzer to study HCAL behavior of isolated charged pion and kaon tracks.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Revering
//         Created:  Friday, 29 Oct 2021 10:30:42 GMT
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

//Jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h" 

//Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DarkPhoton/MuAnalyzer/interface/PileupHistograms.h"
#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"
#include "DarkPhoton/MuAnalyzer/interface/ECAL.h"
#include "DarkPhoton/MuAnalyzer/interface/Jets.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PileupAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PileupAnalyzer(const edm::ParameterSet&);
      ~PileupAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, PileupHistograms);
    bool MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, SigEventInfo info);
    bool MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, SigEventInfo info, double FmuE, double minCSCdr);
    int PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetToken m_simTracksToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label;
    edm::EDGetToken m_dbremweighttoken;
    edm::EDGetToken m_dbremInitialEToken;
    edm::EDGetToken m_dbremFinalEToken;
    edm::EDGetToken m_dbremAngleToken;
    edm::EDGetToken m_dbremXToken;
    edm::EDGetToken m_dbremYToken;
    edm::EDGetToken m_dbremZToken;
    edm::EDGetToken m_trigResultsToken;
    edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
    std::vector<std::string> m_muonPathsToPass;
    const reco::Track* selectedTrack;
    const reco::Muon* selectedMuon;
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label;
    edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label;
    edm::EDGetToken m_caloJet_label;
    edm::EDGetToken m_theSTAMuonLabel;
    edm::EDGetToken m_pfJetCollection_label;
    bool m_isMC;
    bool m_isSig;
    bool m_hasDpho;
    double m_scaleFactor;
    bool m_runRandomTrackEfficiency;
    double standaloneE;
    double weight_;
    double pairVtxChi;
    double diMuonMass;

    PileupHistograms allEvents;
    PileupHistograms preSelectedEvents;
    PileupHistograms failSig;
    PileupHistograms passTight;
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
PileupAnalyzer::PileupAnalyzer(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  reducedEndcapRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
  reducedBarrelRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
  m_caloJet_label(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("CaloJetSource"))),
  m_theSTAMuonLabel(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("StandAloneTracks"))),
  m_pfJetCollection_label(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PFJets"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",true)),
  m_isSig (iConfig.getUntrackedParameter<bool>("isSig",false)),
  m_hasDpho (iConfig.getUntrackedParameter<bool>("hasDpho",false)),
  m_scaleFactor (iConfig.getUntrackedParameter<double>("scaleFactor",2000)),
  m_runRandomTrackEfficiency (iConfig.getUntrackedParameter<bool>("runRandomTrackEfficiency",false))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
    m_simTracksToken = consumes<edm::SimTrackContainer> (iConfig.getParameter<edm::InputTag>("g4SimHits"));
    m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    if(m_hasDpho)
    {
       m_dbremweighttoken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremweighttag"));
       m_dbremInitialEToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremInitialETag"));
       m_dbremFinalEToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremFinalETag"));
       m_dbremAngleToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremAngleTag"));
       m_dbremXToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremXTag"));
       m_dbremYToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremYTag"));
       m_dbremZToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremZTag"));
    }
  }
  else
  {
     m_trigResultsToken = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("trigResults"));
     m_muonPathsToPass   = iConfig.getParameter<std::vector<std::string> >("muonPathsToPass");
  }
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  allEvents.book(fs->mkdir("allEvents"));
  preSelectedEvents.book(fs->mkdir("preSelected"));
  failSig.book(fs->mkdir("failSig"));
  passTight.book(fs->mkdir("passTight"));
}


PileupAnalyzer::~PileupAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void PileupAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  Muons myMuons;
  Jets myJets;
  EventInfo myEventInfo;
  Tracks myTracks;
  CSC myCSCs;
  HCAL myHCAL;
  ECAL myECAL;
  SigEventInfo info;
  allEvents.ResetCutFlow();
  edm::Handle<SimTrackContainer> simtracks;
  bool genAccept = false;
  weight_ = 1.;
  info.eventWeight=weight_;
  edm::Handle<reco::TrackCollection> staTracks;
  iEvent.getByToken(m_theSTAMuonLabel, staTracks);
  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);
  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  int trackProgress = myTracks.SelectTracks(iEvent, trackCollection_label, vtxHandle);
  for(int i=0;i<trackProgress;i++){allEvents.IncCutFlow(info.eventWeight);}
  if(trackProgress<4){return;}
  //Find gen k or pi  
  double particleEta, particlePhi;
  int cutpro=0;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);
  double minMuDr=-1;
  /*for (const auto& particle : *(genParticles.product()))
  {
     if(std::abs(particle.pdgId())!=211&&std::abs(particle.pdgId())!=321) continue;
     if(cutpro<1) cutpro=1;
     if(std::abs(particle.eta())>2.4||std::abs(particle.eta())<1.6) continue;
     if(cutpro<2) cutpro=2;
     if(particle.pt()<15) continue;
     cutpro=3;
     particleEta = particle.eta();
     particlePhi = particle.phi();
     info.particlePDG=particle.pdgId();
     genAccept=true;
  }
  for (const auto& particle : *(genParticles.product()))
  {
     if(std::abs(particle.pdgId())!=13) continue;
     if(particle.pt()<15) continue;
     double dR=deltaR(particleEta,particlePhi,particle.eta(),particle.phi());
     if(minMuDr<0||dR<minMuDr){minMuDr=dR;} 
  }
  info.minMuDr=minMuDr;
  
  for(int i=0;i<cutpro;i++){allEvents.IncCutFlow(info.eventWeight);}
  if(!genAccept){return;}*/
  if(info.bremR==0){info.bremR=525;}
  //Find deflected muon or the initial muon if there was no deflection

  bool TrackMatched = false;
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);  
  double minPartDr= -1;
  int minPartPdg=0;
  cutpro=0;
  //Pair our dark brem candidate track with a tagging muon to reconstruct a Z
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
  {
     bool muClose=false;
     for (const auto& particle : *(genParticles.product()))
     {
        if(particle.pt()<15) continue;
        double dR=deltaR((*iTrack)->eta(),(*iTrack)->phi(),particle.eta(),particle.phi());
        if(std::abs(particle.pdgId())==13)
        {
           if(minMuDr<0||dR<minMuDr){minMuDr=dR;}
           if(dR<0.2){muClose=true;} 
        }
        else if(dR<minPartDr||minPartDr<0)
        {
           minPartDr=dR;
           minPartPdg=particle.pdgId();
        }
     }
     if(muClose){continue;}
     //if(deltaR(particleEta,particlePhi,(*iTrack)->eta(),(*iTrack)->phi())>0.2) continue;
     if(cutpro<1) cutpro=1;
     double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label,0.3,vtxHandle,*iTrack);
     double probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
     if((probeTrackIso/(*iTrack)->pt())>0.05) continue;
     if(cutpro<2) cutpro=2;
     //if(probeEcalIso>40.) continue;
     if(cutpro<3) cutpro=3;
     selectedTrack = (*iTrack); 
     TrackMatched=true;
  }
  info.minMuDr=minMuDr;
  info.particlePDG=minPartPdg;
  for(int i=0;i<cutpro;i++){allEvents.IncCutFlow(info.eventWeight);}
  if(TrackMatched)
  {
     info.probeTrackIso = myTracks.GetIsolation(iEvent,trackCollection_label,0.3,vtxHandle,selectedTrack)/selectedTrack->pt();
     info.probeReducedChi = selectedTrack->chi2()/selectedTrack->ndof();
     info.probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(selectedTrack));
     info.probeTrackEta = selectedTrack->eta();
     info.probeTrackPt = selectedTrack->pt();
     info.probeTrackPhi = selectedTrack->phi();
     info.dxy = selectedTrack->dxy();
     info.dsz = selectedTrack->dsz();
     GlobalPoint FillerPoint; 
     myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, selectedTrack,transientTrackBuilder->build(selectedTrack),  FillerPoint);
     GlobalPoint TrackGlobalPoint = GlobalPoint(GlobalPoint::Polar(selectedTrack->theta(),selectedTrack->phi(),info.bremR));
     info.minCSCImpactParameter=myCSCs.minTotalImpactParameter;
     info.cscMinDR=myCSCs.minDR;
     double staMinDr = 7.;
     double stadEoverE = 10.;
     double staE = -10;
     double closestapproach = 100.;

     if(staTracks->size()>0)
     {
       for(reco::TrackCollection::const_iterator staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack)
       {
          //reco::TransientTrack track(*staTrack, &*theMGField, theTrackingGeometry);
          reco::TransientTrack track = transientTrackBuilder->build(*staTrack); 
          double dR = deltaR(track.impactPointTSCP().momentum().eta(),track.impactPointTSCP().momentum().phi(),selectedTrack->eta(),selectedTrack->phi());
          if(dR<staMinDr)
          {
            staMinDr=dR;
            closestapproach=track.stateAtBeamLine().transverseImpactParameter().value();
          }
          if(dR<0.5&&std::sqrt(track.impactPointTSCP().momentum().mag2()>staE))
          {
            staE = std::sqrt(track.impactPointTSCP().momentum().mag2());
            stadEoverE=(std::sqrt(track.impactPointTSCP().momentum().mag2())-std::sqrt(selectedTrack->momentum().mag2()))/std::sqrt(selectedTrack->momentum().mag2());
          }
       }
     } 
     info.closestApproach=closestapproach;
     if(staMinDr<0.5) 
     {
        info.standaloneE = staE;
        info.standaloneDEoverE = stadEoverE;
     }
     info.staMinDr=staMinDr;
     if(!myHCAL.FindMuonHits(iEvent, iSetup, HBHERecHit_Label, TrackGlobalPoint, info, info.standaloneE, weight_, selectedTrack->dsz(), selectedTrack->charge(), transientTrackBuilder->build(selectedTrack))){return;}
     allEvents.IncCutFlow(info.eventWeight);
     info.nCellsFound=myHCAL.CellsFound;
     //Signal categories
     bool passSigCuts = false;
//     if(myHCAL.m_HitsOverThresh<2&&fabs(myHCAL.ieta)!=17){return;}
//     allEvents.IncCutFlow(info.eventWeight);
     //Generic HCAL hits in cone
     std::pair<double, int> coneHits = myHCAL.EnergyInCone(iEvent,iSetup,HBHERecHit_Label, TrackGlobalPoint);
     //Calo Jet Stuff
     edm::Handle<CaloJetCollection> caloJets;
     iEvent.getByToken(m_caloJet_label, caloJets);
     double nearProbeCaloJetHadE = -1.;
     double nearProbeCaloJetEmE = -1.;
     double minProbeDr=-1;
     for(auto CaloJet = caloJets->begin(); CaloJet!= caloJets->end(); ++CaloJet)
     {
        double probeDr = deltaR(selectedTrack->eta(),selectedTrack->phi(), CaloJet->eta(), CaloJet->phi());
        double hadE = CaloJet->hadEnergyInHE()+CaloJet->hadEnergyInHB();
        double ecalE = CaloJet->emEnergyInEB()+CaloJet->emEnergyInEE();
        if(probeDr<minProbeDr||minProbeDr<0){minProbeDr=probeDr;}
        if(probeDr<0.4)
        {
           
           if((hadE+ecalE)>(nearProbeCaloJetHadE+nearProbeCaloJetEmE))
           {
              nearProbeCaloJetHadE = hadE;
              nearProbeCaloJetEmE = ecalE;
           }
        }
     }
     for(int k=0;k<7;k++)
     {
        info.hitEnergies[k]=myHCAL.m_hitEnergies[k];
     }
     info.minCaloJetDr=minProbeDr;
     info.hitsOverThresh = myHCAL.m_HitsOverThresh;
     info.coneHits = myHCAL.m_coneHits;
     info.coneEnergy = myHCAL.m_coneEnergy;
     info.caloJetHcalE = nearProbeCaloJetHadE;
     info.caloJetEcalE = nearProbeCaloJetEmE;
     info.caloJetTotalE = nearProbeCaloJetHadE+nearProbeCaloJetEmE;
     info.caloConeE = coneHits.first;
     info.caloConeHits = coneHits.second;
     //PF Jets stuff
     myJets.SelectJets(iEvent,m_pfJetCollection_label);
     myJets.FindNearestJet(selectedTrack);
     info.nearestJetDr=myJets.nearestJetDr;
     info.nearestJetE=myJets.nearestJetE;
     info.nearJetMaxE=myJets.nearJetMaxE;
     //Signal cut flow
     if(nearProbeCaloJetHadE<0)
     {
        allEvents.IncCutFlow(info.eventWeight);
        if(staMinDr>1.0)
        {
           allEvents.IncCutFlow(info.eventWeight);
           if(myHCAL.m_HitsOverThresh>2)
           {
              allEvents.IncCutFlow(info.eventWeight);
              if(myHCAL.m_coneEnergy<10)
              {
                 allEvents.IncCutFlow(info.eventWeight);
                 if(myCSCs.minDR>0.05)
                 {
                    allEvents.IncCutFlow(info.eventWeight);
                 }
              }
           }
        }
     }
     bool passPreSelection=false;
     if(nearProbeCaloJetHadE<0&&myHCAL.m_HitsOverThresh>2) passPreSelection=true;
     if(myHCAL.m_coneEnergy<10&&nearProbeCaloJetHadE<0)
     {
        if(staMinDr>1.0&&myCSCs.minDR>0.05)
	{
          info.sigCat = 3;
	  passSigCuts = true;
          passTight.FillHists(info);
	}
        else if(staMinDr<1.0)
        {
           if(stadEoverE<-0.8&&myHCAL.m_coneEnergy<1.)
           {
              info.sigCat = 2;
           }
           else if(stadEoverE<-0.6&&myHCAL.m_coneEnergy<2.)
           {
              info.sigCat = 1;
	   }
           else{info.sigCat=0;}
        }
        else{info.sigCat = 0;}
     }
     else{info.sigCat=0;}
     if(!passSigCuts)
     {
        failSig.FillHists(info);
     }
     if(passPreSelection)
     {
        preSelectedEvents.FillHists(info);
     }
     allEvents.FillHists(info);
  }
}

bool PileupAnalyzer::MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, SigEventInfo info)
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
    info.genProbeMatchdE=mindE;
    info.genProbeMatchdR=mindR;
  }

  return matched;
}

int PileupAnalyzer::PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label)
{
  int cutpro = 0;
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
  {
     if(deltaR(Tag->eta(),Tag->phi(),(*iTrack)->eta(),(*iTrack)->phi())<0.2) continue;
     if(cutpro<1){cutpro=1;}
     if(Tag->charge()==(*iTrack)->charge()) continue;
     if(cutpro<2){cutpro=2;}
     int pairingResult = myTracks.PairTracks(iTrack, Tag->globalTrack(), transientTrackBuilder);
     if(pairingResult>(cutpro-2)){cutpro=pairingResult+2;}
     if(pairingResult==3)
     {
       double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label,0.3,vtxHandle,*iTrack);
       double probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
       if((probeTrackIso/(*iTrack)->pt())>0.05||probeEcalIso>10.) continue;
       if(cutpro<6){cutpro=6;}
       if(probeTrackIso<0){continue;}
       if(cutpro<7){cutpro=7;}
       selectedTrack = (*iTrack); 
       selectedMuon = Tag;
       diMuonMass = myTracks.MuonTrackMass;
       pairVtxChi = myTracks.pairvertexchi;
     } 
  }
  return cutpro;
}

bool PileupAnalyzer::MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, SigEventInfo info, double FmuE, double minCSCdr)
{
  bool matched = false;
  double mindE = 1.;
  double globalmuonE;
  double mindR = 0.1;
  if(myMuons.selectedMuons.size()==0){return matched;}
  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
  {
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
        }
        else
        {
          standaloneE = -1.;
        }
        globalmuonE = (*iMuon)->p();
     }
  }
  if(matched)
  {
    info.probeMuondE = mindE;
    info.probeMuondR = mindR;
    info.standaloneE = standaloneE;
    info.globalMuonE = globalmuonE;
  }
  else
  {
    standaloneE=0;
  }
  return matched;
}
// ------------ method called once each job just before starting event loop  ------------
void 
PileupAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PileupAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PileupAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PileupAnalyzer);
