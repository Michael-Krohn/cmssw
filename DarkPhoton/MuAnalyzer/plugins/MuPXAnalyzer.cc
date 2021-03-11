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
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1Trigger/interface/CaloSpare.h"
#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
//Triggers
#include "FWCore/Common/interface/TriggerNames.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

//Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

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
    bool MatchTrackToMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, Muons myMuons, edm::EDGetToken m_recoMuonToken);
    double GetPuWeight(edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct);

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
    edm::EDGetToken m_caloTower_label;
    edm::EDGetToken m_caloJet_label;
    edm::EDGetToken m_PUInfoToken;
    bool m_isMC;
    double weight_;
    double nPUmean_;
    double standaloneE;
    TH1F* fPUDataHist_;
    TH1F* fPUMCHist_;
    TH1F* puWeightHist_;
    MuPXHistograms myHistograms;
    MuPXHistograms muProbe;
    MuPXHistograms nonMuonProbe;
    MuPXHistograms badTagIso;
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
  m_caloTower_label(consumes<l1t::CaloTowerBxCollection>(iConfig.getParameter<edm::InputTag>("TowerSource"))),
  m_caloJet_label(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("CaloJetSource"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",false))
{
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
    m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    m_PUInfoToken = consumes< std::vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("edmPileupInfo"));
  }
  else
  {
     m_trigResultsToken = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("trigResults"));
     m_muonPathsToPass   = iConfig.getParameter<std::vector<std::string> >("muonPathsToPass");
  }
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  TFile *puFile = TFile::Open("/hdfs/cms/user/revering/dphoton/MuPlusXSkim/muPxDataPileup.root");
  fPUDataHist_ = (TH1F*)puFile->Get("pileup");
  fPUDataHist_->SetDirectory(0);
  puFile->Close();
  fPUDataHist_->Scale(1./fPUDataHist_->Integral());
  fPUMCHist_=(TH1F*)fPUDataHist_->Clone();
  float fPUMCValues[100] = {4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05, 3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473, 0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138, 0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411, 0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554, 0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895, 0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877, 0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612, 0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551, 0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934, 0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915, 0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932, 0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885, 0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012, 0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05, 2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06, 3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07, 5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07, 1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08, 6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08};
  for(int i=1;i<101; i++)
  {
     fPUMCHist_->SetBinContent(i,fPUMCValues[i-1]);
  }
  fPUMCHist_->Scale(1./fPUMCHist_->Integral());
  fPUDataHist_->Divide(fPUMCHist_);
  puWeightHist_=(TH1F*)fPUDataHist_->Clone();

  myHistograms.book(fs->mkdir("allEvents"));
  muProbe.book(fs->mkdir("muProbe"));
  nonMuonProbe.book(fs->mkdir("nonMuonProbe"));
  badTagIso.book(fs->mkdir("badTagIso"));
}


    double nPUmean_;
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
     edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct;
     iEvent.getByToken(m_PUInfoToken,hPileupInfoProduct);
     assert(hPileupInfoProduct.isValid());
     double PUweight = GetPuWeight(hPileupInfoProduct);
     weight_  = eventInfo->weight()*PUweight;
     info.nPUmean=nPUmean_;
     info.pileupWeight=PUweight;
  }
  else
  {
     weight_=1;
     nPUmean_=0;
  }
  //if(!founddpho){return;}
  info.cutProgress++;
 
  myTracks.SelectTracks(iEvent, trackCollection_label);
  myMuons.SelectMuons(iEvent, m_recoMuonToken);
  //myJets.SelectJets(iEvent, m_pfJetCollection_label);
//  myHistograms.m_NPassingTag->Fill(myMuons.selectedMuons.size(),weight_);
 
  edm::Handle<std::vector<reco::Track> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_label,thePATTrackHandle);  

  bool Paired = false;
  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  edm::Handle <reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);

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
                double probeTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,(*iTrack)->momentum().eta(),(*iTrack)->momentum().phi(),0.3,vtxHandle,*iTrack);
                double probeEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
                if(probeTrackIso/(*iTrack)->pt()<0.2&&probeEcalIso<30) continue;
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
          /* else if((*iMuon)->isTrackerMuon())
           {
              if(myTracks.PairTracks(iTrack, (*iMuon)->innerTrack(), transientTrackBuilder))
              {
                double probeTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,(*iTrack)->momentum().eta(),(*iTrack)->momentum().phi(),0.3,vtxHandle,*iTrack);
                double probeEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
                if(probeTrackIso/(*iTrack)->pt()<0.2&&probeEcalIso<30) continue;
 
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
           }*/
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
 
  if(Paired)
  {
     info.cutProgress++;
     info.muonTrackMass=muTrackMass;
     info.probeTrack=selectedTrack;
     info.tagMuon=selectedMuon;
     info.probeTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,selectedTrack->momentum().eta(),selectedTrack->momentum().phi(),0.3,vtxHandle,selectedTrack);
     info.probeEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*selectedTrack));
     info.probeHcalIso=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*selectedTrack));
     info.nPUmean=nPUmean_;
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
     info.tagVtx=myTracks.tagVtx;
     info.nVtx=myTracks.NVertices;
     info.probeVtx=myTracks.probeVtx;
     info.tagEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build((*selectedMuon).globalTrack()));
     info.tagHcalIso=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build((*selectedMuon).globalTrack()));

    //Jet analysis
    edm::Handle<reco::PFJetCollection> pfjetHandle;
    iEvent.getByToken(m_pfJetCollection_label, pfjetHandle);
    if(!pfjetHandle.isValid()){return;}
    const reco::PFJetCollection pfjets = *(pfjetHandle.product());

    int njets=0;
    double nearestJetDr=-1;
    double nearestJetE;
    double nearJetMaxE;
    //for(std::vector<const reco::PFJet*>::const_iterator it = myJets.selectedJets.begin(); it!=myJets.selectedJets.end();++it)
    for(reco::PFJetCollection::const_iterator it = pfjets.begin(); it!=pfjets.end();++it)
    {
      double NHF, NEMF, CHF, MUF, CEMF;
      int CHM,NumNeutralParticles,NumConst;
      NHF = it->neutralHadronEnergyFraction();
      NEMF = it->neutralEmEnergyFraction();
      CHF = it->chargedHadronEnergyFraction();
      MUF = it->muonEnergyFraction();
      CEMF = it->chargedEmEnergyFraction();
      NumConst = it->chargedMultiplicity()+it->neutralMultiplicity();
      NumNeutralParticles =it->neutralMultiplicity();
      CHM = it->chargedMultiplicity();
      bool JetId;
      if(fabs(it->eta())<2.6){JetId=(CEMF<0.8&&CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9);}
      else if(fabs(it->eta())<2.7){JetId= (CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 );}
      else if(fabs(it->eta())<3.0){JetId=  ( NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>1);}
      else {JetId=(NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10);}
      if(!JetId){continue;}
      njets++;
      double jetDr=deltaR(selectedTrack->eta(),selectedTrack->phi(),it->eta(),it->phi());
      if(jetDr<0.3&&it->energy()>nearJetMaxE){nearJetMaxE=it->energy();}
      if(jetDr<nearestJetDr||nearestJetDr<0)
      {
        nearestJetDr=jetDr;
        nearestJetE=it->energy();
      }
      
      //myHistograms.m_JetPt->Fill((*it)->pt(),weight_);
    }
    info.nearestJetE=nearestJetE;
    info.nearestJetDr=nearestJetDr;
    info.nearJetMaxE=nearJetMaxE;
    info.nJets=njets;
    //Tower Isolation
    edm::Handle<BXVector<l1t::CaloTower> > towerHandle;
    edm::Handle<CaloJetCollection> caloJets;
    iEvent.getByToken(m_caloJet_label, caloJets);
    double nearMuCaloJetE = 0;
    double nearProbeCaloJetHadE = 0;
    double nearProbeCaloJetEmE = 0;
    bool foundmu = false;
    bool foundprobe = false;

    for(auto CaloJet = caloJets->begin(); CaloJet!= caloJets->end(); ++CaloJet)
    {
       double muonDr = deltaR(selectedMuon->eta(), selectedMuon->phi(), CaloJet->eta(), CaloJet->phi());
       double probeDr = deltaR(selectedTrack->eta(), selectedTrack->phi(), CaloJet->eta(), CaloJet->phi());
       double hadE = CaloJet->hadEnergyInHE()+CaloJet->hadEnergyInHB();
       double ecalE = CaloJet->emEnergyInEB()+CaloJet->emEnergyInEE();
       if(muonDr<0.2)
       {
          foundmu=true;
          if((hadE+ecalE)>nearMuCaloJetE){nearMuCaloJetE = hadE+ecalE;}
       }
       if(probeDr<0.2)
       {
          foundprobe=true;
          if((hadE+ecalE)>(nearProbeCaloJetHadE+nearProbeCaloJetEmE))
          {
            nearProbeCaloJetHadE = hadE;
            nearProbeCaloJetEmE = ecalE;
          }
       }
    }
    if(foundmu){info.tagCaloJetE=nearMuCaloJetE;}
    if(foundprobe)
    {
       info.probeCaloJetE = nearProbeCaloJetHadE+nearProbeCaloJetEmE;
       info.probeCaloJetEcal = nearProbeCaloJetEmE;
       info.probeCaloJetHcal = nearProbeCaloJetHadE;
    }
/*
    edm::Handle<l1t::CaloTowerBxCollection> towers;
    iEvent.getByToken(m_caloTower_label,towers);
    double TowerE = 0;
    for(int ibx=towers->getFirstBX(); ibx<=towers->getLastBX(); ++ibx)
    {
      for(auto tow=towers->begin(ibx); tow!=towers->end(ibx); tow++)
      {
         if(tow->hwPt()<=0) {continue;}
         printf("Tower :\n  BX= %i, ipt= %i, ieta= %i, iphi= %i, etEm = %i, etHad = %i.\n", ibx, tow->hwPt(), tow->hwEta(), tow->hwPhi(), tow->hwEtEm(), tow->hwEtHad());
         double dR = deltaR(selectedTrack->eta(), selectedTrack->phi(), tow->hwEta(), tow->hwPhi());
         if(dR<0.5){TowerE+=tow->hwPt();} 
      }
    }*/
  }
  myHistograms.FillHists(info);
  bool matched = false;
  if(Paired)
  {
    matched=MatchTrackToMuon(iEvent, selectedTrack, myMuons,m_recoMuonToken);
    if(info.tagTrackIso<0||info.tagVtx!=info.probeVtx){badTagIso.FillHists(info);}
    else if(matched)
    {
      muProbe.FillHists(info);
    }
    else{nonMuonProbe.FillHists(info);} 
  }
}

bool MuPXAnalyzer::MatchTrackToMuon(const edm::Event& iEvent,const reco::Track* selectedTrack, Muons myMuons, edm::EDGetToken m_recoMuonToken)
{
  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  bool matched = false;
  double mindPt = 1.;
  if(myMuons.selectedMuons.size()==0){return matched;}
  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); ++iMuon)
  {
     if(!iMuon->isGlobalMuon()) {continue;}
     double dR = deltaR(iMuon->eta(),iMuon->phi(),selectedTrack->eta(),selectedTrack->phi());
     double dPtOverPt =  std::abs((iMuon->pt()-selectedTrack->pt())/selectedTrack->pt());
     if(dR > 0.2) continue;
     if(dPtOverPt<mindPt)
     {
        mindPt=dPtOverPt;
	matched=true;
     }
  }
  return matched;
}

double MuPXAnalyzer::GetPuWeight(edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct){
   const std::vector<PileupSummaryInfo> *inPUInfos = hPileupInfoProduct.product();
   for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = inPUInfos->begin(); itPUInfo!=inPUInfos->end(); ++itPUInfo) 
   {
      if(itPUInfo->getBunchCrossing()==0) 
      {
         //nPU      = itPUInfo->getPU_NumInteractions();
         nPUmean_  = itPUInfo->getTrueNumInteractions();
      }
   }
   float lNPVW = float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(nPUmean_)));
   if(nPUmean_>74){lNPVW=float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(74)));}
   if(nPUmean_<1){lNPVW=float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(0)));}

   return lNPVW;
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
