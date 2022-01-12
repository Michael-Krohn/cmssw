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
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

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
#include "DataFormats/JetReco/interface/CaloJet.h"
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
    bool MatchTrackToMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, Muons myMuons, edm::EDGetToken m_recoMuonToken, MuPXEventInfo& info);
    bool MatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double *genInfo);
    bool LooseMatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double *genInfo);

    double GetPuWeight(edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetToken m_recoPhotonToken;
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
    bool m_signalRegion;
    double weight_;
    double nPUmean_;
    double standaloneE;
    TH1F* fPUDataHist_;
    TH1F* fPUMCHist_;
    TH1F* puWeightHist_;
    MuPXHistograms myHistogramsNoCuts;
    MuPXHistograms myHistograms;
    MuPXHistograms muProbe;
    MuPXHistograms nonMuonProbe;
    MuPXHistograms genMatched;
    MuPXHistograms nonGenMatched;
    MuPXHistograms passSig;
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
  m_recoPhotonToken (consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  reducedEndcapRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
  reducedBarrelRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
  m_pfJetCollection_label(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PFJets"))),
  m_caloJet_label(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("CaloJetSource"))),
  m_isMC (iConfig.getUntrackedParameter<bool>("isMC",false)),
  m_signalRegion(iConfig.getUntrackedParameter<bool>("signalRegion",false))
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
  TFile *MCPuFile = TFile::Open("/hdfs/cms/user/revering/dphoton/MuPlusXSkim/TTToSemiLeptonic_MCWeights.root");
  fPUMCHist_= (TH1F*)MCPuFile->Get("MCWeights/TrueNumInteractions");
  fPUDataHist_->Scale(1./fPUDataHist_->Integral());
  fPUMCHist_->Scale(1./fPUMCHist_->Integral());
  fPUDataHist_->Divide(fPUMCHist_);
  puWeightHist_=(TH1F*)fPUDataHist_->Clone();

  myHistogramsNoCuts.book(fs->mkdir("allEventsNoCuts"));
  myHistograms.book(fs->mkdir("allEvents"));
  muProbe.book(fs->mkdir("muProbe"));
  nonMuonProbe.book(fs->mkdir("nonMuonProbe"));
  genMatched.book(fs->mkdir("genMatched"));
  nonGenMatched.book(fs->mkdir("nonGenMatched"));
  passSig.book(fs->mkdir("passSig"));
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
  CSC myCSCs;

  myCSCs.minTotalImpactParameter = 1000;
  myCSCs.minDR = 1000;
  myCSCs.minTotalImpactParameter_Muon = 1000;
  myCSCs.minDR_Muon = 1000;

  info.cutProgress=0;
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
  std::vector<double> HCALisolationInfo;

  GlobalPoint MatchedGlobalPoint, VertexPosition;

  myHistogramsNoCuts.FillHists(info);

  for(std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin(); iMuon != myMuons.selectedMuons.end(); ++iMuon)
     {
        double hightrackpt=0;
        int nPairs=0;
        //if((*iMuon)->pt()<highpt) continue;
        for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
        {
	   std::cout << "Looping over tracks" << std::endl;
           if((*iMuon)->charge()==(*iTrack)->charge()) continue;
           if((*iMuon)->pt()<highpt){continue;}
           if(fabs((*iTrack)->dxy())>0.05) continue;
           if((*iMuon)->isGlobalMuon())
	   std::cout << "Pass basic selection" << std::endl;
           {
              if(myTracks.PairTracks(iTrack, (*iMuon)->globalTrack(), transientTrackBuilder))
              {
	        std::cout << "Paired track and muon" << std::endl;
                double probeTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,0.3,vtxHandle,*iTrack);
                double probeEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
                HCALisolationInfo=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*iTrack));
                if(m_signalRegion)
                {
//                   if((probeTrackIso/(*iTrack)->pt())>0.2||probeEcalIso>30) continue;
                   if((probeTrackIso/(*iTrack)->pt())>0.05) continue;
                   if(probeEcalIso>10) continue;
                   if(probeTrackIso<0) continue;
                   if(HCALisolationInfo[0] < 30) continue;
	           std::cout << "Passes HardBrem selection" << std::endl;
                }
                if(!m_signalRegion)
                {
                   if((probeTrackIso/(*iTrack)->pt())<0.2&&probeEcalIso<30) continue;
                }
                Paired=true;
		std::cout << "Paired: " << Paired << std::endl;
                if(nPairs==0){pairedTracks.clear();}
                nPairs++;
                pairedTracks.push_back((*iTrack));
		std::cout << "(*iTrack)->pt(): " << (*iTrack)->pt() << std::endl;
                if((*iTrack)->pt()>hightrackpt)
                {
		   std::cout << "myTracks.MuonTrackMass: " << myTracks.MuonTrackMass << std::endl;
                   highpt=(*iMuon)->pt();
                   pairVtxChi=myTracks.pairvertexchi;
                   muTrackMass=myTracks.MuonTrackMass;
                   selectedMuon = (*iMuon);
                   selectedTrack = (*iTrack); 
                   //hightrackpt = (*iTrack)->pt();
                   myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, myTracks.one_momentum, myTracks.tracksToVertex, myTracks.fittedVertex.position());
                }
              }
           } 
        }
        if(nPairs>0){nPairedTracks=nPairs;} 
     }

  //Check for other tracks that would have passed the global selection.
  std::cout << "2nd check, Paired: " << Paired << std::endl;
  if(Paired==true)
  {
     for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedTracks.begin(); iTrack != myTracks.selectedTracks.end(); ++iTrack)
     {  
        if(deltaR((*iTrack)->eta(),(*iTrack)->phi(),selectedTrack->eta(),selectedTrack->phi())<0.1) continue;

        if(selectedMuon->charge()==(*iTrack)->charge()) continue;
        if(myTracks.PairTracks(iTrack, selectedMuon->globalTrack(), transientTrackBuilder))
        {
	   std::cout << "Setting Paired to false" << std::endl;
           Paired=false;   
        }
     }
  }
  std::cout << "3rd check, Paired: " << Paired << std::endl;
  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
//  if(myMuons.highmuonpt<26.){return;}
  if(myMuons.highmuonpt>26.){info.cutProgress++;}
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
     info.probeTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,0.3,vtxHandle,selectedTrack);
     info.probeEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*selectedTrack));
     std::cout << "check #2" << std::endl;
     HCALisolationInfo=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*selectedTrack));
     std::cout << "HCALisolationInfo[1]: " << HCALisolationInfo[1] << std::endl;
     std::cout << "check #3" << std::endl;
     info.probeHcalIso=HCALisolationInfo[0];
     info.nHitsInCone=HCALisolationInfo[1];
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
     info.tagTrackIso=myTracks.GetIsolation(iEvent,trackCollection_label,0.3,vtxHandle,selectedMuon->innerTrack());
     info.tagVtx=myTracks.tagVtx;
     info.nVtx=myTracks.NVertices;
     info.probeVtx=myTracks.probeVtx;
     info.tagEcalIso=myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build((*selectedMuon).globalTrack()));
     std::vector<double> HCALisolationInfo2;
     HCALisolationInfo2=myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build((*selectedMuon).globalTrack()));
     info.tagHcalIso=HCALisolationInfo2[0];

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
      
    }
    info.nearestJetE=nearestJetE;
    info.nearestJetDr=nearestJetDr;
    info.nearJetMaxE=nearJetMaxE;
    info.nJets=njets;
    //Tower Isolation
    edm::Handle<CaloJetCollection> caloJets;
    iEvent.getByToken(m_caloJet_label, caloJets);
    double nearMuCaloJetE = 0;
    double nearProbeCaloJetHadE = -1;
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
       //if(info.probeCaloJetE>80){info.paired=false;}
    }
    if(selectedTrack->eta()>-2.35&&nearProbeCaloJetHadE<0)
    {
       info.sigCat=1;
    }
    else{info.sigCat=0;}
    //Nearby Photons
    edm::Handle<reco::PhotonCollection> phoHandle;
    iEvent.getByToken(m_recoPhotonToken, phoHandle);
    const reco::PhotonCollection& photon = *(phoHandle.product());
    double maxPhoEt=0;
    math::XYZVector phoMomentum;
    if(!photon.empty())
    {
       for(reco::PhotonCollection::const_iterator itr = photon.begin(); itr != photon.end(); ++itr)
       {
          if(deltaR(itr->eta(),itr->phi(),selectedTrack->eta(),selectedTrack->phi())<0.2)
          {
             if(itr->et()>maxPhoEt)
             {
                maxPhoEt=itr->et();
                phoMomentum = itr->momentum();
             }   
          }
       }
    }
    info.adjMuonTrackMass = info.muonTrackMass;
    if(maxPhoEt>0)
    {
       math::XYZVector probeMomentum = selectedTrack->momentum();
       math::XYZVector tagMomentum = selectedMuon->momentum();
       probeMomentum = probeMomentum+phoMomentum;
       double total_energy = sqrt(probeMomentum.mag2() + 0.106*0.106) + sqrt(tagMomentum.mag2() + 0.106*0.106);
       double total_px = tagMomentum.x() + probeMomentum.x();
       double total_py = tagMomentum.y() + probeMomentum.y();
       double total_pz = tagMomentum.z() + probeMomentum.z();
       info.adjMuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
    }
  }
  if(Paired)
  {
    bool matched=MatchTrackToMuon(iEvent, selectedTrack, myMuons,m_recoMuonToken, info);
    if(m_isMC)
    {
      double genInfo[4];
      bool genmatched = MatchTrackToGenMuon(iEvent, selectedTrack, selectedMuon, m_genParticleToken, genInfo);
      info.farMuE = genInfo[0];
      info.farMuDr = genInfo[1];
      info.farMuEta = genInfo[2];
      info.genTagMuDr = genInfo[3];
      if(genmatched){genMatched.FillHists(info);}
      else if(info.adjMuonTrackMass>80&&info.adjMuonTrackMass<100)
      {
        double genInfoTwo[8];
        bool loosegenmatched = LooseMatchTrackToGenMuon(iEvent,selectedTrack,selectedMuon,m_genParticleToken,genInfoTwo);
        info.farMuE = genInfoTwo[0];
        info.farMuDr = genInfoTwo[1];
        info.farMuEta = genInfoTwo[2];
        info.genTagMuDr = genInfoTwo[3];
        info.genPdgId = genInfoTwo[4];
        info.nLooseMuons = genInfoTwo[5];
        info.tagPdgId = genInfoTwo[6];
        info.nonMuPdgId = genInfoTwo[7];
        info.passLooseMatch=loosegenmatched;
        if(!info.passLooseMatch)
        {
           nonGenMatched.FillHists(info);   
           if(info.sigCat==1){passSig.FillHists(info);}
        }
      }
    } 
    if(matched)
    {
      info.globalMuonEplusCaloJetHadEplusEcalE = info.matchedProbeGlobalMuonE + info.probeHcalIso + info.probeEcalIso;
      muProbe.FillHists(info);
    }
    else{nonMuonProbe.FillHists(info);} 

    MatchedGlobalPoint = myCSCs.TrackGlobalPoint;
    VertexPosition = myCSCs.TrackVertex;
    double randphi = MatchedGlobalPoint.phi()+2.0;
    GlobalPoint RandGlobalPoint(GlobalPoint::Polar(MatchedGlobalPoint.theta(),randphi,MatchedGlobalPoint.mag()));
    bool GoodRand=true;

    std::cout << "Going into myHCAL.HitsPlots" << std::endl;
    info.ProbeConeEnergy = myHCAL.HitsPlots(iEvent, iSetup, HBHERecHit_Label, MatchedGlobalPoint, RandGlobalPoint, GoodRand, selectedTrack->pt(), VertexPosition);
    std::cout << "selectedTrack->eta(): " << selectedTrack->eta() << std::endl;

  }
  std::cout << "Filling allEvents" << std::endl;
  myHistograms.FillHists(info);
}

bool MuPXAnalyzer::MatchTrackToMuon(const edm::Event& iEvent,const reco::Track* selectedTrack, Muons myMuons, edm::EDGetToken m_recoMuonToken, MuPXEventInfo& info)
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
        info.matchedProbeGlobalMuonE = iMuon->energy();
    }
  }
  return matched;
}

bool MuPXAnalyzer::MatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double* genInfo)
{
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(m_genParticleToken, genParticles);
   int nmuons=0;
   int ntotalmuons=0;
   double farMuE=0;
   double tagMuDr=1000;
   double farMuDr=0;
   double farMuEta, farMuProbeDr;
   for (const auto& particle : *(genParticles.product()))
   {
      if (std::abs(particle.pdgId()) != 13)
         continue;
      ntotalmuons++;
      if (particle.mother()->pdgId() != 23)
         continue;
      nmuons++;
      double dR = deltaR(particle.eta(), particle.phi(), selectedTrack->eta(), selectedTrack->phi());
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

bool MuPXAnalyzer::LooseMatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selectedTrack, const reco::Muon* selectedMuon, edm::EDGetToken m_genParticleToken, double* genInfo)
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
   double farMuEta, farMuProbeDr;
   genInfo[7]=0;
   for (const auto& particle : *(genParticles.product()))
   {
      if (std::abs(particle.pdgId()) != 13)
      {
         double dR = deltaR(particle.eta(), particle.phi(), selectedTrack->eta(), selectedTrack->phi());
         if(dR<0.1){genInfo[7]=std::abs(particle.pdgId());}
         continue;
      }
   //   if (fabs(particle.mother()->pdgId()) != 13)
  //      continue;
      if (fabs(particle.mother()->pdgId())==23)
         nmuons++;
      ntotalmuons++;
      double dR = deltaR(particle.eta(), particle.phi(), selectedTrack->eta(), selectedTrack->phi());
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
   if(nPUmean_>70){lNPVW=float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(70)));}
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
