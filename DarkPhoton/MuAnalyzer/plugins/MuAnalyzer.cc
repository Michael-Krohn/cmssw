// -*- C++ -*-
//
// Package:    DarkPhoton/MuAnalyzer
// Class:      MuAnalyzer
//
/**\class MuAnalyzer MuAnalyzer_AOD.cc DarkPhoton/MuAnalyzer/plugins/MuAnalyzer_AOD.cc

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
#include <TH2.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

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

//Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Jets
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
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

class MuAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MuAnalyzer(const edm::ParameterSet&);
  ~MuAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, Histograms);
  int PairTagAndProbe(const edm::Event& iEvent,
                      const edm::EventSetup& iSetup,
                      Tracks myTracks,
                      const reco::Muon* Tag,
                      ECAL myECAL,
                      edm::ESHandle<TransientTrackBuilder> transientTrackBuilder,
                      edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                      edm::Handle<reco::VertexCollection> vtxHandle,
                      edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label,
                      edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label);
  double GetPuWeight(edm::Handle<std::vector<PileupSummaryInfo>> hPileupInfoProduct);

  edm::EDGetToken m_recoMuonToken;
  edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
  edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
  edm::EDGetToken m_genParticleToken;
  edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label;
  //    edm::EDGetTokenT<CSCRecHit2DCollection> CSCRecHitToken;
  edm::EDGetToken m_trigResultsToken;
  edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
  edm::EDGetToken m_PUInfoToken;
  std::vector<std::string> m_muonPathsToPass;
  const reco::Track* selectedTrack;
  const reco::Muon* selectedMuon;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label;
  edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label;
  edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label;
  edm::EDGetToken m_caloJet_label;
  edm::EDGetToken m_simTracksToken;

  bool m_isMC;
  bool m_runRandomTrackEfficiency;
  double weight_;
  double nPUmean_;
  double muIdWeight_;
  double muIsoWeight_;
  double muTrigWeight_;
  double muRocWeight_;
  GlobalPoint vertexPosition;
  double diMuonMass;
  bool paired_;
  TH1F* fPUDataHist_;
  TH1F* fPUMCHist_;
  TH1F* puWeightHist_;
  edm::Handle<reco::CaloJetCollection> caloJets;
  Histograms myHistograms;
  Histograms passNeighborCheck;
  Histograms minTwoHits;
  Histograms missMultiple;
  Histograms noCrossDepths;
  Histograms crossDepths;
  Histograms onEdge;
  Histograms notOnEdge;
  Histograms notOnEdgeMissingHits;
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
MuAnalyzer::MuAnalyzer(const edm::ParameterSet& iConfig)
    : m_recoMuonToken(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("recoMuons"))),
      trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      primaryVertices_Label(
          consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      CSCSegment_Label(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
      //  CSCRecHitToken(consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("CSCRecHitTag"))),
      HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>>(
          iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
      reducedEndcapRecHitCollection_Label(
          consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
      reducedBarrelRecHitCollection_Label(
          consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
      m_caloJet_label(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("CaloJetSource"))),
      m_isMC(iConfig.getUntrackedParameter<bool>("isMC", false)),
      m_runRandomTrackEfficiency(iConfig.getUntrackedParameter<bool>("runRandomTrackEfficiency", false)) {
  //now do what ever initialization is needed
  if (m_isMC) {
    m_genParticleToken = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
    m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    m_PUInfoToken = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("edmPileupInfo"));
  } else {
    m_trigResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigResults"));
    m_muonPathsToPass = iConfig.getParameter<std::vector<std::string>>("muonPathsToPass");
  }
  if (m_isMC) {
    TFile* puFile = TFile::Open(
        "/local/cms/user/revering/datafiles/PileupHists/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root");
    fPUDataHist_ = (TH1F*)puFile->Get("pileup");
    fPUDataHist_->SetDirectory(0);
    puFile->Close();
    TFile* MCPuFile = TFile::Open("/local/cms/user/revering/datafiles/PileupHists/DYJets_Private_MCWeights.root");
    fPUMCHist_ = (TH1F*)MCPuFile->Get("MCWeights/TrueNumInteractions");
    fPUDataHist_->Scale(1. / fPUDataHist_->Integral());
    fPUMCHist_->Scale(1. / fPUMCHist_->Integral());
    fPUDataHist_->Divide(fPUMCHist_);
    puWeightHist_ = (TH1F*)fPUDataHist_->Clone();
  }

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  myHistograms.book(fs->mkdir("allEvents"));
  passNeighborCheck.book(fs->mkdir("passNeighborCheck"));
  minTwoHits.book(fs->mkdir("minTwoHits"));
  missMultiple.book(fs->mkdir("missMultiple"));
  noCrossDepths.book(fs->mkdir("noCrossDepths"));
  crossDepths.book(fs->mkdir("crossDepths"));
  onEdge.book(fs->mkdir("onEdge"));
  notOnEdge.book(fs->mkdir("notOnEdge"));
  notOnEdgeMissingHits.book(fs->mkdir("notOnEdgeMissingHits"));
}

MuAnalyzer::~MuAnalyzer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void MuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  Muons myMuons;
  EventInfo info;
  Tracks myTracks;
  CSC myCSCs;
  HCAL myHCAL;
  ECAL myECAL;
  paired_ = false;
  myHistograms.m_eventCount->Fill(0.5);
  int cutProgress = 0;
  myHistograms.cutProgress = 0;
  cutProgress++;
  if (!m_isMC) {
    if (!info.passTriggers(iEvent, m_trigResultsToken, m_muonPathsToPass))
      return;
    weight_ = 1;
  } else {
    edm::Handle<GenEventInfoProduct> eventInfo;
    iEvent.getByToken(m_genInfoToken, eventInfo);
    weight_ = eventInfo->weight();
    edm::Handle<std::vector<PileupSummaryInfo>> hPileupInfoProduct;
    iEvent.getByToken(m_PUInfoToken, hPileupInfoProduct);
    assert(hPileupInfoProduct.isValid());
    double PUWeight = GetPuWeight(hPileupInfoProduct);
    weight_ = weight_ * PUWeight;
  }
  info.eventWeight = weight_;
  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);

  cutProgress++;
  int muProgress = myMuons.SelectMuons(iEvent, m_recoMuonToken);
  myHistograms.IncCutFlow(weight_, muProgress);
  if (muProgress < 4) {
    return;
  }
  int maxTrackPairProgress = -1;
  iEvent.getByToken(m_caloJet_label, caloJets);
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  for (std::vector<const reco::Muon*>::const_iterator muon = myMuons.selectedMuons.begin();
       muon != myMuons.selectedMuons.end();
       muon++) {
    if (paired_) {
      continue;
    }
    if (fabs((*muon)->eta()) > 2.4 || fabs((*muon)->eta()) < 1.45)
      continue;
    myTracks.selectedTracks.clear();
    myTracks.selectedEndcapTracks.clear();
    myTracks.trackSelected = false;
    int trackProgress = myTracks.SelectTracks(iEvent, trackCollection_label, vtxHandle, *muon, transientTrackBuilder);
    if (myTracks.trackSelected) {
      trackProgress = trackProgress + PairTagAndProbe(iEvent,
                                                      iSetup,
                                                      myTracks,
                                                      *muon,
                                                      myECAL,
                                                      transientTrackBuilder,
                                                      trackCollection_label,
                                                      vtxHandle,
                                                      reducedEndcapRecHitCollection_Label,
                                                      reducedBarrelRecHitCollection_Label);
    }
    if (trackProgress > maxTrackPairProgress) {
      maxTrackPairProgress = trackProgress;
      if (m_isMC) {
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByToken(m_genParticleToken, genParticles);
        muRocWeight_ = myMuons.MuonRoccorWeight((*muon), genParticles);
        muIdWeight_ = myMuons.MuonTightIDWeight((*muon)->pt(), (*muon)->eta())[0];
        muIsoWeight_ = myMuons.MuonTightIsoWeight((*muon)->pt(), (*muon)->eta())[0];
        muTrigWeight_ = myMuons.MuonTrigWeight((*muon)->pt(), (*muon)->eta())[0];
      }
    }
  }
  if (m_isMC) {
    weight_ = weight_ * muRocWeight_ * muIdWeight_ * muIsoWeight_ * muTrigWeight_;
  }
  myHistograms.IncCutFlow(weight_, maxTrackPairProgress);
  if (!paired_) {
    return;
  }
  //Fail events with additional tracks pairing to the selected tag muon
  for (std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedTracks.begin();
       iTrack != myTracks.selectedTracks.end();
       ++iTrack) {
    if (deltaR((*iTrack)->eta(), (*iTrack)->phi(), selectedTrack->eta(), selectedTrack->phi()) < 0.3)
      continue;
    if (selectedMuon->charge() == (*iTrack)->charge())
      continue;
    if (myTracks.PairTracks(iTrack, selectedMuon->globalTrack(), transientTrackBuilder)) {
      return;
    }
  }
  if (myTracks.MatchMuonToTrack(iEvent, trackCollection_label, vtxHandle, selectedMuon, transientTrackBuilder) < 6) {
    return;
  }
  info.probeHcalIso =
      myHCAL.GetIsolation(iEvent, iSetup, HBHERecHit_Label, transientTrackBuilder->build(*selectedTrack));
  info.tagHcalIso =
      myHCAL.GetIsolation(iEvent, iSetup, HBHERecHit_Label, transientTrackBuilder->build(myTracks.MuonMatchedTrack));
  if (info.tagHcalIso < 0 || info.probeHcalIso < 0) {
    return;
  }
  if (info.tagHcalIso < info.probeHcalIso) {
    selectedTrack = myTracks.MuonMatchedTrack;
    double temp = info.tagHcalIso;
    info.tagHcalIso = info.probeHcalIso;
    info.probeHcalIso = temp;
  }
  if (info.probeHcalIso > 20)
    return;
  myCSCs.minTotalImpactParameter = 1000;
  myCSCs.minDR = 1000;
  myCSCs.minTotalImpactParameter_Muon = 1000;
  myCSCs.minDR_Muon = 1000;
  info.diMuonMass = diMuonMass;

  myCSCs.ExtrapolateTrackToCSC(
      iEvent, iSetup, CSCSegment_Label, selectedTrack, transientTrackBuilder->build(selectedTrack), vertexPosition);
  if (myCSCs.minTotalImpactParameter > 10) {
    return;
  }
  info.minDR = myCSCs.minDR;
  for (int i = 0; i < 4; i++) {
    info.minDrByDepth[i] = myCSCs.minDrByDepth[i];
  }
  info.ProbeEta = selectedTrack->eta();
  info.ProbePhi = selectedTrack->phi();
  cutProgress++;
  myHistograms.PlotTrackDisappearance(myCSCs.TrackP,
                                      myCSCs.TrackEta,
                                      myCSCs.TrackPhi,
                                      myCSCs.minDR,
                                      myCSCs.minTotalImpactParameter,
                                      myCSCs.TrackP_dR,
                                      myCSCs.TrackEta_dR,
                                      myCSCs.TrackPhi_dR);
  //if(myCSCs.minDR > 0.1){return;}

  myHCAL.CheckHCAL(iEvent, iSetup, HBHERecHit_Label);
  GlobalPoint trackGlobalPoint = GlobalPoint(GlobalPoint::Polar(selectedTrack->theta(), selectedTrack->phi(), 525));
  std::pair<double, int> heConeEnergy = myHCAL.EnergyInCone(iEvent, iSetup, HBHERecHit_Label, trackGlobalPoint);
  /*  if(heConeEnergy.first>30)
  {
     myHistograms.m_failHEEnergy->Fill(heConeEnergy.first);
     return;
  }*/
  info.heConeEnergy = heConeEnergy.first;
  myHistograms.m_histogram_TrackerTack_P->Fill(myCSCs.TrackP);

  //  double MatchedMinDr;
  //  MatchedMinDr = 10;
  GlobalPoint MatchedGlobalPoint;
  bool TrackHCAL = true;
  if (TrackHCAL) {
    //     MatchedMinDr = myCSCs.minDR;
    MatchedGlobalPoint = myCSCs.TrackGlobalPoint;
  } else {
    //     MatchedMinDr = myCSCs.minDR_Muon;
    MatchedGlobalPoint = myCSCs.MuonGlobalPoint;
  }

  //if(MatchedMinDr < 0.1)
  if (true) {
    cutProgress++;
    double randphi = MatchedGlobalPoint.phi() + 2.0;
    if (randphi > ROOT::Math::Pi()) {
      randphi -= 2 * ROOT::Math::Pi();
    }
    if (MatchedGlobalPoint.eta() < 0 && randphi < -0.9 && randphi > -1.6) {
      randphi = randphi - 4.0 + 2 * ROOT::Math::Pi();
    }
    GlobalPoint RandGlobalPoint(GlobalPoint::Polar(MatchedGlobalPoint.theta(), randphi, MatchedGlobalPoint.mag()));
    bool GoodRand = true;
    double minDR_MuonHCAL = myHCAL.MuonMindR(iEvent, iSetup, HBHERecHit_Label, MatchedGlobalPoint);
    bool foundHits = myHCAL.HitsPlots(iEvent,
                                      iSetup,
                                      HBHERecHit_Label,
                                      MatchedGlobalPoint,
                                      RandGlobalPoint,
                                      GoodRand,
                                      myHistograms,
                                      transientTrackBuilder->build(selectedTrack));
    if (!foundHits) {
      return;
    }
    info.probeiEta = myHCAL.TrackiEta;
    info.cellsFound = myHCAL.CellsFound;
    info.etaEdgeMin = myHCAL.etaEdgeMin;
    info.phiEdgeMin = myHCAL.phiEdgeMin;
    for (int depth = 0; depth < 7; depth++) {
      info.foundDepths[depth]=myHCAL.foundDepths[depth];
      if(!info.foundDepths[depth]){continue;}
      info.expectedHits++;
      info.hitEnergies[depth] = myHCAL.m_hitEnergies[depth];
      info.neighborHitEnergies[depth] = myHCAL.m_neighborHitEnergies[depth];
      info.neighborHitEnergiesPhi[depth] = myHCAL.m_neighborHitEnergiesPhi[depth];
      info.cellCenterDEta[depth] = myHCAL.cellDetaByDepth[depth];
      info.cellCenterDPhi[depth] = myHCAL.cellDphiByDepth[depth];
      info.cellEtaEdgeDistance[depth] = myHCAL.cellEtaEdgeDistance[depth];
      info.cellPhiEdgeDistance[depth] = myHCAL.cellPhiEdgeDistance[depth];
    }
    info.centerInCone = myHCAL.centerInCone;
    info.centerOutOfCone = myHCAL.centerOutOfCone;
    info.minDR_muonHCAL = minDR_MuonHCAL;
    info.coneHits = myHCAL.m_coneHits;
    info.crossDepths = myHCAL.m_crossDepths;
    info.coneEnergy = myHCAL.m_coneEnergy;
    info.hitsOverThresh = myHCAL.m_HitsOverThresh;
    info.probePt = selectedTrack->pt();
    info.probeCharge = selectedTrack->charge();
    myHistograms.PlotCSCHits(iEvent, iSetup, CSCSegment_Label);
    myHistograms.FillHists(info);
    int nFoundHits=0;
    for(int depth=0;depth<7;depth++)
    {
       if(!info.foundDepths[depth]){continue;}
       if(info.hitEnergies[depth]>0.2)
       {
          nFoundHits++;
       }
    }
    if (info.hitsOverThresh > 1) {minTwoHits.FillHists(info);}
    if(info.expectedHits-nFoundHits>2){missMultiple.FillHists(info);}
    if(!info.crossDepths){noCrossDepths.FillHists(info);}
    else{crossDepths.FillHists(info);}
    if(myHCAL.onEdge){onEdge.FillHists(info);}
    else
    {
       notOnEdge.FillHists(info);
       if(info.expectedHits-nFoundHits>2){notOnEdgeMissingHits.FillHists(info);}
    }
    bool passAdjCheck = true;
    for (int i = 0; i < 7; i++) {
      if (info.hitEnergies[i] < 0.1) {
        if (info.neighborHitEnergies[i] > 0.1) {
          passAdjCheck = false;
        }
        if (info.neighborHitEnergiesPhi[i] > 0.1) {
          passAdjCheck = false;
        }
      }
    }
    if (passAdjCheck) {
      passNeighborCheck.FillHists(info);
    }
    //myHistograms.PlotHCALHits(iEvent,iSetup,HBHERecHit_Label);
  }
}

int MuAnalyzer::PairTagAndProbe(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup,
                                Tracks myTracks,
                                const reco::Muon* Tag,
                                ECAL myECAL,
                                edm::ESHandle<TransientTrackBuilder> transientTrackBuilder,
                                edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                                edm::Handle<reco::VertexCollection> vtxHandle,
                                edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label,
                                edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label) {
  int cutpro = 0;
  for (std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin();
       iTrack != myTracks.selectedEndcapTracks.end();
       ++iTrack) {
    if (Tag->charge() == (*iTrack)->charge())
      continue;
    if (cutpro < 1) {
      cutpro = 1;
    }
    int pairingResult = myTracks.PairTracks(iTrack, Tag->globalTrack(), transientTrackBuilder);
    if (pairingResult > (cutpro - 1)) {
      cutpro = pairingResult + 1;
    }
    if (pairingResult == 3) {
      double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label, 0.3, vtxHandle, *iTrack);
      double probeEcalIso = myECAL.GetIsolation(iEvent,
                                                iSetup,
                                                reducedEndcapRecHitCollection_Label,
                                                reducedBarrelRecHitCollection_Label,
                                                transientTrackBuilder->build(*iTrack));
      double minProbeDr = 5;
      for (auto CaloJet = caloJets->begin(); CaloJet != caloJets->end(); ++CaloJet) {
        double probeDr = deltaR((*iTrack)->eta(), (*iTrack)->phi(), CaloJet->eta(), CaloJet->phi());
        double hadE = CaloJet->hadEnergyInHE() + CaloJet->hadEnergyInHB();
        double ecalE = CaloJet->emEnergyInEB() + CaloJet->emEnergyInEE();
        if (probeDr < minProbeDr) {
          if ((hadE + ecalE) > 15) {
            minProbeDr = probeDr;
          }
        }
      }
      if ((probeTrackIso / (*iTrack)->pt()) > 0.05)
        continue;
      if (cutpro < 5) {
        cutpro = 5;
      }
      if (probeEcalIso > 10.)
        continue;
      if (cutpro < 6) {
        cutpro = 6;
      }
      //if(minProbeDr<0.4) continue;
      if (cutpro < 7) {
        cutpro = 7;
      }
      selectedTrack = (*iTrack);
      selectedMuon = Tag;
      paired_ = true;
      diMuonMass = myTracks.MuonTrackMass;
      vertexPosition = myTracks.fittedVertex.position();
    }
  }
  return cutpro;
}

bool MuAnalyzer::isTrackMatchedToMuon(const edm::Event& iEvent,
                                      std::vector<const reco::Track*>::const_iterator& Track,
                                      Histograms myHistograms) {
  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);

  bool matched = false;

  for (std::vector<reco::GenParticle>::const_iterator iGENparticle = genParticles->begin();
       iGENparticle != genParticles->end();
       iGENparticle++) {
    if (std::abs(iGENparticle->pdgId()) != 13)
      continue;

    if (deltaR((*Track)->eta(), (*Track)->phi(), iGENparticle->eta(), iGENparticle->phi()) > 0.1)
      continue;
    matched = true;
  }
  return matched;
}
// Get pileup weights
double MuAnalyzer::GetPuWeight(edm::Handle<std::vector<PileupSummaryInfo>> hPileupInfoProduct) {
  const std::vector<PileupSummaryInfo>* inPUInfos = hPileupInfoProduct.product();
  for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = inPUInfos->begin(); itPUInfo != inPUInfos->end();
       ++itPUInfo) {
    if (itPUInfo->getBunchCrossing() == 0) {
      nPUmean_ = itPUInfo->getTrueNumInteractions();
    }
  }
  float lNPVW = float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(nPUmean_)));
  if (nPUmean_ > 70) {
    lNPVW = float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(70)));
  }
  if (nPUmean_ < 1) {
    lNPVW = float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(0)));
  }
  return lNPVW;
}
// ------------ method called once each job just before starting event loop  ------------
void MuAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void MuAnalyzer::endJob() { myHistograms.Normalize(); }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuAnalyzer);
