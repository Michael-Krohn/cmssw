// -*- C++ -*-
//
// Package:    DarkPhoton/SigAnalyzer
// Class:      SigAnalyzer
//
/**\class SigAnalyzer SigAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/SigAnalyzer.cc

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
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/transform.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

//for Standalone Muon Tracking
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

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
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

//Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DarkPhoton/MuAnalyzer/interface/MCHistograms.h"
#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"
#include "DarkPhoton/MuAnalyzer/interface/ECAL.h"
#include "DarkPhoton/MuAnalyzer/interface/Jets.h"
#include "DarkPhoton/MuAnalyzer/interface/ZptCorrection.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SigAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SigAnalyzer(const edm::ParameterSet&);
  ~SigAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, MCHistograms);
  bool MatchMotherToTrack(const edm::Event& iEvent,
                          const reco::Muon* mother,
                          edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                          SigEventInfo info);
  std::vector<int> GetGenParticlesInCone(const edm::Event& iEvent,
						    const reco::Track* selTrack,
						    edm::EDGetToken m_genParticleToken);

  double MatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selTrack, edm::EDGetToken m_genParticleToken);
  bool MatchTrackToMuon(const edm::Event& iEvent,
                        math::XYZTLorentzVectorD mother,
                        Muons myMuons,
                        SigEventInfo info,
                        double FmuE,
                        double minCSCdr);
  int PairTagAndProbe(
      const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      Tracks myTracks,
      const reco::Muon* Tag,
      ECAL myECAL,
      HCAL myHCAL,
      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
      edm::ESHandle<TransientTrackBuilder> transientTrackBuilder,
      edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
      edm::Handle<reco::VertexCollection> vtxHandle,
      edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label,
      edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label);
  bool PairMuons(const reco::Muon* Tag,
			   Muons myMuons,
			   Tracks myTracks,
                           edm::ESHandle<TransientTrackBuilder> transientTrackBuilder);
  std::vector<double> GetPuWeight(edm::Handle<std::vector<PileupSummaryInfo>> hPileupInfoProduct);

  bool passesMetFilters(const edm::Event& iEvent);
  edm::EDGetToken m_recoMuonToken;
  edm::EDGetToken m_simTracksToken;
  edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
  edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
  edm::EDGetToken m_genParticleToken;
  edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label;
  //edm::EDGetTokenT<CSCRecHit2DCollection> CSCRecHit_Label;
  edm::EDGetTokenT<DTRecSegment4DCollection> DTSegment_Label;
  edm::EDGetToken m_cscSimHits;
  //    edm::EDGetTokenT<CSCRecHit2DCollection> CSCRecHitToken;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken;
  edm::EDGetToken m_dbremweighttoken;
  edm::EDGetToken m_dbremInitialEToken;
  edm::EDGetToken m_dbremFinalEToken;
  edm::EDGetToken m_dbremAngleToken;
  edm::EDGetToken m_dbremXToken;
  edm::EDGetToken m_dbremYToken;
  edm::EDGetToken m_dbremZToken;
  edm::EDGetToken m_trigResultsToken;
  edm::EDGetToken m_trigEventToken;
  edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
  std::vector<std::string> m_muonPathsToPass;
  const reco::Track* selectedTrack;
  const reco::Track* lastSelectedTrack;
  reco::Track lastSelObj;
  const reco::Track* selectedTagTrack;
  const reco::Muon* selectedMuon;
  const reco::Muon* selectedGlobalMuon;
  const reco::Muon* matchedMuon;
  double hcalIso;
  double hcalIso0p3;
  double hcalIso0p4;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label;
  edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label;
  edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label;
  edm::EDGetToken m_caloJet_label;
  std::vector<edm::EDGetTokenT<reco::JetTagCollection>> m_jetTagTokens;
  edm::EDGetTokenT<reco::GsfElectronCollection> m_electronsToken;
  edm::EDGetToken m_PFMetToken;
  edm::EDGetToken m_PUInfoToken;
  edm::EDGetToken m_theSTAMuonLabel;
  edm::EDGetToken m_displacedSTAMuonLabel;
  edm::EDGetToken m_refittedSTAMuonLabel;
  edm::EDGetToken m_pfJetCollection_label;
  bool m_isMC;
  bool m_hasDpho;
  bool m_ecalBrem;
  bool m_sameSign;
  bool m_offPeak;
  bool m_ttbarCR;
  bool m_missingCSCFilter;
  bool m_highHcalCR;
  bool m_tagAligned;
  bool m_doubleMuon;
  bool m_studyBarrel;
  double m_scaleFactor;
  std::string m_binSFsfname;
  bool m_runRandomTrackEfficiency;
  bool noBrem;
  double standaloneE;
  double weight_;
  double nPUmean_;
  double pairVtxChi;
  double diMuonMass;
  double probeTIP;
  double probeLIP;
  double doubleMuIdWeight_;
  double doubleMuIsoWeight_;
  double doubleMuTrigWeight_;
  double doubleMuRocWeight_;
  std::map<int, std::vector<int>> metFailingEvents;
  //DBrem binning scale factors
  std::map<int, TH1F*> binningSFs;
  //pileup reweighting
  TH1F* fPUDataHist_;
  TH1F* fPUDataUpHist_;
  TH1F* fPUDataDownHist_;
  TH1F* fPUMCHist_;
  TH1F* puWeightHist_;
  TH1F* puWeightUpHist_;
  TH1F* puWeightDownHist_;
  //Met filter checking

  edm::Handle<reco::CaloJetCollection> caloJets;
  MCHistograms allEvents;
  MCHistograms noBremEvents;
  MCHistograms gen;
  MCHistograms genAccepted;
  MCHistograms preSelectedEvents;
  MCHistograms lowHits;
  MCHistograms failSig;
  MCHistograms totalDisappear;
  MCHistograms failLast;
  MCHistograms totalDisappearLowEcal;
  MCHistograms totalDisappearNoHcalIso;
  MCHistograms passSigNearbyCSC;
  MCHistograms partialDisappear;
  MCHistograms partialDisappearLowEcal;
  MCHistograms partialDisappear_largeCSCDefl;
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
SigAnalyzer::SigAnalyzer(const edm::ParameterSet& iConfig)
    : m_recoMuonToken(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("recoMuons"))),
      trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      primaryVertices_Label(
          consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      CSCSegment_Label(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
      //CSCRecHit_Label(consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("CSCRecHitsLabel"))),
      DTSegment_Label(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("DTSegmentLabel"))),
      //  CSCRecHitToken(consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("CSCRecHitTag"))),
      rpcRecHitToken(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("RPCRecHitTag"))),
      HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>>(
          iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
      reducedEndcapRecHitCollection_Label(
          consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
      reducedBarrelRecHitCollection_Label(
          consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
      m_caloJet_label(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("CaloJetSource"))),
      m_electronsToken(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
      m_PFMetToken(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("met"))),
      m_theSTAMuonLabel(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("StandAloneTracks"))),
      m_displacedSTAMuonLabel(
          consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("DisplacedStandaloneTracks"))),
      m_refittedSTAMuonLabel(
          consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("RefittedStandaloneTracks"))),
      m_pfJetCollection_label(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PFJets"))),
      m_isMC(iConfig.getUntrackedParameter<bool>("isMC", true)),
      m_hasDpho(iConfig.getUntrackedParameter<bool>("hasDpho", false)),
      m_ecalBrem(iConfig.getUntrackedParameter<bool>("ecalBrem", false)),
      m_sameSign(iConfig.getUntrackedParameter<bool>("sameSign", false)),
      m_offPeak(iConfig.getUntrackedParameter<bool>("offPeakCR", false)),
      m_ttbarCR(iConfig.getUntrackedParameter<bool>("ttbarCR", false)),
      m_missingCSCFilter(iConfig.getUntrackedParameter<bool>("missingCSC", false)),
      m_highHcalCR(iConfig.getUntrackedParameter<bool>("highHcalCr", false)),
      m_tagAligned(iConfig.getUntrackedParameter<bool>("tagAlignedCr", false)),
      m_doubleMuon(iConfig.getUntrackedParameter<bool>("doubleMuonCr", false)),
      m_studyBarrel(iConfig.getUntrackedParameter<bool>("studyBarrelMuons", false)),
      m_scaleFactor(iConfig.getUntrackedParameter<double>("scaleFactor", 2000)),
      m_binSFsfname(iConfig.getUntrackedParameter<std::string>("binSFsfname", "")),
      m_runRandomTrackEfficiency(iConfig.getUntrackedParameter<bool>("runRandomTrackEfficiency", false)) {
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  //now do what ever initialization is needed
  if (m_isMC) {
    m_genParticleToken = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
    m_simTracksToken = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("g4SimHits"));
    m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    m_PUInfoToken = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("edmPileupInfo"));
    m_cscSimHits = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("cscSimHits"));
    if (m_hasDpho) {
      m_dbremweighttoken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremweighttag"));
      m_dbremInitialEToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremInitialETag"));
      m_dbremFinalEToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremFinalETag"));
      m_dbremAngleToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremAngleTag"));
      m_dbremXToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremXTag"));
      m_dbremYToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremYTag"));
      m_dbremZToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremZTag"));
      TFile* dbremBinSFs = TFile::Open(m_binSFsfname.c_str(), "READ");
      TTree* energyBins = (TTree*)dbremBinSFs->Get("energyBins");
      int energy;
      energyBins->SetBranchAddress("energies", &energy);
      for (int iEntry = 0; energyBins->LoadTree(iEntry) >= 0; ++iEntry) {
        energyBins->GetEntry(iEntry);
        std::string histName = "energySF_ebeam" + std::to_string(energy);
        TH1F* energyBinSFs = (TH1F*)dbremBinSFs->Get(histName.c_str());
        energyBinSFs->SetDirectory(0);
        binningSFs[energy] = energyBinSFs;
      }
      dbremBinSFs->Close();
    }
  }
  m_trigResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigResults"));
  m_trigEventToken = consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD"));
  m_muonPathsToPass = iConfig.getParameter<std::vector<std::string>>("muonPathsToPass");
  m_jetTagTokens =
      edm::vector_transform(iConfig.getParameter<std::vector<edm::InputTag>>("btagAlgos"),
                            [this](edm::InputTag const& tag) { return mayConsume<reco::JetTagCollection>(tag); });

  if (m_isMC) {
    TFile* puFile = TFile::Open(
        "/local/cms/user/revering/datafiles/PileupHists/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root");
    fPUDataHist_ = (TH1F*)puFile->Get("pileup");
    fPUDataHist_->SetDirectory(0);
    puFile->Close();
    TFile* puUpFile = TFile::Open(
        "/local/cms/user/revering/datafiles/PileupHists/PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root");
    fPUDataUpHist_ = (TH1F*)puUpFile->Get("pileup");
    fPUDataUpHist_->SetDirectory(0);
    puUpFile->Close();
    TFile* puDownFile = TFile::Open(
        "/local/cms/user/revering/datafiles/PileupHists/PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root");
    fPUDataDownHist_ = (TH1F*)puDownFile->Get("pileup");
    fPUDataDownHist_->SetDirectory(0);
    puDownFile->Close();

    TFile* MCPuFile = TFile::Open("/local/cms/user/revering/datafiles/PileupHists/DYJets_Private_MCWeights.root");
    fPUMCHist_ = (TH1F*)MCPuFile->Get("MCWeights/TrueNumInteractions");
    fPUDataHist_->Scale(1. / fPUDataHist_->Integral());
    fPUDataUpHist_->Scale(1. / fPUDataUpHist_->Integral());
    fPUDataDownHist_->Scale(1. / fPUDataDownHist_->Integral());
    fPUMCHist_->Scale(1. / fPUMCHist_->Integral());
    fPUDataHist_->Divide(fPUMCHist_);
    fPUDataUpHist_->Divide(fPUMCHist_);
    fPUDataDownHist_->Divide(fPUMCHist_);
    puWeightHist_ = (TH1F*)fPUDataHist_->Clone();
    puWeightUpHist_ = (TH1F*)fPUDataUpHist_->Clone();
    puWeightDownHist_ = (TH1F*)fPUDataDownHist_->Clone();
  }

  noBrem=false;
  lastSelectedTrack=NULL;

  std::string metFname = iConfig.getUntrackedParameter<std::string>("metFailName", "");
  if (metFname != "") {
    TFile* metFailingFile = TFile::Open(iConfig.getUntrackedParameter<std::string>("metFailName", "").c_str());
    TTree* metFailingTree = (TTree*)metFailingFile->Get("myMetFilter/failingEvents/failingEventInfo");
    int metFailLumi, metFailEvent;
    metFailingTree->SetBranchAddress("lumi", &metFailLumi);
    metFailingTree->SetBranchAddress("event", &metFailEvent);
    for (int eventNum = 0; eventNum < metFailingTree->GetEntries(); eventNum++) {
      metFailingTree->GetEntry(eventNum);
      if (metFailingEvents.count(metFailLumi) == 0) {
        metFailingEvents[metFailLumi];
      }
      metFailingEvents[metFailLumi].push_back(metFailEvent);
    }
    metFailingFile->Close();
  }

  allEvents.book(fs->mkdir("allEvents"), m_isMC);
  noBremEvents.book(fs->mkdir("noBremEvents"), m_isMC);
  genAccepted.book(fs->mkdir("genAccepted"), m_isMC);
  gen.book(fs->mkdir("genInfo"), m_isMC);
  preSelectedEvents.book(fs->mkdir("preSelected"), m_isMC);
  failSig.book(fs->mkdir("failSig"), m_isMC);
  lowHits.book(fs->mkdir("lowHits"), m_isMC);
  totalDisappear.book(fs->mkdir("totalDisappear"), m_isMC, 1);
  failLast.book(fs->mkdir("failLastHE"), m_isMC, 1);
  totalDisappearLowEcal.book(fs->mkdir("totalDisappearLowEcal"), m_isMC);
  totalDisappearNoHcalIso.book(fs->mkdir("totalDisappearNoHcalIso"), m_isMC);
  passSigNearbyCSC.book(fs->mkdir("totalDisappearSigNearbyCSC"), m_isMC, 1);
  partialDisappear.book(fs->mkdir("partialDisappear"), m_isMC);
  partialDisappearLowEcal.book(fs->mkdir("oartialDisappearLowEcal"), m_isMC);
  partialDisappear_largeCSCDefl.book(fs->mkdir("partialDisappear_largeCSCDefl"), m_isMC);
}

SigAnalyzer::~SigAnalyzer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void SigAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  Muons myMuons;
  Jets myJets;
  EventInfo myEventInfo;
  Tracks myTracks;
  CSC myCSCs(m_missingCSCFilter);
  HCAL myHCAL;
  ECAL myECAL;
  SigEventInfo info;
  //if(!m_isMC){
  if (!info.passTriggers(iEvent, m_trigResultsToken, m_trigEventToken, m_muonPathsToPass))
    return;
  //}
  noBrem=false;
  allEvents.ResetCutFlow();
  noBremEvents.ResetCutFlow();
  genAccepted.ResetCutFlow();
  totalDisappear.ResetCutFlow();
  passSigNearbyCSC.ResetCutFlow();
  edm::Handle<float> dbremweight, dbremInitialE, dbremFinalE, dbremAngle, dbremX, dbremY, dbremZ;
  bool genAccept = false;
  if (m_isMC) {
    if (m_hasDpho) {
      iEvent.getByToken(m_dbremweighttoken, dbremweight);
      edm::Handle<GenEventInfoProduct> eventInfo;
      iEvent.getByToken(m_genInfoToken, eventInfo);
      weight_ = eventInfo->weight();
      allEvents.IncCutFlow(weight_);
      noBremEvents.IncCutFlow(weight_);
      if (*dbremweight>0) 
      {
        iEvent.getByToken(m_dbremInitialEToken, dbremInitialE);
        iEvent.getByToken(m_dbremAngleToken, dbremAngle);
        iEvent.getByToken(m_dbremXToken, dbremX);
        iEvent.getByToken(m_dbremYToken, dbremY);
        iEvent.getByToken(m_dbremZToken, dbremZ);
        iEvent.getByToken(m_dbremFinalEToken, dbremFinalE);
        info.dbremWeight = *dbremweight * m_scaleFactor;
	if(*dbremweight>0){weight_ = *dbremweight * weight_;}
        double hyp = sqrt(*dbremX * (*dbremX) + *dbremY * (*dbremY) + *dbremZ * (*dbremZ));
        double dt = sqrt(*dbremX * (*dbremX) + *dbremY * (*dbremY));
        double sinTheta = dt / hyp;
        double et = *dbremInitialE * sinTheta;
        info.dbremInitialEt = et;
        info.dbremFinalE = *dbremFinalE;
        info.dbremInitialE = *dbremInitialE;
	TLorentzVector initialMu;
	double initialP = sqrt(info.dbremInitialE*info.dbremInitialE-0.01);
	TVector3 initialDir(*dbremX,*dbremY,*dbremZ);
	initialMu.SetPxPyPzE(0,0,initialP,info.dbremInitialE);
	//initialMu.RotateUz(initialDir.Unit());
	//info.dbremInitialPt = initialMu.Pt();
	//info.dbremInitialEta = initialMu.Eta();
	TLorentzVector finalMu;
	double randPhi = rand()*2.*3.141;
	double dbremP = 0;
	if(info.dbremFinalE>0.1){dbremP = sqrt(info.dbremFinalE*info.dbremFinalE-0.01);}

	finalMu.SetPxPyPzE(dbremP*std::sin(*dbremAngle)*std::cos(randPhi),dbremP*std::sin(*dbremAngle)*std::sin(randPhi),dbremP*std::cos(*dbremAngle),info.dbremFinalE);
	//finalMu.RotateUz(initialDir.Unit());
	//info.dbremFinalPt = finalMu.Pt();
	//info.dbremFinalEta = finalMu.Eta();
        GlobalPoint vtx(*dbremX / 10., *dbremY / 10., *dbremZ / 10.);
        info.vtxZ = fabs(vtx.z());
        info.vtxPerp = vtx.perp();
        //if(fabs(vtx.z())>388&&vtx.perp()<268&&vtx.perp()<(1.2857*fabs(vtx.z()-302.8)))
        //info.bremR = vtx.mag();
        if(m_ecalBrem){
          bool goodZone=false;
          if(info.vtxZ>220&&info.vtxZ<300&&info.vtxPerp<150&&info.vtxPerp>125){goodZone=true;}
          else if(info.vtxZ>299&&info.vtxZ<309&&info.vtxPerp>50&&info.vtxPerp<126){goodZone=true;}
          else if(info.vtxZ>315&&info.vtxZ<345&&info.vtxPerp>55&&info.vtxPerp<160){goodZone=true;}
          if(!goodZone)
          {
    	   info.vtxPerp=info.vtxPerp/10.;
               info.vtxZ=info.vtxZ/10.;
    	   failSig.FillHists(info);
    	   return;
          }
        }
        for (auto const& iter : binningSFs) {
          //Use the first bin with energy > initial energy.
          if (iter.first > *dbremInitialE) {
            info.usedBin = iter.first;
            info.energyBinSF = iter.second->GetBinContent(iter.second->FindBin(*dbremFinalE / *dbremInitialE));
            break;
          }
          gen.FillHists(info);
        }
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByToken(m_genParticleToken, genParticles);
        bool maybeTag = false;
        bool maybeProbe = false;
        for (const auto& particle : *(genParticles.product())) {
          if (std::abs(particle.pdgId()) != 13) {
            continue;
          }
          if (particle.mother()->pdgId() != 23) {
            continue;
          }
          if (fabs(particle.eta()) > 2.4) {
            continue;
          }
          if (particle.pt() < 20) {
            continue;
          }
          if (maybeProbe && particle.pt() > 26) {
            genAccept = true;
            break;
          }
          if (fabs(particle.eta() > 1.6)) {
            maybeProbe = true;
            if (maybeTag) {
              genAccept = true;
              break;
            }
          }
          if (particle.pt() > 26) {
            maybeTag = true;
          }
        }
      } 
      else{noBrem=true;}
    } 
    else 
    {
      weight_ = 1.;
      allEvents.IncCutFlow(weight_);
    }
  } else {
    weight_ = 1.;
    allEvents.IncCutFlow(weight_);
  }
  //Zpt Correction
  if(m_isMC){
     std::string ZptPath = "/local/cms/user/revering/dphoton/slc7/CMSSW_10_6_17_patch1/src/DarkPhoton/MuAnalyzer/ZptCorrection/ZptWeight_DYJets_MG.root";
     ZptCorrection myZptCorrection(ZptPath);

     edm::Handle<reco::GenParticleCollection> genParticles;
     iEvent.getByToken(m_genParticleToken, genParticles);
     bool foundZ = false;
     info.foundZ=false;
     double Zpt = 0;
     for (const auto& particle : *(genParticles.product())) 
     {
       if (std::abs(particle.pdgId()) != 23) 
       {
         continue;
       }
       foundZ=true;
       if(particle.pt()>Zpt){Zpt=particle.pt();}
     }
     if(foundZ)
     {
        info.ZptWeight=myZptCorrection.GetZptWeight(Zpt);
        info.Zpt = Zpt;
        info.eventWeight = info.eventWeight * info.ZptWeight;
	info.foundZ=true;
     }
  }

  if (!passesMetFilters(iEvent)) {
    return;
  }
  if(!noBrem){allEvents.IncCutFlow(weight_);}  //Passes met filters
  noBremEvents.IncCutFlow(weight_);

  if (m_isMC) {
    edm::Handle<std::vector<PileupSummaryInfo>> hPileupInfoProduct;
    iEvent.getByToken(m_PUInfoToken, hPileupInfoProduct);
    assert(hPileupInfoProduct.isValid());
    std::vector<double> PUWeight = GetPuWeight(hPileupInfoProduct);
    weight_ = weight_ * PUWeight[0];
    info.PuUpDownWeights.push_back(PUWeight[1]);
    info.PuUpDownWeights.push_back(PUWeight[2]);
    info.nPUmean = nPUmean_;
    info.pileupWeight = PUWeight[0];
  }
  info.eventWeight = weight_;
  info.genEventWeight = weight_;
  if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
  noBremEvents.IncCutFlow(info.eventWeight);
  if (genAccept) {
    genAccepted.IncCutFlow(info.eventWeight);
  }
  edm::Handle<reco::TrackCollection> staTracks;
  iEvent.getByToken(m_theSTAMuonLabel, staTracks);
  edm::Handle<reco::TrackCollection> displacedStaTracks;
  iEvent.getByToken(m_displacedSTAMuonLabel, displacedStaTracks);
  edm::Handle<reco::TrackCollection> refittedStaTracks;
  iEvent.getByToken(m_refittedSTAMuonLabel, refittedStaTracks);

  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);
  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  int muProgress = myMuons.SelectMuons(iEvent, m_recoMuonToken);
  for (int i = 0; i < muProgress; i++) {
    if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
    noBremEvents.IncCutFlow(info.eventWeight);

    if (genAccept) {
      genAccepted.IncCutFlow(info.eventWeight);
    }
  }
  if (muProgress < 4) {
    return;
  }
  bool SiblingMatched = false;
  int maxTrackProgress = -1;
  iEvent.getByToken(m_caloJet_label, caloJets);
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  for (std::vector<const reco::Muon*>::const_iterator muon = myMuons.selectedMuons.begin();
       muon != myMuons.selectedMuons.end();
       muon++) {
    if (maxTrackProgress == 13) {
      continue;
    }
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
                                                      myHCAL,
                                                      HBHERecHit_Label,
                                                      transientTrackBuilder,
                                                      trackCollection_label,
                                                      vtxHandle,
                                                      reducedEndcapRecHitCollection_Label,
                                                      reducedBarrelRecHitCollection_Label);
    }
    if (trackProgress > maxTrackProgress) {
      maxTrackProgress = trackProgress;
      if (m_isMC) {
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByToken(m_genParticleToken, genParticles);
        info.filledWeights = true;
        info.muRocWeight = myMuons.MuonRoccorWeight((*muon), genParticles);
        info.muIdWeights = myMuons.MuonTightIDWeight((*muon)->pt(), (*muon)->eta());
        info.muIsoWeights = myMuons.MuonTightIsoWeight((*muon)->pt(), (*muon)->eta());
        info.muTrigWeights = myMuons.MuonTrigWeight((*muon)->pt(), (*muon)->eta());
      }
    }
  }
  bool pairedMuons_ = false;
  if (m_doubleMuon){
    for (std::vector<const reco::Muon*>::const_iterator muon = myMuons.selectedMuons.begin();
         muon != myMuons.selectedMuons.end();
         muon++) {
      if (pairedMuons_) {
        continue;
      }
      pairedMuons_ = PairMuons(*muon,myMuons,myTracks,transientTrackBuilder);
      if (pairedMuons_){
        if (m_isMC) {
          edm::Handle<reco::GenParticleCollection> genParticles;
          iEvent.getByToken(m_genParticleToken, genParticles);
          doubleMuRocWeight_ = myMuons.MuonRoccorWeight((*muon), genParticles);
          doubleMuIdWeight_ = myMuons.MuonTightIDWeight((*muon)->pt(), (*muon)->eta())[0];
          doubleMuIsoWeight_ = myMuons.MuonTightIsoWeight((*muon)->pt(), (*muon)->eta())[0];
          doubleMuTrigWeight_ = myMuons.MuonTrigWeight((*muon)->pt(), (*muon)->eta())[0];
        }
      }
    }
  }
  if (m_isMC) {
    if (m_doubleMuon){
      info.eventWeight = info.eventWeight * doubleMuRocWeight_ * doubleMuIdWeight_ * doubleMuIsoWeight_ * doubleMuTrigWeight_;
    }else{
      info.eventWeight =
        info.eventWeight * info.muIdWeights[0] * info.muIsoWeights[0] * info.muTrigWeights[0] * info.muRocWeight;
    }
  }
  if (maxTrackProgress == 13) {
    SiblingMatched = true;
    if (m_isMC) {
      std::vector<double> probeIsoWeights = myMuons.MuonTightIsoWeight(selectedTrack->pt(), selectedTrack->eta());
      info.eventWeight = info.eventWeight * probeIsoWeights[0];
      info.muIsoWeights[0] = info.muIsoWeights[0] * probeIsoWeights[0];
      info.muIsoWeights[1] = info.muIsoWeights[1] * probeIsoWeights[1];
      info.muIsoWeights[2] = info.muIsoWeights[2] * probeIsoWeights[2];
    }
  }
  for (int i = 0; i < maxTrackProgress; i++) {
    if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
    noBremEvents.IncCutFlow(info.eventWeight);
    if (genAccept) {
      genAccepted.IncCutFlow(info.eventWeight);
    }
  }
  if(maxTrackProgress >=12){
    if(diMuonMass>80&&diMuonMass<100){
      allEvents.FillBaseEvents(info);}
  }
  if (maxTrackProgress < 13) {
    return;
  }
  //Make sure the tag muon passed the trigger
  if (!info.matchTagAndTrigger(selectedMuon)) {
    return;
  }
  if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
  noBremEvents.IncCutFlow(info.eventWeight);

  math::XYZTLorentzVectorD mother;
  //Find deflected muon or the initial muon if there was no deflection
  double FmuE = 0;
  if (m_isMC && m_hasDpho &&!noBrem) {
    info.finalMuE = *dbremFinalE;
    FmuE = *dbremFinalE;
  }

  info.nTagMuons = myMuons.selectedMuons.size();
  //Check if another track also passes the pairing requirements.
  if (SiblingMatched) {
    for (std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedTracks.begin();
         iTrack != myTracks.selectedTracks.end();
         ++iTrack) {
      if (deltaR((*iTrack)->eta(), (*iTrack)->phi(), selectedTrack->eta(), selectedTrack->phi()) < 0.3)
        continue;
      if (selectedMuon->charge() == (*iTrack)->charge())
        continue;
      //if(fabs((*iTrack)->dxy())>0.05) continue;
      if (myTracks.PairTracks(iTrack, selectedMuon->globalTrack(), transientTrackBuilder)) {
        double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label, 0.3, vtxHandle, *iTrack);
        double probeEcalIso = myECAL.GetIsolation(iEvent,
                                                  iSetup,
                                                  reducedEndcapRecHitCollection_Label,
                                                  reducedBarrelRecHitCollection_Label,
                                                  transientTrackBuilder->build(*iTrack),
                                                  0.4);
        if ((probeTrackIso / (*iTrack)->pt()) > 0.2 || probeEcalIso > 30.)
          continue;
        if (probeTrackIso < 0) {
          continue;
        }
        SiblingMatched = false;
      }
    }
    //Increment cut flow if no other tracks paired
    if (SiblingMatched) {
      if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
      noBremEvents.IncCutFlow(info.eventWeight);
    }
    //Check if another global muon also passes the pairing requirements
    for (std::vector<const reco::Muon*>::const_iterator muon = myMuons.globalMuons.begin();
         muon != myMuons.globalMuons.end();
         muon++) {
      if (deltaR(selectedTrack->eta(), selectedTrack->phi(), (*muon)->eta(), (*muon)->phi()) > 1.0) {
        if (myTracks.PairMuons((*muon)->muonBestTrack(), selectedMuon->muonBestTrack(), transientTrackBuilder)) {
          SiblingMatched = false;
        }
      }
    }
    mother.SetPxPyPzE(selectedTrack->momentum().x(),
                      selectedTrack->momentum().y(),
                      selectedTrack->momentum().z(),
                      selectedTrack->p());
  }
  if(!info.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
  if (SiblingMatched) {
    //Accept events that pass preselection w/o passing gen acceptance
    if (!genAccept) {
      genAccept = true;
      for (int k = 0; k < 16; k++) {
        genAccepted.IncCutFlow(info.eventWeight);
      }
    }
    //Increment the cut flow if no other global muons paired
    if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
    noBremEvents.IncCutFlow(info.eventWeight);

    if (genAccept) {
      genAccepted.IncCutFlow(info.eventWeight);
    }
    if(m_tagAligned)
    {
       if(!MatchMotherToTrack(iEvent,
                          selectedMuon,
                          trackCollection_label,
                          info)){return;}
       selectedTrack = selectedTagTrack;
    }else if(m_doubleMuon)
    {
       if(!MatchMotherToTrack(iEvent,
                          selectedGlobalMuon,
                          trackCollection_label,
                          info)){return;}
       selectedTrack = selectedTagTrack;
    }
    info.probeTrackIso =
        myTracks.GetIsolation(iEvent, trackCollection_label, 0.3, vtxHandle, selectedTrack) / selectedTrack->pt();
    info.isLoose = selectedTrack->quality(reco::Track::qualityByName("loose"));
    info.isTight = selectedTrack->quality(reco::Track::qualityByName("tight"));
    info.isHighPurity = selectedTrack->quality(reco::Track::qualityByName("highPurity"));
    info.probeReducedChi = selectedTrack->chi2() / selectedTrack->ndof();
    info.tagMuonEta = selectedMuon->eta();
    info.pairVtxChi = pairVtxChi;
    info.probeEcalIso = myECAL.GetIsolation(iEvent,
                                            iSetup,
                                            reducedEndcapRecHitCollection_Label,
                                            reducedBarrelRecHitCollection_Label,
                                            transientTrackBuilder->build(selectedTrack),
                                            0.4);
    info.probeEcalIso0p1 = myECAL.GetIsolation(iEvent,
                                               iSetup,
                                               reducedEndcapRecHitCollection_Label,
                                               reducedBarrelRecHitCollection_Label,
                                               transientTrackBuilder->build(selectedTrack),
                                               0.1);
    info.probeEcalIso0p2 = myECAL.GetIsolation(iEvent,
                                               iSetup,
                                               reducedEndcapRecHitCollection_Label,
                                               reducedBarrelRecHitCollection_Label,
                                               transientTrackBuilder->build(selectedTrack),
                                               0.2);
    info.probeEcalIso0p05 = myECAL.GetIsolation(iEvent,
                                                iSetup,
                                                reducedEndcapRecHitCollection_Label,
                                                reducedBarrelRecHitCollection_Label,
                                                transientTrackBuilder->build(selectedTrack),
                                                0.05);
    info.tagMuonPt = selectedMuon->pt();
    info.probeTrackEta = selectedTrack->eta();
    info.probeCharge = selectedTrack->charge();
    info.probeTrackP = selectedTrack->p();
    info.probeTrackPt = selectedTrack->pt();
    info.probeTrackPhi = selectedTrack->phi();
    info.dxy = selectedTrack->dxy();
    info.dca = probeTIP;
    info.dcal = probeLIP;
    info.dsz = selectedTrack->dsz();
    info.diMuonMass = diMuonMass;
    GlobalPoint FillerPoint;
    myCSCs.ExtrapolateTrackToCSC(
         iEvent, iSetup, CSCSegment_Label, selectedTrack, transientTrackBuilder->build(selectedTrack), FillerPoint);
/*    myCSCs.ExtrapolateTrackToCSC(iEvent,
                                 iSetup,
                                 CSCSegment_Label,
				 CSCRecHit_Label,
                                 rpcRecHitToken,
                                 selectedTrack,
                                 transientTrackBuilder->build(selectedTrack),
                                 FillerPoint);*/
    myCSCs.ExtrapolateTrackToDTs(iEvent,
			                 iSetup,
			                 DTSegment_Label,
			                 rpcRecHitToken,
			                 selectedTrack,
			                 transientTrackBuilder->build(selectedTrack),
			                 FillerPoint);
    info.cscRecHitMinDR = myCSCs.recMinDR;
    info.cscRpcHitMinDR = myCSCs.rpcMinDR;
    info.CSCdEta = myCSCs.CSCdEta;
    info.CSCdPhi = myCSCs.CSCdPhi;
    int nCloseHits = 0;
/*    for(auto closeDr : *myCSCs.m_closeDrs)
    {
        info.closeDrs->push_back(closeDr);
	nCloseHits++;
    }*/
    //Find number of good depths for 1 depth -> 0 depths MC correction
    int nGoodDepths = 0;
    int goodDepth = 0;
    for (int i = 0; i < 4; i++) {
      info.nCloseCSCHits[i] = myCSCs.nCloseHits[i];
      info.minDrByStation[i] = myCSCs.minDrByDepth[i];
      info.minDtDrByStation[i] = myCSCs.minDtDrByDepth[i];
      info.minDPhiByStation[i] = myCSCs.minDPhiByDepth[i];
      info.secMinDPhiByStation[i] = myCSCs.secMinDPhiByDepth[i];
      info.secMinDrByStation[i] = myCSCs.secMinDrByDepth[i];
      info.minDtDPhiByStation[i] = myCSCs.minDtDPhiByDepth[i];
      info.minDEtaByStation[i] = myCSCs.minDEtaByDepth[i];
      info.minDtDEtaByStation[i] = myCSCs.minDtDEtaByDepth[i];
      info.minDZByStation[i] = myCSCs.minDZByDepth[i];
      info.minDtDZByStation[i] = myCSCs.minDtDZByDepth[i];

      if (myCSCs.minDrByDepth[i] < 0.05 && myCSCs.minDrByDepth[i] > 0) {
        nGoodDepths++;
        goodDepth = i;
      }
    }
    double cscScaleFactor = 0;
    if (m_missingCSCFilter) {
      if (nGoodDepths == 1) {
        cscScaleFactor = myCSCs.GetMissingHitWeight(goodDepth, selectedTrack->eta(), selectedTrack->phi());
      }
    }
    GlobalPoint TrackGlobalPoint =
        GlobalPoint(GlobalPoint::Polar(selectedTrack->theta(), selectedTrack->phi(), 525));
    if (m_isMC && m_hasDpho && !noBrem) {
      info.deflectingAngle = *dbremAngle;
      info.dbremFinalE = *dbremFinalE;
    }
    info.minCSCImpactParameter = myCSCs.minTotalImpactParameter;
    info.cscMinDR = myCSCs.minDR;
    if (FmuE == 0) {
      FmuE = selectedTrack->p();
    }
    bool MuonMatched = MatchTrackToMuon(iEvent, mother, myMuons, info, FmuE, info.cscMinDR);
    double staMinDr = 7.;
    if (MuonMatched) {
      info.cscChiSq = myCSCs.CSCChiSq;
      int pairingResult = myTracks.PairTracks(selectedMuon->innerTrack(), matchedMuon->globalTrack(), transientTrackBuilder);
      if(pairingResult==3){
         info.altInvMass = myTracks.MuonTrackMass;
      }
      ROOT::Math::DisplacementVector3D unitmother = mother.Vect().unit();
      GlobalVector unitcsc = myCSCs.CSCTraj.unit();
      double dotpro = unitmother.x() * unitcsc.x() + unitmother.y() * unitcsc.y() + unitmother.z() * unitcsc.z();
      info.cscHitAngleChange = std::acos(dotpro);
    }
    if (m_isMC) {
      info.minGenMuDr = MatchTrackToGenMuon(iEvent, selectedTrack, m_genParticleToken);
      info.pdgIds = GetGenParticlesInCone(iEvent, selectedTrack, m_genParticleToken); 
      myCSCs.ExtrapolateTrackToSimCSC(iEvent, iSetup, m_cscSimHits, selectedTrack);
      info.nSimCSCHits = myCSCs.nSimCSCHits;
      info.minSimCSCDr = myCSCs.minSimCSCDr;
    }
    double stadEoverE = 10.;
    double staE = -10;
    double closestapproach = 100.;

    if (staTracks->size() > 0) {
      for (reco::TrackCollection::const_iterator staTrack = staTracks->begin(); staTrack != staTracks->end();
           ++staTrack) {
        //reco::TransientTrack track(*staTrack, &*theMGField, theTrackingGeometry);
        reco::TransientTrack track = transientTrackBuilder->build(*staTrack);
        double dR = deltaR(track.impactPointTSCP().momentum().eta(),
                           track.impactPointTSCP().momentum().phi(),
                           selectedTrack->eta(),
                           selectedTrack->phi());
        if (dR < staMinDr) {
          staMinDr = dR;
          closestapproach = track.stateAtBeamLine().transverseImpactParameter().value();
        }
        //Pick the largest energy standalone muon within dR of 1.0
        if (dR < 1.0 && std::sqrt(track.impactPointTSCP().momentum().mag2()) > staE) {
          if (track.ndof() < 1) {
            continue;
          }
          //if(track.chi2()/track.ndof()>10){continue;}
          staE = std::sqrt(track.impactPointTSCP().momentum().mag2());
          info.stadEta = track.impactPointTSCP().momentum().eta() - selectedTrack->eta();
          info.stadPhi = track.impactPointTSCP().momentum().phi() - selectedTrack->phi();
          stadEoverE = (std::sqrt(track.impactPointTSCP().momentum().mag2()) - mother.P()) / mother.P();
          info.staChi2 = track.chi2() / track.ndof();
          int nstaHits = 0;
          trackingRecHit_iterator rhbegin = staTrack->recHitsBegin();
          trackingRecHit_iterator rhend = staTrack->recHitsEnd();
          for (trackingRecHit_iterator recHit = rhbegin; recHit != rhend; ++recHit) {
            nstaHits++;
          }
          info.staNHits = nstaHits;
          TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(myTracks.probeVertPoint);
          info.staTransDCA = traj.perigeeParameters().transverseImpactParameter();
          info.staLongDCA = traj.perigeeParameters().longitudinalImpactParameter();
        }
      }
    }
    if (displacedStaTracks->size() > 0) {
      double disStaE = -10;
      for (reco::TrackCollection::const_iterator staTrack = displacedStaTracks->begin();
           staTrack != displacedStaTracks->end();
           ++staTrack) {
        reco::TransientTrack track = transientTrackBuilder->build(*staTrack);
        double dR = deltaR(track.impactPointTSCP().momentum().eta(),
                           track.impactPointTSCP().momentum().phi(),
                           selectedTrack->eta(),
                           selectedTrack->phi());
        if (dR < 1.0 && std::sqrt(track.impactPointTSCP().momentum().mag2() > disStaE)) {
          if (track.ndof() < 1) {
            continue;
          }
          if (track.chi2() / track.ndof() > 6) {
            continue;
          }
          info.disStadEta = track.impactPointTSCP().momentum().eta() - selectedTrack->eta();
          info.disStadPhi = track.impactPointTSCP().momentum().phi() - selectedTrack->phi();
          info.disStaChi2 = track.chi2() / track.ndof();

          info.disStaE = std::sqrt(track.impactPointTSCP().momentum().mag2());
          info.disStaDEoverE = (info.disStaE - mother.P()) / mother.P();
        }
      }
    }
    if (refittedStaTracks->size() > 0) {
      double disStaE = -10;
      for (reco::TrackCollection::const_iterator staTrack = refittedStaTracks->begin();
           staTrack != refittedStaTracks->end();
           ++staTrack) {
        reco::TransientTrack track = transientTrackBuilder->build(*staTrack);
        double dR = deltaR(track.impactPointTSCP().momentum().eta(),
                           track.impactPointTSCP().momentum().phi(),
                           selectedTrack->eta(),
                           selectedTrack->phi());
        if (dR < 1.0 && std::sqrt(track.impactPointTSCP().momentum().mag2() > disStaE)) {
          if (track.ndof() < 1) {
            continue;
          }
          if (track.chi2() / track.ndof() > 6) {
            continue;
          }
          info.refitStaE = std::sqrt(track.impactPointTSCP().momentum().mag2());
          info.refitStaDEoverE = (info.disStaE - mother.P()) / mother.P();
          info.refitStadEta = track.impactPointTSCP().momentum().eta() - selectedTrack->eta();
          info.refitStadPhi = track.impactPointTSCP().momentum().phi() - selectedTrack->phi();
          info.refitStaChi2 = track.chi2() / track.ndof();

          info.refitStaE = std::sqrt(track.impactPointTSCP().momentum().mag2());
          info.refitStaDEoverE = (info.disStaE - mother.P()) / mother.P();
        }
      }
    }

    info.closestApproach = closestapproach;
    if (staMinDr < 1.0) {
      info.standaloneE = staE;
      info.standaloneDEoverE = stadEoverE;

    }
    info.staMinDr = staMinDr;
    /*if (!MatchMotherToTrack(iEvent,
                                     tagVec,
                                     trackCollection_label,
                                     info)) {return;}*/
    if (!myHCAL.FindMuonHits(iEvent,
                             iSetup,
                             HBHERecHit_Label,
                             TrackGlobalPoint,
                             selectedTrack->charge(),
                             transientTrackBuilder->build(selectedTrack))&&!m_studyBarrel) {
      allEvents.FailLocation(myHCAL.ieta, myHCAL.iphi, weight_);
      return;
    }
    if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
    noBremEvents.IncCutFlow(info.eventWeight);
    if (genAccept) {
      genAccepted.IncCutFlow(info.eventWeight);
    }
    info.nCellsFound = myHCAL.CellsFound;
    info.nNeighborCellsFound = myHCAL.NeighborCellsFound;
    for (int depth = 0; depth < 7; depth++) {
      info.foundDepths[depth] = myHCAL.foundDepths[depth];
    }
    info.iEta = myHCAL.TrackiEta;
    info.cellCenterDEta = myHCAL.etaEdgeMin;
    info.cellCenterDPhi = myHCAL.phiEdgeMin;
    if (info.nCellsFound < 4&&!m_studyBarrel) {
      return;
    }
    //Check MET
    edm::Handle<std::vector<reco::PFMET>> pfMETCollection;
    bool ValidPFMET = iEvent.getByToken(m_PFMetToken, pfMETCollection);
    if (!ValidPFMET) {
      printf("Invalid PF Met!\n");
      return;
    }
    info.pfMet_pt = (*pfMETCollection)[0].et();
    double met = (*pfMETCollection)[0].et();
    double phi = (*pfMETCollection)[0].phi();
    TLorentzVector vmet, vlep;
    vlep.SetPxPyPzE(selectedMuon->px(), selectedMuon->py(), selectedMuon->pz(), selectedMuon->pfP4().E());
    double metx = met * std::cos(phi);
    double mety = met * std::sin(phi);
    vmet.SetPxPyPzE(metx, mety, 0.0, met);
    TLorentzVector vw = vlep + vmet;
    info.muMetMT = std::sqrt(2 * vlep.Et() * met * (1 - std::cos(deltaPhi(vlep.Phi(), phi))));
    vlep.SetPxPyPzE(selectedMuon->px() + selectedTrack->px(),
                    selectedMuon->py() + selectedTrack->py(),
                    selectedMuon->pz() + selectedTrack->pz(),
                    selectedMuon->pfP4().E() + selectedTrack->p());
    info.muAndProbeMetMT = std::sqrt(2 * vlep.Et() * met * (1 - std::cos(deltaPhi(vlep.Phi(), phi))));

    //Signal categories
    //bool passSigCuts = false;
    //if(myHCAL.m_HitsOverThresh<2){return;}
    if(!noBrem){allEvents.IncCutFlow(info.eventWeight);}
    noBremEvents.IncCutFlow(info.eventWeight);

    if (genAccept) {
      genAccepted.IncCutFlow(info.eventWeight);
    }
    //Generic HCAL hits in cone
    std::pair<double, int> coneHits = myHCAL.EnergyInCone(iEvent, iSetup, HBHERecHit_Label, TrackGlobalPoint);
    //Calo Jet Stuff
    double nearProbeCaloJetHadE = -1.;
    double nearProbeCaloJetEmE = -1.;
    double minProbeDr = -1;
    double caloJetPt = -1;
    for (auto CaloJet = caloJets->begin(); CaloJet != caloJets->end(); ++CaloJet) {
      double probeDr = deltaR(selectedTrack->eta(), selectedTrack->phi(), CaloJet->eta(), CaloJet->phi());
      double hadE = CaloJet->hadEnergyInHE() + CaloJet->hadEnergyInHB();
      double ecalE = CaloJet->emEnergyInEB() + CaloJet->emEnergyInEE();
      if (probeDr < minProbeDr || minProbeDr < 0) {
        minProbeDr = probeDr;
      }
      if (probeDr < 0.4) {
        if ((hadE + ecalE) > (nearProbeCaloJetHadE + nearProbeCaloJetEmE)) {
          nearProbeCaloJetHadE = hadE;
          nearProbeCaloJetEmE = ecalE;
        }
	caloJetPt = CaloJet->pt();
      }
    }
    for (int k = 0; k < 7; k++) {
      info.hitEnergies[k] = myHCAL.m_hitEnergies[k];
      info.neighborHitEnergies[k] = myHCAL.m_neighborHitEnergies[k];
    }
    bool onEdge = myHCAL.onEdge;
    info.bremDepth = myHCAL.m_bremDepth;
    info.minCaloJetDr = minProbeDr;
    info.caloJetPt = caloJetPt;
    info.hitsOverThresh = myHCAL.m_HitsOverThresh;
    info.coneHits = myHCAL.m_coneHits;
    info.coneEnergy = myHCAL.m_coneEnergy;
    info.caloJetHcalE = nearProbeCaloJetHadE;
    info.caloJetEcalE = nearProbeCaloJetEmE;
    info.caloJetTotalE = nearProbeCaloJetHadE + nearProbeCaloJetEmE;
    //info.caloConeE = coneHits.first;
    info.caloConeE = hcalIso;
    info.hcalIso0p3 = hcalIso0p3;
    info.hcalIso0p4 = hcalIso0p4;
    info.caloConeHits = coneHits.second;
    //PF Jets stuff
    myJets.SelectJets(iEvent, m_pfJetCollection_label);
    myJets.FindNearestJet(selectedTrack);
    myJets.SelectBTagJets(iEvent, m_jetTagTokens);
    info.nJets = myJets.nJets;
    info.nLooseBTagJets = myJets.looseBTagJets;
    info.nMedBTagJets = myJets.medBTagJets;
    info.nTightBTagJets = myJets.tightBTagJets;
    info.nearestJetDr = myJets.nearestJetDr;
    info.nearestJetE = myJets.nearestJetE;
    info.nearJetMaxE = myJets.nearJetMaxE;
    info.run = iEvent.id().run();
    info.lumi = iEvent.id().luminosityBlock();
    info.event = iEvent.id().event();

    if(lastSelectedTrack!=NULL)
    {
        CSC shiftCSCs(m_missingCSCFilter);
        HCAL shiftHCAL;

    	GlobalPoint LastTrackGlobalPoint =
    	    GlobalPoint(GlobalPoint::Polar(lastSelectedTrack->theta(), lastSelectedTrack->phi(), 525));

    	if (shiftHCAL.FindMuonHits(iEvent,
    	                         iSetup,
    	                         HBHERecHit_Label,
    	                         LastTrackGlobalPoint,
    	                         lastSelectedTrack->charge(),
    	                         transientTrackBuilder->build(lastSelectedTrack))) 
    	{
    	   shiftCSCs.ExtrapolateTrackToCSC(iEvent,
                                 iSetup,
                                 CSCSegment_Label,
                                 lastSelectedTrack,
                                 transientTrackBuilder->build(lastSelectedTrack),
                                 FillerPoint);
    	   info.shiftRPCDr = shiftCSCs.rpcMinDR;
	   info.shiftRecDr = shiftCSCs.recMinDR;
	   info.shiftCSCDr = shiftCSCs.minDR;

	   double shiftLastE = -1;
	   if(!shiftHCAL.onEdge)
	   {
		 for(int depth=0;depth<7;depth++)
		 {
		    if(shiftHCAL.foundDepths[depth])
		    {
		       shiftLastE = shiftHCAL.m_hitEnergies[depth];
		       info.shiftHE[depth] = shiftHCAL.m_hitEnergies[depth];
		    }
		 }
	   }
	   info.shiftLastE = shiftLastE;
	}
    } 
    lastSelObj = *selectedTrack;
    lastSelectedTrack = &lastSelObj;
    //Check nearby leptons
    for (const auto& iEle : iEvent.get(m_electronsToken)) {
      if (iEle.pt() < 15)
        continue;
      //if(!(iEle.electronID("cutBasedElectronID-Fall17-94X-V2-veto"))) continue;
      if (iEle.deltaPhiSuperClusterTrackAtVtx() < 0.8 && iEle.deltaEtaSuperClusterTrackAtVtx() < 0.01 &&
          iEle.sigmaIetaIeta() < 0.03 && iEle.hcalOverEcal() < 0.12 &&
          TMath::Abs(iEle.deltaEtaSuperClusterTrackAtVtx()) < 0.007 &&
          TMath::Abs(iEle.deltaPhiSuperClusterTrackAtVtx()) < 0.8) {
        info.nEle++;
        if (deltaR(selectedTrack->eta(), selectedTrack->phi(), iEle.eta(), iEle.phi()) < 0.2)
          continue;
        if (deltaR(selectedMuon->eta(), selectedMuon->phi(), iEle.eta(), iEle.phi()) < 0.2)
          continue;
        info.nExtraEle++;
      }
    }
    edm::Handle<std::vector<reco::Muon>> recoMuons;
    iEvent.getByToken(m_recoMuonToken, recoMuons);
    for (std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {
      if (!(iMuon->isPFMuon() && iMuon->isGlobalMuon()))
        continue;
      if (iMuon->pt() < 15)
        continue;
      if (deltaR(selectedTrack->eta(), selectedTrack->phi(), iMuon->eta(), iMuon->phi()) < 0.2)
        continue;
      if (deltaR(selectedMuon->eta(), selectedMuon->phi(), iMuon->eta(), iMuon->phi()) < 0.2)
        continue;
      info.nExtraMuon++;
    }
    if(noBrem)
    {
       noBremEvents.FillHists(info);
       return;
    }

    //Signal cut flow
    bool passMassCut = (info.diMuonMass > 80 && info.diMuonMass < 100);
    if (m_offPeak) {
      passMassCut = (info.diMuonMass < 70 || info.diMuonMass > 110);
    }
    if (m_ttbarCR) {
      passMassCut = (info.diMuonMass < 76 || info.diMuonMass > 106);
    }
    if (m_sameSign) {
      passMassCut = (info.diMuonMass > 50 && info.diMuonMass < 150);
    }
    if (myHCAL.m_coneEnergy < 10) {
      allEvents.IncCutFlow(info.eventWeight);
      if (genAccept) {
        genAccepted.IncCutFlow(info.eventWeight);
      }
      if (passMassCut) {
        allEvents.IncCutFlow(info.eventWeight);
        if (genAccept) {
          genAccepted.IncCutFlow(info.eventWeight);
        }
        if (info.cscMinDR > 0.05 && info.cscRecHitMinDR > 0.1) {
          allEvents.IncCutFlow(info.eventWeight);
          if (genAccept) {
            genAccepted.IncCutFlow(info.eventWeight);
          }
          if (staMinDr > 1.0) {
            allEvents.IncCutFlow(info.eventWeight);
            if (genAccept) {
              genAccepted.IncCutFlow(info.eventWeight);
            }
          }
        }
      }
    }
    bool passPreSelection = true;
    /*bool passSecCSCDrs = true;
    for(int station=1;station<4;station++)
    {
      if(myCSCs.minDrByDepth[station]<1.0){passSecCSCDrs=false;}
    }*/
    totalDisappear.IncCutFlow(info.eventWeight);
    if (passMassCut) {
      totalDisappear.IncCutFlow(info.eventWeight);
      if (staMinDr > 1.0) {
        totalDisappear.IncCutFlow(info.eventWeight);
        //if ((myCSCs.minDR > 0.05 || (m_missingCSCFilter && nGoodDepths == 1)) && info.cscRecHitMinDR > 0.1 && passSecCSCDrs) 
        if ((info.cscMinDR > 0.05 && info.cscRpcHitMinDR > 0.05 && nCloseHits<3) || (m_missingCSCFilter && nGoodDepths == 1)) 
	{
          if(m_missingCSCFilter&&nGoodDepths==1){info.eventWeight = info.eventWeight * (1 - cscScaleFactor);}

	  totalDisappear.IncCutFlow(info.eventWeight);
          if (info.hitsOverThresh > 1) {
            totalDisappear.IncCutFlow(info.eventWeight);
            if (info.hcalIso0p3<15 || m_highHcalCR) {
	      totalDisappear.IncCutFlow(info.eventWeight);
	      if(!onEdge){
		 double lastE = 10;
		 for(int depth=0;depth<7;depth++)
		 {
		    if(info.foundDepths[depth]){lastE = info.hitEnergies[depth];}
		 }
		 if(lastE<0.1)
		 {
                    totalDisappear.IncCutFlow(info.eventWeight);
                    info.sigCat = 3;
        	         //passSigCuts = true;
                    totalDisappear.FillHists(info);
		    totalDisappear.FillTotalDisappTree(info);
        	    if (info.probeEcalIso < 5) {
                       totalDisappearLowEcal.FillHists(info);
                    }
		 }
		 else{
                    failLast.FillHists(info);
		    failLast.FillTotalDisappTree(info);
		 }
	      }
            }
	    totalDisappearNoHcalIso.FillHists(info);
          }
        } //else if ((myCSCs.minDR > 0.05 || (m_missingCSCFilter && nGoodDepths == 1)) && passSecCSCDrs) 
	else
	{
          passSigNearbyCSC.IncCutFlow(info.eventWeight);
          if (true) {
	    bool steelBrem = true;
	    passSigNearbyCSC.IncCutFlow(info.eventWeight);
	    for(int cscStation = 1; cscStation<4; cscStation++){
	      if(fabs(info.minDrByStation[cscStation])<0.05){steelBrem=false;}
	      if(steelBrem){passSigNearbyCSC.IncCutFlow(info.eventWeight);}
	    }
            if (info.hitsOverThresh > 1) {
	      passSigNearbyCSC.IncCutFlow(info.eventWeight);
              if (info.hcalIso0p3<15 || m_highHcalCR) {
	        passSigNearbyCSC.IncCutFlow(info.eventWeight);
                passSigNearbyCSC.FillHists(info);
                passSigNearbyCSC.FillTotalDisappTree(info);
	      }
            }
	  }
        }
    }	
	if (staMinDr < 1.0) {
        partialDisappear.FillHists(info);
        partialDisappear.FillTree(info);
        if (info.probeEcalIso < 5) {
          partialDisappearLowEcal.FillHists(info);
          partialDisappearLowEcal.FillTree(info);
        }
        if (info.probeCharge > 0) {
          if (info.stadPhi > 0.1) {
            partialDisappear_largeCSCDefl.FillHists(info);
          }
        } else {
          if (info.stadPhi < -0.1) {
            partialDisappear_largeCSCDefl.FillHists(info);
          }
        }
        if (stadEoverE < -0.8 && myHCAL.m_coneEnergy < 1.) {
          info.sigCat = 2;
        } else if (stadEoverE < -0.6 && myHCAL.m_coneEnergy < 2.) {
          info.sigCat = 1;
        } else {
          info.sigCat = 0;
        }
      }
    } else {
      info.sigCat = 0;
    }
    //if (!passSigCuts) {
    //  failSig.FillHists(info);
   // }
    if (passPreSelection) {
      preSelectedEvents.FillHists(info);
    }
    if (info.hitsOverThresh < 2) {
      lowHits.FillHists(info);
    }
    allEvents.FillHists(info);
    if(!onEdge){allEvents.FillTotalDisappTree(info);}
    if(m_tagAligned){allEvents.FillTree(info);}
    if(m_doubleMuon){allEvents.FillTree(info);}

    if (genAccept) {
      genAccepted.FillHists(info);
    }
  }
}

bool SigAnalyzer::MatchMotherToTrack(const edm::Event& iEvent,
                                     const reco::Muon* mother,
                                     edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                                     SigEventInfo info) {
  edm::Handle<std::vector<reco::Track>> thePATTrackHandle;
  iEvent.getByToken(trackCollection_label, thePATTrackHandle);

  bool matched = false;
  double mindR;
  double mindE = 0.2;
  for (std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end();
       ++iTrack) {
    if (iTrack->pt() < 26)
      continue;
    double dR = deltaR(iTrack->eta(), iTrack->phi(), mother->eta(), mother->phi());
    double dEOverE = std::abs(iTrack->pt() - mother->pt()) / mother->pt();
    if (dR > 0.1)
      continue;
    if (dEOverE < mindE) {
      mindR = dR;
      mindE = dEOverE;
      selectedTagTrack = (&(*iTrack));
      matched = true;
      break;
    }
  }
  if (matched) {
    info.genProbeMatchdE = mindE;
    info.genProbeMatchdR = mindR;
  }

  return matched;
}

int SigAnalyzer::PairTagAndProbe(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    Tracks myTracks,
    const reco::Muon* Tag,
    ECAL myECAL,
    HCAL myHCAL,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
    edm::ESHandle<TransientTrackBuilder> transientTrackBuilder,
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
    edm::Handle<reco::VertexCollection> vtxHandle,
    edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label,
    edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label) {
  int cutpro = 0;
  for (std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin();
       iTrack != myTracks.selectedEndcapTracks.end();
       ++iTrack) {
    if (!m_sameSign) {
      if (Tag->charge() == (*iTrack)->charge())
        continue;
    } else {
      if (Tag->charge() != (*iTrack)->charge())
        continue;
    }
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
      double probeHcalIso0p2 =
          myHCAL.GetIsolation(iEvent, iSetup, HBHERecHit_Label, transientTrackBuilder->build(*iTrack), 0.2);
      double probeHcalIso0p3 = 
          myHCAL.GetIsolation(iEvent, iSetup, HBHERecHit_Label, transientTrackBuilder->build(*iTrack), 0.3);
      double probeHcalIso0p4 = 
          myHCAL.GetIsolation(iEvent, iSetup, HBHERecHit_Label, transientTrackBuilder->build(*iTrack), 0.4);
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
      if (probeEcalIso > 10. || probeEcalIso < 0.)
        continue;
      printf("Track passed iso, ecal iso is %f and eta/phi is %f, %f.\n",
             probeEcalIso,
             (*iTrack)->eta(),
             (*iTrack)->phi());
      if (cutpro < 6) {
        cutpro = 6;
      }
      //if(minProbeDr<0.4) continue;
      selectedTrack = (*iTrack);
      selectedMuon = Tag;
      hcalIso = probeHcalIso0p2;
      hcalIso0p3 = probeHcalIso0p3;
      hcalIso0p4 = probeHcalIso0p4;
      diMuonMass = myTracks.MuonTrackMass;
      pairVtxChi = myTracks.pairvertexchi;
      reco::TransientTrack tk = transientTrackBuilder->build(*iTrack);
      TrajectoryStateClosestToPoint traj = tk.trajectoryStateClosestToPoint(myTracks.probeVertPoint);
      probeTIP = traj.perigeeParameters().transverseImpactParameter();
      probeLIP = traj.perigeeParameters().longitudinalImpactParameter();
     //
      if (!m_highHcalCR) {
        if (probeHcalIso0p3 > 30.)
          continue;
      } else if (probeHcalIso0p2 < 30.) {
        continue;
      }
      if (cutpro < 7) {
        cutpro = 7;
      }
    }
  }
  return cutpro;
}

bool SigAnalyzer::PairMuons(const reco::Muon* Tag,
			   Muons myMuons,
			   Tracks myTracks,
                           edm::ESHandle<TransientTrackBuilder> transientTrackBuilder){
  bool paired = false;
  for (std::vector<const reco::Muon*>::const_iterator muon = myMuons.globalMuons.begin();
       muon != myMuons.globalMuons.end();
       muon++) {
    if (paired)
      continue;
    if (fabs((*muon)->eta()) > 2.4 || fabs((*muon)->eta()) < 1.45)
      continue;
    if (!(*muon)->isGlobalMuon())
      continue;
    if (deltaR((*muon)->eta(), (*muon)->phi(), Tag->eta(), Tag->phi()) < 0.2)
      continue;
    if ((*muon)->pt() < 20)
      continue;
    if (Tag->charge() == (*muon)->charge())
      continue;
    int pairingResult = myTracks.PairTracks((*muon)->globalTrack(), Tag->globalTrack(), transientTrackBuilder);
    if (pairingResult == 3) {
      if (myTracks.MuonTrackMass < 80 || myTracks.MuonTrackMass > 100){
	continue;
      }else{ 
	paired = true;
	selectedGlobalMuon = (*muon);
        selectedMuon = Tag;
        diMuonMass = myTracks.MuonTrackMass;
        pairVtxChi = myTracks.pairvertexchi;
      }
    }
  }
  return paired;
}



std::vector<int> SigAnalyzer::GetGenParticlesInCone(const edm::Event& iEvent,
						    const reco::Track* selTrack,
						    edm::EDGetToken m_genParticleToken){
  double minDR = 0.01;
  std::vector<int> particleIds;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);
  for (const auto& particle : *(genParticles.product())) {
    double dR = deltaR(particle.eta(), particle.phi(), selTrack->eta(), selTrack->phi());
    if (dR < minDR) {
       particleIds.push_back(particle.pdgId());
    }
  }
  return particleIds;
}

double SigAnalyzer::MatchTrackToGenMuon(const edm::Event& iEvent,
                                        const reco::Track* selTrack,
                                        edm::EDGetToken m_genParticleToken) {
  double minDR = -1;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(m_genParticleToken, genParticles);
  for (const auto& particle : *(genParticles.product())) {
    if (std::abs(particle.pdgId()) != 13)
      continue;
    double dR = deltaR(particle.eta(), particle.phi(), selTrack->eta(), selTrack->phi());
    if (dR < minDR || minDR < 0) {
      minDR = dR;
    }
  }
  return minDR;
}

bool SigAnalyzer::MatchTrackToMuon(const edm::Event& iEvent,
                                   math::XYZTLorentzVectorD mother,
                                   Muons myMuons,
                                   SigEventInfo info,
                                   double FmuE,
                                   double minCSCdr) {
  bool matched = false;
  double mindE = 1.;
  double globalmuonE;
  double mindR = 0.1;
  if (myMuons.selectedMuons.size() == 0) {
    return matched;
  }
  for (std::vector<const reco::Muon*>::const_iterator iMuon = myMuons.selectedMuons.begin();
       iMuon != myMuons.selectedMuons.end();
       ++iMuon) {
    //if(!(*iMuon)->isGlobalMuon()) {continue;}
    double dR = deltaR((*iMuon)->eta(), (*iMuon)->phi(), mother.eta(), mother.phi());
    double dEOverE = std::abs(std::sqrt(pow((*iMuon)->p(), 2) + pow(0.1056, 2)) - mother.E()) / mother.E();
    if (dR > 0.2)
      continue;
    if (dEOverE < mindE) {
      mindR = dR;
      mindE = dEOverE;
      matched = true;
      matchedMuon = (*iMuon);

      if ((*iMuon)->isAValidMuonTrack(reco::Muon::OuterTrack)) {
        standaloneE = (*iMuon)->outerTrack()->p();
      } else {
        standaloneE = -1.;
      }
      globalmuonE = (*iMuon)->p();
    }
  }
  if (matched) {
    info.probeMuondE = mindE;
    info.probeMuondR = mindR;
    info.standaloneE = standaloneE;
    info.globalMuonE = globalmuonE;
  } else {
    standaloneE = 0;
  }
  return matched;
}

std::vector<double> SigAnalyzer::GetPuWeight(edm::Handle<std::vector<PileupSummaryInfo>> hPileupInfoProduct) {
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
  float puUpWeight = float(puWeightUpHist_->GetBinContent(puWeightUpHist_->FindBin(nPUmean_)));
  if (nPUmean_ > 70) {
    puUpWeight = float(puWeightUpHist_->GetBinContent(puWeightUpHist_->FindBin(70)));
  }
  if (nPUmean_ < 1) {
    puUpWeight = float(puWeightUpHist_->GetBinContent(puWeightUpHist_->FindBin(0)));
  }
  float puDownWeight = float(puWeightDownHist_->GetBinContent(puWeightDownHist_->FindBin(nPUmean_)));
  if (nPUmean_ > 70) {
    puDownWeight = float(puWeightDownHist_->GetBinContent(puWeightDownHist_->FindBin(70)));
  }
  if (nPUmean_ < 1) {
    puDownWeight = float(puWeightDownHist_->GetBinContent(puWeightDownHist_->FindBin(0)));
  }

  std::vector<double> puWeights;
  puWeights.push_back(lNPVW);
  puWeights.push_back(puUpWeight);
  puWeights.push_back(puDownWeight);
  return puWeights;
}

bool SigAnalyzer::passesMetFilters(const edm::Event& iEvent) {
  int thisLumi = iEvent.id().luminosityBlock();
  int thisEvent = iEvent.id().event();
  if (metFailingEvents.count(thisLumi) > 0) {
    for (uint16_t event = 0; event < metFailingEvents[thisLumi].size(); event++) {
      if (thisEvent == metFailingEvents[thisLumi].at(event)) {
        return false;
      }
    }
  }
  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void SigAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SigAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SigAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SigAnalyzer);
