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
#include "DataFormats/PatCandidates/interface/Electron.h"
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

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SigAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
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
    bool MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, SigEventInfo info);
    bool MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, SigEventInfo info, double FmuE, double minCSCdr);
    int PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label);
    double GetPuWeight(edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct);

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
    edm::EDGetTokenT<reco::GsfElectronCollection> m_electronsToken;  
    edm::EDGetToken m_PUInfoToken;
    edm::EDGetToken m_theSTAMuonLabel;
    edm::EDGetToken m_displacedSTAMuonLabel;
    edm::EDGetToken m_pfJetCollection_label;
    bool m_isMC;
    bool m_isSig;
    bool m_hasDpho;
    double m_scaleFactor;
    bool m_runRandomTrackEfficiency;
    double standaloneE;
    double weight_;
    double nPUmean_;
    double pairVtxChi;
    double diMuonMass;
    double probeTIP;
    double probeLIP;
   
    //pileup reweighting
    TH1F* fPUDataHist_;
    TH1F* fPUMCHist_;
    TH1F* puWeightHist_;

    edm::Handle<reco::CaloJetCollection> caloJets;
    MCHistograms allEvents;
    MCHistograms genAccepted;
    MCHistograms preSelectedEvents;
    MCHistograms fourHits;
    MCHistograms failSig;
    MCHistograms totalDisappear;
    MCHistograms partialDisappear;
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
SigAnalyzer::SigAnalyzer(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  reducedEndcapRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
  reducedBarrelRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
  m_caloJet_label(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("CaloJetSource"))),
  m_electronsToken(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  m_theSTAMuonLabel(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("StandAloneTracks"))),
  m_displacedSTAMuonLabel(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("DisplacedStandaloneTracks"))),
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
     m_PUInfoToken = consumes< std::vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("edmPileupInfo"));  
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

//  TFile *puFile = TFile::Open("/hdfs/cms/user/revering/dphoton/PileupHists/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root");
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

  allEvents.book(fs->mkdir("allEvents"));
  genAccepted.book(fs->mkdir("genAccepted"));
  preSelectedEvents.book(fs->mkdir("preSelected"));
  failSig.book(fs->mkdir("failSig"));
  fourHits.book(fs->mkdir("fourHits"));
  totalDisappear.book(fs->mkdir("totalDisappear"));
  partialDisappear.book(fs->mkdir("partialDisappear"));
}


SigAnalyzer::~SigAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void SigAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  genAccepted.ResetCutFlow();
  edm::Handle<SimTrackContainer> simtracks;
  edm::Handle<float> dbremweight, dbremInitialE, dbremFinalE, dbremAngle, dbremX, dbremY, dbremZ;
  bool genAccept = false;
  if(m_hasDpho)
  {
     iEvent.getByToken(m_dbremweighttoken,dbremweight);
     if(*dbremweight>0)
     {
       iEvent.getByToken(m_dbremFinalEToken,dbremFinalE);
       iEvent.getByToken(m_dbremInitialEToken,dbremInitialE);
       iEvent.getByToken(m_dbremAngleToken,dbremAngle);
       iEvent.getByToken(m_dbremXToken,dbremX);
       iEvent.getByToken(m_dbremYToken,dbremY);
       iEvent.getByToken(m_dbremZToken,dbremZ);
       weight_ = *dbremweight*m_scaleFactor;
       if(weight_>10){return;}
       edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct;
       iEvent.getByToken(m_PUInfoToken,hPileupInfoProduct); 
       assert(hPileupInfoProduct.isValid()); 
       double PUWeight = GetPuWeight(hPileupInfoProduct);
       weight_ = weight_*PUWeight;
       info.nPUmean=nPUmean_;
       info.pileupWeight=PUWeight;
       if(*dbremInitialE<20){return;}
       edm::Handle<reco::GenParticleCollection> genParticles;
       iEvent.getByToken(m_genParticleToken, genParticles);
       bool maybeTag = false;
       bool maybeProbe = false; 
       for (const auto& particle : *(genParticles.product()))
       {
          if (std::abs(particle.pdgId()) != 13){continue;}
          if (particle.mother()->pdgId() != 23){continue;} 
          if (fabs(particle.eta())>2.4){continue;}
          if (particle.pt()<20){continue;}
          if(maybeProbe&&particle.pt()>26)
          {
             genAccept = true;  
             break;
          }
          if (fabs(particle.eta()>1.6))
          {
             maybeProbe=true;
             if(maybeTag)
             { 
                genAccept=true;
                break;
             }
          }
          if (particle.pt()>26)
          {
            maybeTag=true;
          }
       }

     }
     else{return;}
  }
  else{weight_ = 1.;}
  info.eventWeight=weight_;
  allEvents.IncCutFlow(info.eventWeight); 
  if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
  edm::Handle<reco::TrackCollection> staTracks;
  iEvent.getByToken(m_theSTAMuonLabel, staTracks);
  edm::Handle<reco::TrackCollection> displacedStaTracks;
  iEvent.getByToken(m_displacedSTAMuonLabel, displacedStaTracks);
  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(primaryVertices_Label, vtxHandle);
  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  int muProgress = myMuons.SelectMuons(iEvent, m_recoMuonToken);
  for(int i=0;i<muProgress;i++)
  {
     allEvents.IncCutFlow(info.eventWeight);
     if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
  }
  if(muProgress<4){return;}
  bool SiblingMatched = false;
  int maxTrackProgress=0;
  iEvent.getByToken(m_caloJet_label, caloJets);

  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  for(std::vector<const reco::Muon*>::const_iterator muon = myMuons.selectedMuons.begin(); muon!=myMuons.selectedMuons.end();muon++)
  { 
     if(maxTrackProgress==13){continue;}
     myTracks.selectedTracks.clear();
     myTracks.selectedEndcapTracks.clear();
     int trackProgress = myTracks.SelectTracks(iEvent, trackCollection_label, vtxHandle, *muon, transientTrackBuilder);
     if(trackProgress==6){trackProgress=trackProgress+PairTagAndProbe(iEvent,iSetup,myTracks,*muon,myECAL,transientTrackBuilder,trackCollection_label, vtxHandle,reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label);}
     if(trackProgress>maxTrackProgress){maxTrackProgress=trackProgress;}
  }
  if(maxTrackProgress==13)
  {
       SiblingMatched=true;
  }
  for(int i=0;i<maxTrackProgress;i++)
  {
     allEvents.IncCutFlow(info.eventWeight);
     if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
  }
  if(maxTrackProgress<13){return;}
  math::XYZTLorentzVectorD mother;
  //find dark brem vertex
  if(m_isMC&&m_hasDpho)
  {
     GlobalPoint vtx(*dbremX/10.,*dbremY/10.,*dbremZ/10.);
//     if(fabs(vtx.z())>388&&vtx.perp()<268&&vtx.perp()<(1.2857*fabs(vtx.z()-302.8)))
     info.bremR = vtx.mag();
     info.vtxZ = vtx.z();
     info.vtxPerp = vtx.perp();
     info.dbremInitialE=*dbremInitialE;
  //   }
  }
  if(info.bremR==0){info.bremR=525;}
  //Find deflected muon or the initial muon if there was no deflection
  double FmuE = 0;
  if(m_isMC&&m_hasDpho)
  {
     info.finalMuE=*dbremFinalE;
     FmuE=*dbremFinalE;
  }

  info.nTagMuons = myMuons.selectedMuons.size();
  if(SiblingMatched)
  {
     for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedTracks.begin(); iTrack != myTracks.selectedTracks.end(); ++iTrack)
     {
        if(deltaR((*iTrack)->eta(),(*iTrack)->phi(),selectedTrack->eta(),selectedTrack->phi())<0.3) continue;
        if(selectedMuon->charge()==(*iTrack)->charge()) continue;
        //if(fabs((*iTrack)->dxy())>0.05) continue;
        if(myTracks.PairTracks(iTrack, selectedMuon->globalTrack(), transientTrackBuilder))
        {       
          double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label,0.3,vtxHandle,*iTrack);
          double probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
          if((probeTrackIso/(*iTrack)->pt())>0.2||probeEcalIso>30.) continue;
          if(probeTrackIso<0){continue;}
          SiblingMatched=false;
        } 
     }
    mother.SetPxPyPzE(selectedTrack->momentum().x(),selectedTrack->momentum().y(),selectedTrack->momentum().z(),selectedTrack->p());
  }
//  if(!myEventInfo.goodPrimaryVertex(iEvent, primaryVertices_Label)) return;
  if(SiblingMatched)
  {
     //Accept events that pass preselection w/o passing gen acceptance
     if(!genAccept)
     {
        genAccept=true;
        for(int k=0;k<16;k++){genAccepted.IncCutFlow(info.eventWeight);}
     }
     allEvents.IncCutFlow(info.eventWeight); 
     if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
     info.probeTrackIso = myTracks.GetIsolation(iEvent,trackCollection_label,0.3,vtxHandle,selectedTrack)/selectedTrack->pt();
     info.isLoose = selectedTrack->quality(reco::Track::qualityByName("loose"));
     info.isTight = selectedTrack->quality(reco::Track::qualityByName("tight"));
     info.isHighPurity = selectedTrack->quality(reco::Track::qualityByName("highPurity"));
     info.probeReducedChi = selectedTrack->chi2()/selectedTrack->ndof();
     info.tagMuonEta = selectedMuon->eta();
     info.pairVtxChi = pairVtxChi;
     info.probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(selectedTrack));
     info.tagMuonPt = selectedMuon->pt();
     info.probeTrackEta = selectedTrack->eta();
     info.probeTrackPt = selectedTrack->pt();
     info.probeTrackPhi = selectedTrack->phi();
     info.dxy = selectedTrack->dxy();
     info.dca = probeTIP;
     info.dcal = probeLIP;
     info.dsz = selectedTrack->dsz();
     info.diMuonMass = diMuonMass;
   // this wraps tracks with additional methods that are used in vertex-calculation
     GlobalPoint FillerPoint; 
     myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, selectedTrack,transientTrackBuilder->build(selectedTrack),  FillerPoint);
     GlobalPoint TrackGlobalPoint = GlobalPoint(GlobalPoint::Polar(selectedTrack->theta(),selectedTrack->phi(),info.bremR));
     if(m_isMC&&m_hasDpho)
     {
        info.deflectingAngle=*dbremAngle;
        info.dbremFinalE=*dbremFinalE;
     }
     info.minCSCImpactParameter=myCSCs.minTotalImpactParameter;
     info.cscMinDR=myCSCs.minDR;
     if(FmuE==0){FmuE=selectedTrack->p();}
     bool MuonMatched = MatchTrackToMuon(iEvent,mother,myMuons,info,FmuE,myCSCs.minDR);
     double staMinDr = 7.;
     if(MuonMatched)
     {
        info.cscChiSq=myCSCs.CSCChiSq;
        ROOT::Math::DisplacementVector3D unitmother = mother.Vect().unit();
	GlobalVector unitcsc = myCSCs.CSCTraj.unit();
        double dotpro = unitmother.x()*unitcsc.x()+unitmother.y()*unitcsc.y()+unitmother.z()*unitcsc.z();
        info.cscHitAngleChange=std::acos(dotpro);
     }    
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
            if(track.ndof()<1){continue;}
            if(track.chi2()/track.ndof()>10){continue;}
            staE = std::sqrt(track.impactPointTSCP().momentum().mag2());
            stadEoverE=(std::sqrt(track.impactPointTSCP().momentum().mag2())-mother.P())/mother.P();
            if(track.ndof()>0){info.staChi2=track.chi2()/track.ndof();}
            int nstaHits=0;
            trackingRecHit_iterator rhbegin = staTrack->recHitsBegin();
            trackingRecHit_iterator rhend = staTrack->recHitsEnd();
            for (trackingRecHit_iterator recHit = rhbegin; recHit != rhend; ++recHit) {
               nstaHits++;
            }
            info.staNHits=nstaHits;
            TrajectoryStateClosestToPoint  traj = track.trajectoryStateClosestToPoint(myTracks.probeVertPoint);
            info.staTransDCA=traj.perigeeParameters().transverseImpactParameter();
            info.staLongDCA=traj.perigeeParameters().longitudinalImpactParameter();
          }
       }
     } 
     if(displacedStaTracks->size()>0)
     {
       double disStaE = -10;
       for(reco::TrackCollection::const_iterator staTrack = displacedStaTracks->begin(); staTrack != displacedStaTracks->end(); ++staTrack)
       {
          reco::TransientTrack track = transientTrackBuilder->build(*staTrack); 
          double dR = deltaR(track.impactPointTSCP().momentum().eta(),track.impactPointTSCP().momentum().phi(),selectedTrack->eta(),selectedTrack->phi());
          if(dR<0.5&&std::sqrt(track.impactPointTSCP().momentum().mag2()>disStaE))
          {
            if(track.ndof()<1){continue;}
            if(track.chi2()/track.ndof()>6){continue;}
            info.disStaE = std::sqrt(track.impactPointTSCP().momentum().mag2());
            info.disStaDEoverE = (info.disStaE-mother.P())/mother.P();
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
     if(!myHCAL.FindMuonHits(iEvent, iSetup, HBHERecHit_Label, TrackGlobalPoint, info, info.standaloneE, weight_, selectedTrack->dsz(), selectedTrack->charge(), transientTrackBuilder->build(selectedTrack)))
     {
        allEvents.FailLocation(myHCAL.ieta,myHCAL.iphi, weight_);
        return;
     }
     allEvents.IncCutFlow(info.eventWeight);
     if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
     info.nCellsFound=myHCAL.CellsFound;
     if(info.nCellsFound<4){return;}
     //Signal categories
     bool passSigCuts = false;
     //if(myHCAL.m_HitsOverThresh<2){return;}
     allEvents.IncCutFlow(info.eventWeight);
     if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
     //Generic HCAL hits in cone
     std::pair<double, int> coneHits = myHCAL.EnergyInCone(iEvent,iSetup,HBHERecHit_Label, TrackGlobalPoint);
     //Calo Jet Stuff
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
     info.bremDepth = myHCAL.m_bremDepth;
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
     //Check nearby leptons
     for (const auto& iEle : iEvent.get(m_electronsToken)) 
     {
        if(iEle.pt()<15) continue;
        //if(!(iEle.electronID("cutBasedElectronID-Fall17-94X-V2-veto"))) continue;
        if(iEle.deltaPhiSuperClusterTrackAtVtx() < 0.8 && iEle.deltaEtaSuperClusterTrackAtVtx() < 0.01 && iEle.sigmaIetaIeta() < 0.03 && iEle.hcalOverEcal() < 0.12 && TMath::Abs(iEle.deltaEtaSuperClusterTrackAtVtx()) < 0.007 && TMath::Abs(iEle.deltaPhiSuperClusterTrackAtVtx()) < 0.8)
        {
           info.nEle++;
           if(deltaR(selectedTrack->eta(),selectedTrack->phi(),iEle.eta(),iEle.phi())<0.2) continue;
           if(deltaR(selectedMuon->eta(),selectedMuon->phi(),iEle.eta(),iEle.phi())<0.2) continue;
           info.nExtraEle++;
        }
     }                   
     edm::Handle<std::vector<reco::Muon>> recoMuons;
     iEvent.getByToken(m_recoMuonToken, recoMuons);
     for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end();iMuon++)
     {
        if(!(iMuon->isPFMuon() && iMuon->isGlobalMuon())) continue;
        if(iMuon->pt()<15) continue;
        if(deltaR(selectedTrack->eta(),selectedTrack->phi(),iMuon->eta(),iMuon->phi())<0.2) continue;
        if(deltaR(selectedMuon->eta(),selectedMuon->phi(),iMuon->eta(),iMuon->phi())<0.2) continue;
        info.nExtraMuon++;
     }                                      
     //Signal cut flow

     if(myHCAL.m_coneEnergy<10)
     {
        allEvents.IncCutFlow(info.eventWeight);
        if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
        if(myCSCs.minDR>0.05)
        {
           allEvents.IncCutFlow(info.eventWeight);
           if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
           if(staMinDr>1.0)
           {
              allEvents.IncCutFlow(info.eventWeight);
              if(genAccept){genAccepted.IncCutFlow(info.eventWeight);}
           }
        }
     }
     bool passPreSelection=false;
     passPreSelection=true;
     if(staMinDr>1.0&&myCSCs.minDR>0.05)
     {
       if(myHCAL.m_coneEnergy<10){info.sigCat=3;}
       passSigCuts = true;
       totalDisappear.FillHists(info);
     }
     else if(staMinDr<1.0)
     {
        partialDisappear.FillHists(info);
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
     if(!passSigCuts)
     {
        failSig.FillHists(info);
     }
     if(passPreSelection)
     {
        preSelectedEvents.FillHists(info);
     }
     if(info.coneHits==4){fourHits.FillHists(info);}
     allEvents.FillHists(info);
     if(genAccept){genAccepted.FillHists(info);}
  }
}

bool SigAnalyzer::MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, SigEventInfo info)
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

int SigAnalyzer::PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label)
{
  int cutpro = 0;
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin(); iTrack != myTracks.selectedEndcapTracks.end(); ++iTrack)
  {
     if(Tag->charge()==(*iTrack)->charge()) continue;
     if(cutpro<1){cutpro=1;}
     int pairingResult = myTracks.PairTracks(iTrack, Tag->globalTrack(), transientTrackBuilder);
     if(pairingResult>(cutpro-1)){cutpro=pairingResult+1;}
     if(pairingResult==3)
     {
       double probeTrackIso = myTracks.GetIsolation(iEvent, trackCollection_label,0.3,vtxHandle,*iTrack);
       double probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(*iTrack));
       double minProbeDr=5;
       for(auto CaloJet = caloJets->begin(); CaloJet!= caloJets->end(); ++CaloJet)
       {
          double probeDr = deltaR((*iTrack)->eta(),(*iTrack)->phi(), CaloJet->eta(), CaloJet->phi());
          double hadE = CaloJet->hadEnergyInHE()+CaloJet->hadEnergyInHB();
          double ecalE = CaloJet->emEnergyInEB()+CaloJet->emEnergyInEE();
          if(probeDr<minProbeDr)
          {
             if((hadE+ecalE)>15){minProbeDr=probeDr;}
          }
       }

       if((probeTrackIso/(*iTrack)->pt())>0.05) continue;
       if(cutpro<5){cutpro=5;}
       if(probeEcalIso>10.) continue;
       if(cutpro<6){cutpro=6;}
       if(minProbeDr<0.4) continue;
       if(cutpro<7){cutpro=7;}
       selectedTrack = (*iTrack); 
       selectedMuon = Tag;
       diMuonMass = myTracks.MuonTrackMass;
       pairVtxChi = myTracks.pairvertexchi;
       reco::TransientTrack tk = transientTrackBuilder->build(*iTrack);
       TrajectoryStateClosestToPoint  traj = tk.trajectoryStateClosestToPoint(myTracks.probeVertPoint);
       probeTIP = traj.perigeeParameters().transverseImpactParameter();
       probeLIP = traj.perigeeParameters().longitudinalImpactParameter();
     } 
  }
  return cutpro;
}

bool SigAnalyzer::MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, SigEventInfo info, double FmuE, double minCSCdr)
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

double SigAnalyzer::GetPuWeight(edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct)
{
   const std::vector<PileupSummaryInfo> *inPUInfos = hPileupInfoProduct.product();
   for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = inPUInfos->begin(); itPUInfo!=inPUInfos->end(); ++itPUInfo)
   {
      if(itPUInfo->getBunchCrossing()==0)
      nPUmean_  = itPUInfo->getTrueNumInteractions();
   }
   float lNPVW = float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(nPUmean_)));
   if(nPUmean_>70){lNPVW=float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(70)));}
   if(nPUmean_<1){lNPVW=float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(0)));}

   return lNPVW;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SigAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SigAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SigAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SigAnalyzer);
