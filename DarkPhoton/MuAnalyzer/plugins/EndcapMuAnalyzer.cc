// -*- C++ -*-
//
// Package:    DarkPhoton/EndcapMuAnalyzer
// Class:      EndcapMuAnalyzer
// 
/**\class EndcapMuAnalyzer EndcapMuAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/EndcapMuAnalyzer.cc

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
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

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

class EndcapMuAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit EndcapMuAnalyzer(const edm::ParameterSet&);
      ~EndcapMuAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    bool isTrackMatchedToMuon(const edm::Event&, std::vector<const reco::Track*>::const_iterator&, MCHistograms);
    bool MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, SigEventInfo info);
    bool MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, SigEventInfo info, double FmuE, double minCSCdr);
    double MatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selTrack, edm::EDGetToken m_genParticleToken, double minGenDE);
    int PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, HCAL myHCAL, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label);
    std::vector<double> GetPuWeight(edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct);

    edm::EDGetToken m_recoMuonToken;
    edm::EDGetToken m_simTracksToken;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label;
    edm::EDGetToken m_genParticleToken;
    edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label;
//    edm::EDGetTokenT<CSCRecHit2DCollection> CSCRecHitToken;
    edm::EDGetTokenT<RPCRecHitCollection > rpcRecHitToken;
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
    edm::EDGetToken m_cscSimHits;
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

    edm::Handle<reco::CaloJetCollection> caloJets;
    MCHistograms allEvents;
    MCHistograms genAccepted;
    MCHistograms preSelectedEvents;
    MCHistograms lowHits;
    MCHistograms failSig;
    MCHistograms totalDisappear;
    MCHistograms partialDisappear;
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
EndcapMuAnalyzer::EndcapMuAnalyzer(const edm::ParameterSet& iConfig):
  m_recoMuonToken (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("recoMuons"))),
  trackCollection_label(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  primaryVertices_Label(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  CSCSegment_Label(consumes<CSCSegmentCollection > (iConfig.getParameter<edm::InputTag>("CSCSegmentLabel"))),
//  CSCRecHitToken(consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("CSCRecHitTag"))),
  rpcRecHitToken(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("RPCRecHitTag"))),
  HBHERecHit_Label(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("HBHERecHits"))),
  reducedEndcapRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"))),
  reducedBarrelRecHitCollection_Label(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"))),
  
  m_cscSimHits(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("cscSimHits"))),
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
  usesResource("TFileService");
  edm::Service<TFileService> fs;
   //now do what ever initialization is needed
  if (m_isMC){
    m_genParticleToken = consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"));
    m_simTracksToken = consumes<edm::SimTrackContainer> (iConfig.getParameter<edm::InputTag>("g4SimHits"));
    m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    m_PUInfoToken = consumes< std::vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("edmPileupInfo"));  
  }
  else
  {
     m_trigResultsToken = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("trigResults"));
     m_muonPathsToPass   = iConfig.getParameter<std::vector<std::string> >("muonPathsToPass");
  }
  allEvents.book(fs->mkdir("allEvents"));
  genAccepted.book(fs->mkdir("genAccepted"));
  preSelectedEvents.book(fs->mkdir("preSelected"));
  failSig.book(fs->mkdir("failSig"));
  lowHits.book(fs->mkdir("lowHits"));
  totalDisappear.book(fs->mkdir("totalDisappear"));
  partialDisappear.book(fs->mkdir("partialDisappear"));
  partialDisappear_largeCSCDefl.book(fs->mkdir("partialDisappear_largeCSCDefl"));
}


EndcapMuAnalyzer::~EndcapMuAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void EndcapMuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  Muons myMuons(true);
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
  bool genAccept = false;
  weight_=1.;
  info.eventWeight=weight_;
  info.genEventWeight=weight_;
  allEvents.IncCutFlow(info.eventWeight); 
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
  int maxTrackProgress=-1;
  iEvent.getByToken(m_caloJet_label, caloJets);
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
  myTracks.selectedTracks.clear();
  myTracks.selectedEndcapTracks.clear();
  int muProgress = myMuons.SelectMuons(iEvent,m_recoMuonToken);
  int trackProgress = myTracks.SelectUnpairedTracks(iEvent, trackCollection_label, vtxHandle);
  for(std::vector<const reco::Track*>::const_iterator iTrack = myTracks.selectedEndcapTracks.begin();iTrack != myTracks.selectedEndcapTracks.end();++iTrack)
  {
     selectedTrack = (*iTrack); 
     if(info.bremR==0){info.bremR=525;}
     //Find deflected muon or the initial muon if there was no deflection
     double FmuE = 0;
     math::XYZTLorentzVectorD mother;
     mother.SetPxPyPzE(selectedTrack->momentum().x(),selectedTrack->momentum().y(),selectedTrack->momentum().z(),selectedTrack->p());
     info.probeTrackIso = myTracks.GetIsolation(iEvent,trackCollection_label,0.3,vtxHandle,selectedTrack)/ selectedTrack->pt();
     info.isLoose = selectedTrack->quality(reco::Track::qualityByName("loose"));
     info.isTight = selectedTrack->quality(reco::Track::qualityByName("tight"));
     info.isHighPurity = selectedTrack->quality(reco::Track::qualityByName("highPurity"));
     info.probeReducedChi = selectedTrack->chi2()/ selectedTrack->ndof();
     info.pairVtxChi = pairVtxChi;
     info.probeEcalIso = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(selectedTrack));
     info.probeEcalIso0p1 = myECAL.GetIsolation(iEvent,iSetup,reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(selectedTrack),0.1);
     info.probeEcalIso0p2 = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(selectedTrack),0.2);
     info.probeEcalIso0p05 = myECAL.GetIsolation(iEvent,iSetup, reducedEndcapRecHitCollection_Label, reducedBarrelRecHitCollection_Label, transientTrackBuilder->build(selectedTrack),0.05);
     info.probeTrackEta = selectedTrack->eta();
     info.probeCharge = selectedTrack->charge();
     info.probeTrackPt = selectedTrack->pt();
     info.probeTrackPhi = selectedTrack->phi();
     GlobalPoint FillerPoint; 
     myCSCs.ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, rpcRecHitToken, selectedTrack,transientTrackBuilder->build(selectedTrack),  FillerPoint);
     info.cscRecHitMinDR=myCSCs.recMinDR;
     info.CSCdEta=myCSCs.CSCdEta;
     info.CSCdPhi=myCSCs.CSCdPhi;
     GlobalPoint TrackGlobalPoint = GlobalPoint(GlobalPoint::Polar(selectedTrack->theta(),selectedTrack->phi(),info.bremR));
     info.minCSCImpactParameter=myCSCs.minTotalImpactParameter;
     info.cscMinDR=myCSCs.minDR;
     if(FmuE==0){FmuE=selectedTrack->p();}
     bool MuonMatched = MatchTrackToMuon(iEvent,mother,myMuons,info,FmuE,myCSCs.minDR);
     double minGenDE;
     double minGenMuDr = MatchTrackToGenMuon(iEvent, selectedTrack, m_genParticleToken, minGenDE);
     if(minGenMuDr>0.2){continue;}
     info.minGenMuDr=minGenMuDr;
     info.minGenDE=minGenDE;
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
     edm::Handle<edm::PSimHitContainer> cscSimHits;
     iEvent.getByToken(m_cscSimHits,cscSimHits);
     edm::ESHandle<CSCGeometry> cscGeom;
     iSetup.get<MuonGeometryRecord>().get(cscGeom);
     info.minSimCSCDr=100;
     if(cscSimHits.isValid())
     {
        for(unsigned int simHitCounter= 0; simHitCounter < cscSimHits->size(); simHitCounter++)
        {
           const PSimHit& simHit = (*cscSimHits)[simHitCounter];
           if(abs(simHit.particleType()) == 13) 
           {
              const CSCLayer* csclayer = cscGeom->layer((CSCDetId)simHit.detUnitId());
              const auto& thisPoint = csclayer->toGlobal(simHit.localPosition());
              double SimCSCDr =  deltaR(thisPoint.eta(),thisPoint.phi(),selectedTrack->eta(),selectedTrack->phi());
              if(SimCSCDr<info.minSimCSCDr)
              {
                 info.minSimCSCDr=SimCSCDr;
                 info.minSimCSCELoss=simHit.energyLoss();
              }
           }
        }
     }
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
            if(track.ndof()>0)
            {
               if(track.chi2()/track.ndof()<10)
               {
                  if(std::sqrt(track.impactPointTSCP().momentum().mag2())>staE)
                  {
                     info.stadEta = track.impactPointTSCP().momentum().eta()-selectedTrack->eta();
                     info.stadPhi = track.impactPointTSCP().momentum().phi()-selectedTrack->phi();
                  }
               }
            }
          }
          //Pick the largest energy standalone muon within dR of 1.0
          if(dR<1.0&&std::sqrt(track.impactPointTSCP().momentum().mag2())>staE)
          {
            if(track.ndof()<1){continue;}
            if(track.chi2()/track.ndof()>10){continue;}
            staE = std::sqrt(track.impactPointTSCP().momentum().mag2());
            info.stadEta = track.impactPointTSCP().momentum().eta()-selectedTrack->eta();
            info.stadPhi = track.impactPointTSCP().momentum().phi()-selectedTrack->phi();
            stadEoverE=(std::sqrt(track.impactPointTSCP().momentum().mag2())-mother.P())/mother.P();
            info.staChi2=track.chi2()/track.ndof();
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

     info.closestApproach=closestapproach;
     if(staMinDr<1.0) 
     {
        info.standaloneE = staE;
        info.standaloneDEoverE = stadEoverE;
     }
     info.staMinDr=staMinDr;
     if(!myHCAL.FindMuonHits(iEvent, iSetup, HBHERecHit_Label, TrackGlobalPoint, selectedTrack->charge(), transientTrackBuilder->build(selectedTrack)))
     {
        allEvents.FailLocation(myHCAL.ieta,myHCAL.iphi, weight_);
        continue;
     }
     allEvents.IncCutFlow(info.eventWeight);
     info.nCellsFound=myHCAL.CellsFound;
     info.nNeighborCellsFound=myHCAL.NeighborCellsFound;
     info.iEta=myHCAL.TrackiEta;
     info.cellCenterDEta=myHCAL.Tdeta;
     info.cellCenterDPhi=myHCAL.Tdphi;
     if(info.nCellsFound<4){continue;}
     //Signal categories
     bool passSigCuts = false;
     //if(myHCAL.m_HitsOverThresh<2){return;}
     allEvents.IncCutFlow(info.eventWeight);
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
        info.neighborHitEnergies[k]=myHCAL.m_neighborHitEnergies[k];
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
     if(info.probeEcalIso<0){continue;}
     if(staMinDr>1.0&&myCSCs.minDR>0.05)
     {
       if(myHCAL.m_coneEnergy<10&&info.caloConeE<20&&info.probeEcalIso<5&&info.caloJetTotalE<0)
       {
          info.sigCat=3;
          passSigCuts = true;
          totalDisappear.FillHists(info);
       }
     }
     else if(staMinDr<1.0)
     {
        partialDisappear.FillHists(info);
        if(info.probeCharge>0)
        {
           if(info.stadPhi>0.1){partialDisappear_largeCSCDefl.FillHists(info);}
        }
        else
        {
           if(info.stadPhi<-0.1){partialDisappear_largeCSCDefl.FillHists(info);}
        }
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
     if(info.hitsOverThresh<2){lowHits.FillHists(info);}
     allEvents.FillHists(info);
  }
}

bool EndcapMuAnalyzer::MatchMotherToTrack(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, SigEventInfo info)
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

int EndcapMuAnalyzer::PairTagAndProbe(const edm::Event& iEvent, const edm::EventSetup& iSetup, Tracks myTracks, const reco::Muon* Tag, ECAL myECAL, HCAL myHCAL, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, edm::ESHandle<TransientTrackBuilder> transientTrackBuilder, edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label, edm::Handle<reco::VertexCollection> vtxHandle, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label)
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
       double probeHcalIso = myHCAL.GetIsolation(iEvent,iSetup, HBHERecHit_Label, transientTrackBuilder->build(*iTrack));
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
       if(probeEcalIso>5.||probeEcalIso<0.) continue;
       if(cutpro<6){cutpro=6;}
       //if(minProbeDr<0.4) continue;
       if(probeHcalIso>20.) continue;
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

double EndcapMuAnalyzer::MatchTrackToGenMuon(const edm::Event& iEvent, const reco::Track* selTrack, edm::EDGetToken m_genParticleToken, double minGenDE)
{
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(m_genParticleToken, genParticles);
   double minMuDr=1000;
   for (const auto& particle : *(genParticles.product()))
   {
      if (std::abs(particle.pdgId())!= 13)
         continue;
      double dR = deltaR(particle.eta(), particle.phi(), selTrack->eta(), selTrack->phi());
      if(dR<minMuDr)
      {
         minMuDr=dR;
         minGenDE=(selTrack->pt()-particle.pt())/particle.pt(); 
      }
   }
   return minMuDr;
}


bool EndcapMuAnalyzer::MatchTrackToMuon(const edm::Event& iEvent, math::XYZTLorentzVectorD mother, Muons myMuons, SigEventInfo info, double FmuE, double minCSCdr)
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

std::vector<double> EndcapMuAnalyzer::GetPuWeight(edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct)
{
   const std::vector<PileupSummaryInfo> *inPUInfos = hPileupInfoProduct.product();
   for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = inPUInfos->begin(); itPUInfo!=inPUInfos->end(); ++itPUInfo)
   {
      if(itPUInfo->getBunchCrossing()==0){nPUmean_  = itPUInfo->getTrueNumInteractions();}
   }
   float lNPVW = float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(nPUmean_)));
   if(nPUmean_>70){lNPVW=float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(70)));}
   if(nPUmean_<1){lNPVW=float(puWeightHist_->GetBinContent(puWeightHist_->FindBin(0)));}
   float puUpWeight = float(puWeightUpHist_->GetBinContent(puWeightUpHist_->FindBin(nPUmean_)));
   if(nPUmean_>70){puUpWeight=float(puWeightUpHist_->GetBinContent(puWeightUpHist_->FindBin(70)));}
   if(nPUmean_<1){puUpWeight=float(puWeightUpHist_->GetBinContent(puWeightUpHist_->FindBin(0)));}
   float puDownWeight = float(puWeightDownHist_->GetBinContent(puWeightDownHist_->FindBin(nPUmean_)));
   if(nPUmean_>70){puDownWeight=float(puWeightDownHist_->GetBinContent(puWeightDownHist_->FindBin(70)));}
   if(nPUmean_<1){puDownWeight=float(puWeightDownHist_->GetBinContent(puWeightDownHist_->FindBin(0)));}
 
   std::vector<double> puWeights;
   puWeights.push_back(lNPVW);
   puWeights.push_back(puUpWeight);
   puWeights.push_back(puDownWeight);
   return puWeights;
}

// ------------ method called once each job just before starting event loop  ------------
void 
EndcapMuAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EndcapMuAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EndcapMuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EndcapMuAnalyzer);
