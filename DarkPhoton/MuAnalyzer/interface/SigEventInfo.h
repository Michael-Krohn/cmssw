#ifndef SigEventInfo_h
#define SigEventInfo_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TH1F.h"

class SigEventInfo {
public:
  SigEventInfo();
  bool passTriggers(const edm::Event& iEvent,
                    edm::EDGetToken m_trigResultsToken,
                    edm::EDGetToken m_trigEventToken,
                    std::vector<std::string> m_muonPathsToPass);
  bool goodPrimaryVertex(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label);
  bool matchTagAndTrigger(const reco::Muon* selectedMuon);
  double eventWeight;
  double dbremWeight;
  double genEventWeight;
  bool filledWeights;
  std::vector<double> PuUpDownWeights;
  double muRocWeight;
  double ZptWeight;
  double Zpt;
  bool foundZ;
  std::vector<double> muIdWeights;
  std::vector<double> muIsoWeights;
  std::vector<double> muTrigWeights;
  int sigCat;
  double dbremZ;
  double dbremPerp;
  double dbremInitialEt;
  double bremR;
  double dbremInitialE;
  double dbremFinalE;
  double deflectingAngle;
  int nTagMuons;
  int nPairedProbes;
  double probeTrackIso;
  double probeEcalIso;
  double probeEcalIso0p1;
  double probeEcalIso0p2;
  double probeEcalIso0p05;
  bool isLoose;
  bool isTight;
  bool isHighPurity;
  double probeReducedChi;
  double tagMuonEta;
  double tagMuonPt;
  double pfMet_pt;
  double muMetMT;
  double muAndProbeMetMT;
  double probeTrackEta;
  double probeTrackPhi;
  double probeTrackPt;
  double probeTrackP;
  double pairVtxChi;
  double diMuonMass;
  double altInvMass;
  double minSimCSCDr;
  double minSimCSCELoss;

  //CSC parameters
  double minCSCImpactParameter;
  double cscMinDR;
  int nSimCSCHits;
  int nCloseCSCHits[4];
  double minDrByStation[4];
  double secMinDrByStation[4];
  double minDtDrByStation[4];
  double minDEtaByStation[4];
  double minDtDEtaByStation[4];
  double minDPhiByStation[4];
  double secMinDPhiByStation[4];
  double minDtDPhiByStation[4];
  double minDZByStation[4];
  double minDtDZByStation[4];
  double cscRecHitMinDR;
  double cscRpcHitMinDR;
  double cscChiSq;
  double cscHitAngleChange;
  //Transient Track Studies
  double closestApproach;
  double staMinDr;
  double standaloneDEoverE;
  double dxy;
  double dsz;
  double dca;
  double dcal;
  //HCAL
  int nCellsFound;
  int nNeighborCellsFound;
  int iEta;
  double cellCenterDEta;
  double cellCenterDPhi;
  double coneEnergy;
  int coneHits;
  double hitEnergies[7];
  bool foundDepths[7];
  double neighborHitEnergies[7];
  int hitsOverThresh;
  double vtxZ;
  double vtxPerp;
  int bremDepth;
  //Pf Jets
  double nearestJetDr;
  double nearJetMaxE;
  double nearestJetE;
  int nJets;
  int nLooseBTagJets;
  int nMedBTagJets;
  int nTightBTagJets;
  //Calo Jets
  double minCaloJetDr;
  double caloJetPt;
  double caloJetHcalE;
  double caloJetEcalE;
  double caloJetTotalE;
  double caloConeE;
  double hcalIso0p3;
  double hcalIso0p4;
  int caloConeHits;
  //Gen Matching
  double minGenMuDr;
  double minGenDE;
  double genProbeMatchdE;
  double genProbeMatchdR;
  double probeMuondE;
  double probeMuondR;
  double standaloneE;
  double globalMuonE;
  int particlePDG;
  double minMuDr;
  int nEle;
  int nExtraEle;
  int nExtraMuon;
  int probeCharge;
  double CSCdEta;
  double CSCdPhi;
  double stadEta;
  double stadPhi;
  double staChi2;
  int staNHits;
  double staTransDCA;
  double staLongDCA;
  double finalMuE;
  double disStaE;
  double disStaDEoverE;
  double disStadEta;
  double disStadPhi;
  double disStaChi2;
  double refitStaE;
  double refitStaDEoverE;
  double refitStadEta;
  double refitStadPhi;
  double refitStaChi2;


  int lumi;
  int event;
  int run;
  //pileup
  double pileupWeight;
  int nPUmean;
  //Energy bin SFs
  int usedBin;
  double energyBinSF;
  std::vector<int> pdgIds;
  //Trigger objects
  std::vector<trigger::TriggerObject> passingMuons;
};
#endif
