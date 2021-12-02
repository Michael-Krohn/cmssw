#ifndef MCHistograms_h
#define MCHistograms_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DarkPhoton/MuAnalyzer/interface/SigEventInfo.h"

class TH1F;
class TH2F;

class MCHistograms{
  public:
    MCHistograms();

    void book(TFileDirectory histFolder);
    void FillHists(SigEventInfo info);
    void Normalize();
    void IncCutFlow(double weight);
    void FailLocation(double ieta, double iphi, double weight);
    void ResetCutFlow();
  public:
    TH1F* m_eventCount;
    TH1F* m_cutProgress;
    TH1F* m_signalCategory;
    TH1F* m_ProbeTrackPt;
    TH1F* m_ProbeTrackEta;
    TH1F* m_ProbeMuonEta;
    TH1F* m_ProbeMuonPt;
    TH1F* m_ConeHits;
    TH1F* m_ConeEnergy;
    TH1F* m_nMatchedMuons;
    TH1F* m_HitsOverThresh;
    TH2F* m_DBremLocation;
    TH1F* m_DeflectingAngle;
    TH1F* m_initialMuE;
    TH1F* m_finalMuE;
    TH1F* m_DBremScalingDistance;
    TH1F* m_ProbeTrackIsolation;
    TH1F* m_ProbeEcalIsolation;
    TH2F* m_IsolationVariables;
    TH1F* m_TrackMotherdE;
    TH1F* m_TrackMotherdR;
    TH1F* m_FractionalELost;
    int cutProgress;
    TH1F* m_BremDepth;
    TH1F* m_CSC_dR;
    TH1F* m_MuonMotherdE;
    TH1F* m_MuonMotherdR;
    TH1F* m_CSCHitChiSq;
    TH1F* m_CSCHitAngleChange;
    TH1F* m_StandaloneMuonE;
    TH1F* m_GlobalMuonE;
    TH1F* m_StandalonedE;
    TH2F* m_MatchedLocation;
    TH1F* m_NStandalone;
    TH2F* m_NMatchNHitLocation;
    TH1F* m_NMatchNHitHCALEnergy;
    TH1F* m_AllVertexOffset;
    TH1F* m_EventWeights;
    TH1F* m_SignalSelectionCuts;
    TH1F* m_DBremR;
    TH1F* m_HEHitByDepth;
    TH1F* m_DrtoStandalone;
    TH1F* m_VtxZ;
    TH1F* m_ProbePFJetE;
    TH1F* m_ProbePFJetDr;
    TH1F* m_ProbePFJetMaxE;
    TH1F* m_CaloJetMinDr;
    TH1F* m_CaloJetEcalE;
    TH1F* m_CaloJetHcalE;
    TH1F* m_CaloJetTotalE;
    TH1F* m_CaloConeEnergy;
    TH1F* m_CaloConeHits;
    TH1F* m_HECellsFound;
    TH2F* m_HEEnergyStandaloneDEoverE;
    TH2F* m_FailLocation;
/*    TH2F* m_BigEtaAdjFaildEtadPhiMuPlusEtaPlus;
    TH2F* m_BigEtaAdjFaildEtadPhiMuPlusEtaMinus;
    TH2F* m_BigEtaAdjFaildEtadPhiMuMinusEtaPlus;
    TH2F* m_BigEtaAdjFaildEtadPhiMuMinusEtaMinus;
    TH2F* m_BigEtaAlldEtadPhiMuPlusEtaPlus;
    TH2F* m_BigEtaAlldEtadPhiMuPlusEtaMinus;
    TH2F* m_BigEtaAlldEtadPhiMuMinusEtaPlus;
    TH2F* m_BigEtaAlldEtadPhiMuMinusEtaMinus;
    TH2F* m_SmallEtaAdjFaildEtadPhiMuPlusEtaPlus;
    TH2F* m_SmallEtaAdjFaildEtadPhiMuPlusEtaMinus;
    TH2F* m_SmallEtaAdjFaildEtadPhiMuMinusEtaPlus;
    TH2F* m_SmallEtaAdjFaildEtadPhiMuMinusEtaMinus;
    TH2F* m_SmallEtaAlldEtadPhiMuPlusEtaPlus;
    TH2F* m_SmallEtaAlldEtadPhiMuPlusEtaMinus;
    TH2F* m_SmallEtaAlldEtadPhiMuMinusEtaPlus;
    TH2F* m_SmallEtaAlldEtadPhiMuMinusEtaMinus;*/
    TH1F* m_MinCSCImpactParameter;
    TH1F* m_DiMuonMass;
    TH1F* m_TagMuonEta;
    TH1F* m_TagMuonPt;
    TH1F* m_nSelectedMuons;
    TH1F* m_PairVtxChi;
    TH1F* m_TrackIsolation;
    TH1F* m_TrackQuality;
    TH1F* m_TrackReducedChi;
    TH1F* m_HcalPlusEcal;
    TH1F* m_Probedxy;
    TH1F* m_Probedsz;
};
#endif
