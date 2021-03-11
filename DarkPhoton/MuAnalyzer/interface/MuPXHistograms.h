#ifndef MuPXHistograms_h
#define MuPXHistograms_h

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
#include "DarkPhoton/MuAnalyzer/interface/MuPXEventInfo.h"

class TH1F;
class TH2F;

class MuPXHistograms{
  public:
    MuPXHistograms();

    void book(TFileDirectory histFolder);
    void IncCutFlow();
    void ResetCutFlow();
    void FillHists(MuPXEventInfo);
  public:
    TH1F* m_eventWeight;
    TH1F* m_eventCount;
    TH1F* m_cutProgress;
    TH1F* m_MuonTrackMass;
    TH1F* m_ProbeEta;
    TH1F* m_ProbePt;
    TH1F* m_ProbePhi;
    TH2F* m_ProbeEtaPhi;
    TH1F* m_NPassingProbe;
    TH1F* m_TagEta;
    TH1F* m_TagPt;
    TH1F* m_TagPhi;
    TH2F* m_TagEtaPhi;
    TH1F* m_NPassingTag;
    TH1F* m_TagProbeVtxChi;
    TH1F* m_ProbeTrackIso;
    TH1F* m_ProbeHcalIso;
    TH1F* m_ProbeEcalIso;
    TH2F* m_ProbeCombinedIso;     
    double cutProgress;
    //Multiple Paired Track Hists
    TH1F* m_SmallestCone;
    TH1F* m_AverageDr;
    //Study Tag Muon selection
    TH1F* m_TagProbeDr;
    TH1F* m_TagTrackIso;
    TH1F* m_TagHcalIso;
    TH1F* m_TagEcalIso;
    TH1F* m_NJets;
    TH2F* m_CaloSumJetE;
    TH2F* m_HCALIsoJetDr;
    TH2F* m_ECALIsoJetDr;
    TH1F* m_ProbeJetDr;
    TH1F* m_JetPt;
    //Calo Jet 
    TH1F* m_ProbeJetE;
    TH1F* m_TagCaloJetE;
    TH2F* m_CaloJetEcalE;
    TH2F* m_CaloJetHcalE;
    //Vtx Stuff
    TH1F* m_ProbeVtx;
    TH1F* m_TagVtx;
    TH1F* m_NVertices;
    TH2F* m_TagProbeVtx;
    TH1F* m_TagProbeVtxd0;
    TH1F* m_TagProbeVtxdz;
    TH1F* m_PileupWeights;
    TH1F* m_PUmean;
};
#endif
