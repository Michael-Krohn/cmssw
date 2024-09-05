#ifndef eventHistos_h
#define eventHistos_h

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

class eventHistos
{
   public:
      eventHistos();
      void book(TFileDirectory histFolder);
      void IncCutFlow();
      void ResetCutFlow();
      void FillTree(int runNum, int lumiSection, int eventNum, int sigCat);
      TH1F* m_eventWeight;
      TH1F* m_eventCount;
      TH1F* m_staMinDr;
      TH1F* m_staMinDEoverE;
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
      TH1F* m_NJets;
      TH1F* m_ProbeEcalIso;
      TH2F* m_ProbeCombinedIso;     
      TTree* m_passingEventInfo;
      double cutProgress;
      int run;
      int lumi;
      int event;
      int category;
};
#endif
