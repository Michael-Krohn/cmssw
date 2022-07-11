#ifndef Histograms_h
#define Histograms_h

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
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"

class TH1F;
class TH2F;

class Histograms {
public:
  Histograms();

  void book(TFileDirectory histFolder);
  void PlotTrackDisappearance(double TrackP,
                              double TrackEta,
                              double TrackPhi,
                              double minDR,
                              double minTotalImpactParameter,
                              double TrackP_dR,
                              double TrackEta_dR,
                              double TrackPhi_dR);
  void PlotCSCHits(const edm::Event& iEvent,
                   const edm::EventSetup& iSetup,
                   edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label);
  void PlotHCALHits(
      const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label);
  void Normalize();
  void IncCutFlow(double weight);
  void IncCutFlow(double weight, int nCuts);
  void FillHists(EventInfo);

public:
  int cutProgress;
  TH1F* m_eventCount;
  TH1F* m_cutProgress;
  TH1F* m_MuonTrackMass;
  TH1F* m_failHEEnergy;
  TH1F* m_histogram_TrackerTack_P;
  TH1F* m_probePt;
  TH1F* m_CSCMinDR;
  TH1F* m_CSCMinDepthDr[4];
  TH1F* m_CSCOtherStationsGood[4];
  TH2F* m_CSCMissingDepthDr[4];
  TH2F* m_CSCGoodDepthDr[4];
  TH2F* m_CSCOtherGoodMissingDepth[4];
  //Tracker track location
  TH2F* m_histogram_TrackerTrack_EtaPhi_negEta;
  TH2F* m_histogram_TrackerTrack_EtaPhi_negEta_pT0to50;
  TH2F* m_histogram_TrackerTrack_EtaPhi_negEta_pT50to100;
  TH2F* m_histogram_TrackerTrack_EtaPhi_negEta_pT100to150;
  TH2F* m_histogram_TrackerTrack_EtaPhi_negEta_pT150toInf;
  TH2F* m_histogram_TrackerTrack_EtaPhi_posEta;
  TH2F* m_histogram_TrackerTrack_EtaPhi_posEta_pT0to50;
  TH2F* m_histogram_TrackerTrack_EtaPhi_posEta_pT50to100;
  TH2F* m_histogram_TrackerTrack_EtaPhi_posEta_pT100to150;
  TH2F* m_histogram_TrackerTrack_EtaPhi_posEta_pT150toInf;

  //Tracker track matched to muon CSC segment
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf;

  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP10;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP10;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP5;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP5;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP3;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP3;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP2;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP2;

  TH2F* m_histogram_TrackerTrack_dRMatched_EtaPhi_negEta;
  TH2F* m_histogram_TrackerTrack_dRMatched_EtaPhi_posEta;

  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_dR0p1;
  TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_dR0p1;

  //CSC Segment Location
  TH2F* m_histogram_CSCHits_EtaPhi;

  //HCAL hits histogram
  TH2F* m_HCALHits_EtaPhi;
  TH1F* m_HEConeEnergy;
  TH1F* m_tagHcalIso;
  TH1F* m_probeHcalIso;
  TH1F* m_centerInCone;
  TH1F* m_centerOutOfCone;
  TH1F* m_MinDR_MuonHCAL;
  TH1F* m_MinDR_RandomHCAL;
  TH1F* m_HitDepth_MuonHCAL;
  TH1F* m_HitDepth_RandomHCAL;
  TH1F* m_ConeHits;
  TH1F* m_ConeEnergy;
  TH1F* m_RandomConeHits;
  TH1F* m_RandomConeEnergy;
  TH1F* m_HitsOverThresh;
  TH1F* m_RandomHitsOverThresh;
  TH1F* m_Layer_Spectra[7];
  TH1F* m_RLayer_Spectra[7];
  TH1F* m_MissingHits;
  TH2F* m_DepthPairSpectra[6];
  TH2F* m_RDepthPairSpectra[6];
  TH2F* m_HcalMissing[7];
  TH2F* m_HcalGood[7];
  TH2F* m_Layer_Eta[7];
  TH2F* m_RLayer_Eta[7];
  TH1F* m_ValidIDs;
  TH1F* m_HEHitsDR;
  TH2F* m_HECellEdgeDistance[7];
  TH2F* m_HESmallCellEdgeDistance[7];
  TH2F* m_HEBigCellDEtaDPhi[7];
  TH2F* m_HESmallCellDEtaDPhi[7];
  TH2F* m_HEBigCellDEtaDPhiPosQ[7];
  TH2F* m_HESmallCellDEtaDPhiPosQ[7];
  TH2F* m_HEBigCellDEtaDPhiNegQ[7];
  TH2F* m_HESmallCellDEtaDPhiNegQ[7];
  TH2F* m_etaPhiEdgeMin;
  TH2F* m_smallEtaPhiEdgeMin;
};

#endif
