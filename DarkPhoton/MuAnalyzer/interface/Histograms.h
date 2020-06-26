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

class TH1F;
class TH2F;

class Histograms{
  public:
    Histograms();

    void book(edm::Service<TFileService> fs);
    void PlotTrackDisappearance(double TrackP, double TrackEta, double TrackPhi, double minDR, double minTotalImpactParameter, double TrackP_dR, double TrackEta_dR, double TrackPhi_dR);
    void PlotCSCHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label);
    void PlotHCALHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label);
    void Normalize();

  public:
    TH1F* m_eventCount;
    TH1F* m_cutProgress;
    
    TH1F* m_MuonTrackMass;
    
    TH1F* m_GENMuonTrackdR;
    TH1F* m_nMuonTrackCand;
    TH1F* m_nTracksPairedPerMuon;
    TH1F* m_nMuonsPairedPerEvent;
    TH1F* m_nTracksNoMuon;

    TH1F* m_histogram_TrackerTack_P;
    TH1F* m_histogram_MuonTrack_P;

    TH2F* m_histogram_differenceVertex_xy;
    TH2F* m_histogram_misshitVertex_xy;
    TH1F* m_histogram_differenceVertex_z;
    TH1F* m_histogram_misshitVertex_z;
    TH1F* m_adjacentfailvertex_z;

    TH1F* m_MinTransverseImpactParameter;
    TH1F* m_MinLongiudinalImpactParameter;
    TH1F* m_MinTotalImpactParameter;
    TH1F* m_MinDR;

    TH1F* m_MinDR_Muon;
    TH1F* m_MinTotalImpactParameterMuon;

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
    TH2F* m_histogram_HCALHits_EtaPhi;
    TH2F* m_histogram_BlankHCALHits_EtaPhi;

    TH1F* m_MinDR_MuonHCAL;
    TH1F* m_MinDR_RandomHCAL;
    TH1F* m_HitEnergy_MinDR_MuonHCAL;
    TH1F* m_HitEnergy_RandomHCAL;
    TH1F* m_HitDepth_MuonHCAL;
    TH1F* m_HitDepth_RandomHCAL;
    TH1F* m_ConeHits;
    TH1F* m_ConeEnergy;
    TH1F* m_MissThreshConeEnergy;
    TH1F* m_RandomConeHits;
    TH1F* m_RandomConeEnergy;   
    TH1F* m_HitsOverThresh;
    TH1F* m_RandomHitsOverThresh;
    TH1F* m_Layer_Spectra[7];
    TH1F* m_RLayer_Spectra[7];
    TH2F* m_DepthPairSpectra[6];
    TH2F* m_RDepthPairSpectra[6];
    TH2F* m_Layer_Eta[7];
    TH2F* m_RLayer_Eta[7];
    TH1F* m_ValidIDs;
    TH1F* m_Missing_ValidIDs;
    TH1F* m_MissingHits;
    TH1F* m_RMissingHits;
    TH2F* m_MissingHitsMap;
    TH2F* m_RMissingHitsMap;
    TH1F* m_MissingHitsEnergy;
    TH1F* m_MissingHitsDR;
    TH1F* m_BlankHitsDR;
    TH1F* m_TrackHCALDR_MissHit;
    TH1F* m_TrackHCALDR_GoodHits;
    TH1F* m_BlankDepth;
    TH1F* m_4BlankDepth;
    TH2F* m_BlankCellDetaDphiPosEta;
    TH2F* m_BlankCellSmallDetaDphiPosEta;
    TH2F* m_BlankCellDetaDphiNegEta;
    TH2F* m_BlankCellSmallDetaDphiNegEta;   
    TH1F* m_MuonEtaDist;
    TH1F* m_TrackPt;
    TH1F* m_MissHitTrackPt;
    TH2F* m_DBremLocation;
    TH1F* m_DeflectingAngle;
    TH1F* m_initialMuE;
    TH1F* m_finalMuE;
    TH1F* m_dphoEnergy;
    TH1F* m_NThreshCut;
};

#endif
