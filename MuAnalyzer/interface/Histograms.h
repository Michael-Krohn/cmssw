#ifndef Histograms_h
#define Histograms_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


class TH1F;
class TH2F;

class Histograms{
  public:
    Histograms();

    void book(edm::Service<TFileService> fs);
    void PlotTrackDisappearance(double TrackP, double TrackEta, double TrackPhi, double minDR, double minTotalImpactParameter, double TrackP_dR, double TrackEta_dR, double TrackPhi_dR);

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
    TH1F* m_histogram_differenceVertex_z;

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

    TH2F* m_histogram_TrackerTrackMatched_EtaPhi_posEta_dR0p35;
    TH2F* m_histogram_TrackerTrackMatched_EtaPhi_negEta_dR0p35;

    TH1F* m_MinDR_MuonHCAL;
    TH1F* m_HitEnergy_MinDR_MuonHCAL;

};

#endif
