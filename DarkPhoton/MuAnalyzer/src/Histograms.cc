#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "TH2F.h"
#include "TTree.h"
#include <iostream>

Histograms::Histograms() {}

void Histograms::book(TFileDirectory histFolder) {
  m_eventCount = histFolder.make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
  m_cutProgress = histFolder.make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 10, -.5, 9.5);
  m_MuonTrackMass = histFolder.make<TH1F>("MuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50, 150);
  m_treeVariables = histFolder.make<TTree>("treeVariables", "Tree containing cell edge variables");
  m_treeVariables->Branch("charge",&charge);
  m_treeVariables->Branch("hitsOverThresh", &hitsOverThresh);
  m_treeVariables->Branch("vertexZ", &vertexZ);
  m_treeVariables->Branch("depthsFound", &depthsFound);
  m_treeVariables->Branch("iEta", &iEta);
  m_treeVariables->Branch("eta", &eta);
  m_treeVariables->Branch("phi", &phi);
  m_treeVariables->Branch("pt", &pt);
  for(int i=0;i<7;i++){
    std::string etaName = "cellDetaDepth_" + std::to_string(i);
    std::string phiName = "cellDphiDepth_" + std::to_string(i);
    m_treeVariables->Branch(etaName.c_str(), &cellDetaDepth[i]);
    m_treeVariables->Branch(phiName.c_str(), &cellDphiDepth[i]);
    std::string hitName = "hitEnergiesDepth" + std::to_string(i);
    m_treeVariables->Branch(hitName.c_str(), &hitEnergies[i]);
    m_treeVariables->Branch(("neighborEta"+hitName).c_str(), &neighborHitEnergies[i]);
    m_treeVariables->Branch(("neighborPHi"+hitName).c_str(), &neighborHitEnergiesPhi[i]);
    std::string shiftName = "depth_"+std::to_string(i);
    m_treeVariables->Branch(("xShiftEnergy"+shiftName).c_str(), &xShiftEnergies[i]);
    m_treeVariables->Branch(("yShiftEnergy"+shiftName).c_str(), &yShiftEnergies[i]);
    m_treeVariables->Branch(("zShiftEnergy"+shiftName).c_str(), &zShiftEnergies[i]);
    m_treeVariables->Branch(("xShifts"+shiftName).c_str(), &xShifts[i]);
    m_treeVariables->Branch(("yShifts"+shiftName).c_str(), &yShifts[i]);
    m_treeVariables->Branch(("zShifts"+shiftName).c_str(), &zShifts[i]);
  }
  m_failHEEnergy = histFolder.make<TH1F>("failHEEnergy", ";HE Energy (GeV);Events", 200, 0, 200);
  m_histogram_TrackerTack_P = histFolder.make<TH1F>("TrackerTrack_P", "", 100, 0, 1000);
  m_probePt = histFolder.make<TH1F>("probePt", ";Probe Track Pt (GeV);Events", 100, 0, 100);
  m_CSCMinDR = histFolder.make<TH1F>("CSCMinDR", "", 140, 0, 3);
  for (int i = 0; i < 4; i++) {
    std::string histName = "CSCMinDR_Depth" + std::to_string(i);
    m_CSCMinDepthDr[i] = histFolder.make<TH1F>(histName.c_str(), ";#DeltaR to Nearest CSC Hit;Total Weight", 140, 0, 3);
    m_CSCRecMinDepthDr[i] = histFolder.make<TH1F>((histName+"rec").c_str(), ";#DeltaR to Nearest CSC Hit;Total Weight", 140, 0, 3);
    m_RpcMinDepthDr[i] = histFolder.make<TH1F>((histName+"rpc").c_str(), ";#DeltaR to Nearest CSC Hit;Total Weight", 140, 0, 3);

    histName = histName + "missing";
    m_CSCMissingDepthDr[i] = histFolder.make<TH2F>(
        histName.c_str(), ";Tag Muon #eta; Tag Muon #phi; Total Weight", 100, -2.6, 2.6, 100, -3.15, 3.15);
    m_CSCGoodDepthDr[i] = histFolder.make<TH2F>(
        (histName + "good").c_str(), ";Tag Muon #eta; Tag Muon #phi; Total Weight", 100, -2.6, 2.6, 100, -3.15, 3.15);
    m_CSCRecMissingDepthDr[i] = histFolder.make<TH2F>(
        (histName+"rec").c_str(), ";Tag Muon #eta; Tag Muon #phi; Total Weight", 100, -2.6, 2.6, 100, -3.15, 3.15);
    m_CSCRecGoodDepthDr[i] = histFolder.make<TH2F>(
        (histName + "recgood").c_str(), ";Tag Muon #eta; Tag Muon #phi; Total Weight", 100, -2.6, 2.6, 100, -3.15, 3.15);
    m_RpcMissingDepthDr[i] = histFolder.make<TH2F>(
        (histName+"rpc").c_str(), ";Tag Muon #eta; Tag Muon #phi; Total Weight", 100, -2.6, 2.6, 100, -3.15, 3.15);
    m_RpcGoodDepthDr[i] = histFolder.make<TH2F>(
        (histName + "rpcgood").c_str(), ";Tag Muon #eta; Tag Muon #phi; Total Weight", 100, -2.6, 2.6, 100, -3.15, 3.15);

    histName = "CSCMinDr_Depth" + std::to_string(i) + "_otherDepthsGood";
    m_CSCOtherStationsGood[i] =
        histFolder.make<TH1F>(histName.c_str(), ";#DeltaR to Nearest CSC Hit;Total Weight", 140, 0, 3);
    histName = histName + "missing";
    m_CSCOtherGoodMissingDepth[i] = histFolder.make<TH2F>(
        histName.c_str(), ";Tag Muon #eta; Tag Muon #phi; Total Weight", 100, -2.6, 2.6, 100, -3.15, 3.15);
  }
  m_CSCMinDR = histFolder.make<TH1F>("CSCMinDR", "", 140, 0, 3);
  m_CSCMinDR = histFolder.make<TH1F>("CSCMinDR", "", 140, 0, 3);

  m_histogram_TrackerTrack_EtaPhi_negEta =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT0to50 =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT0to50", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT50to100 =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT50to100", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT100to150 =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT100to150", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT150toInf =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT150toInf", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT0to50 =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT0to50", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT50to100 =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT50to100", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT100to150 =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT100to150", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT150toInf =
      histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT150toInf", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_CSCHits_EtaPhi = histFolder.make<TH2F>(
      "CSCHits_EtaPhi", "; #eta; #phi; Events", 100, -2.4, 2.4, 72, -ROOT::Math::Pi(), ROOT::Math::Pi());
  m_HCALHits_EtaPhi =
      histFolder.make<TH2F>("HCALHits_EtaPhi", "; i#eta; i#phi; Events", 57, -28.5, 28.5, 72, 0.5, 72.5);

  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_IP15_posEta", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP10 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP10", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP10 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP10", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP5 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP5", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP5 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP5", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP3 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP3", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP3 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP3", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP2 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP2", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP2 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP2", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_TrackerTrack_dRMatched_EtaPhi_negEta =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_dR0p1", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_dRMatched_EtaPhi_posEta =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_dR0p1", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_TrackerTrackMatched_EtaPhi_posEta_dR0p1 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_dR0p1", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_dR0p1 =
      histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_dR0p1", "", 100, -2.4, -1.63, 100, -3.2, 3.2);

  m_centerInCone = histFolder.make<TH1F>("CenterInCone", ";Projected cells in 0.2 DR cone", 10, -0.5, 9.5);
  m_centerOutOfCone = histFolder.make<TH1F>("CenterOutOfCone", ";Projected cells out of 0.2 DR cone", 10, -0.5, 9.5);
  m_MinDR_MuonHCAL = histFolder.make<TH1F>("MinDR_MuonHCAL", "; Min #Delta R; Events", 140, 0, 0.4);
  m_MinDR_RandomHCAL = histFolder.make<TH1F>("MinDR_RandomHCAL", "; Min #Delta R; Events", 140, 0, 0.4);
  m_HitDepth_MuonHCAL = histFolder.make<TH1F>("HitDepth_MuonHCAL", "; Depth; Events", 40, 0, 20);
  m_HitDepth_RandomHCAL = histFolder.make<TH1F>("HitDepth_RandomHCAL", "; Depth; Events", 40, 0, 20);
  m_ConeHits = histFolder.make<TH1F>("ConeHits", "; Hits in Cone; Events", 50, -0.5, 49.5);
  m_ConeEnergy = histFolder.make<TH1F>("ConeEnergy", "; Energy in Matched Hits (GeV); Events", 200, 0, 20);
  m_RandomConeHits = histFolder.make<TH1F>("RandomConeHits", "; Hits in Cone; Events", 50, -0.5, 49.5);
  m_RandomConeEnergy = histFolder.make<TH1F>("RandomConeEnergy", "; Energy in Cone (GeV); Events", 200, 0, 20);
  m_HitsOverThresh = histFolder.make<TH1F>("HitsOverThresh", "; Hits Over Threshold; Events", 8, -0.5, 7.5);
  m_MissingHits = histFolder.make<TH1F>("MissingHits", "; Missing Hits; Events", 8, -0.5, 7.5);
  m_HEConeEnergy =
      histFolder.make<TH1F>("HEConeEnergy", ";HCAL Energy in 0.4 #Delta R Cone (GeV); Events", 100, 0, 100);
  m_tagHcalIso = histFolder.make<TH1F>("tagHcalIso", ";HCAL Energy in 0.4 #Delta R Cone (GeV); Events", 100, 0, 100);
  m_probeHcalIso =
      histFolder.make<TH1F>("probeHcalIso", ";HCAL Energy in 0.4 #Delta R Cone (GeV); Events", 100, 0, 100);
  m_RandomHitsOverThresh = histFolder.make<TH1F>("RandomHitsOverThresh", "; Hits Over Threshold; Events", 8, -0.5, 7.5);

  for (int i = 0; i < 7; i++) {
    std::string name = "Layer" + std::to_string(i + 1) + "Spectra";
    std::string plusname = name + "PlusNeighbor";
    std::string rname = "R" + name;
    std::string pairname = "DepthPairSpectra_" + std::to_string(i + 1) + std::to_string(i + 2);
    std::string rpairname = "R" + pairname;
    std::string pairaxis = "; Depth " + std::to_string(i + 1) + " Hit Energy (GeV); Depth " + std::to_string(i + 2) +
                           " Hit Energy (GeV); Events";
    m_Layer_Spectra[i] = histFolder.make<TH1F>(name.c_str(), "; Hit Energy (GeV); Events", 400, 0, 10);
    m_Layer_PlusNeighbor_Spectra[i] = histFolder.make<TH1F>(plusname.c_str(), "; Energy (GeV); Events", 200, 0, 10);
    m_RLayer_Spectra[i] = histFolder.make<TH1F>(rname.c_str(), "; Hit Energy (GeV); Events", 400, 0, 10);
    if (i < 6) {
      m_DepthPairSpectra[i] = histFolder.make<TH2F>(pairname.c_str(), pairaxis.c_str(), 100, 0, 10, 100, 0, 10);
      m_RDepthPairSpectra[i] = histFolder.make<TH2F>(rpairname.c_str(), pairaxis.c_str(), 100, 0, 10, 100, 0, 10);
    }
    std::string twoDname = "Depth" + std::to_string(i + 1) + "HcalMissing";
    m_HcalMissing[i] = histFolder.make<TH2F>(
        twoDname.c_str(), ";Tag Muon #eta;Tag Muon #phi; Events", 100, -2.6, 2.6, 100, -3.15, 3.15);
    std::string goodName = "Depth" + std::to_string(i + 1) + "HcalGood";
    m_HcalGood[i] = histFolder.make<TH2F>(
        goodName.c_str(), ";Tag Muon #eta;Tag Muon #phi; Events", 100, -2.6, 2.6, 100, -3.15, 3.15);
  }
  for (int i = 0; i < 7; i++) {
    std::string name = "Layer" + std::to_string(i + 1) + "Eta";
    std::string rname = "R" + name;
    m_Layer_Eta[i] =
        histFolder.make<TH2F>(name.c_str(), "; Hit Energy (GeV); #eta ; Events", 150, 0, 10, 50, -2.6, 2.6);
    m_RLayer_Eta[i] =
        histFolder.make<TH2F>(rname.c_str(), "; Hit Energy (GeV); #eta ; Events", 150, 0, 10, 50, -2.6, 2.6);
  }
  m_ValidIDs = histFolder.make<TH1F>("ValidCellIds", "; Valid Cell Ids; Events", 101, -0.5, 100.5);
  for (int depth = 0; depth < 7; depth++) {
    std::string histName = "CellCenterDistanceDepth" + std::to_string(depth);
    std::string edgeName = "CellEdgeDistanceDepth" + std::to_string(depth);
    m_HEBigCellDEtaDPhi[depth] = histFolder.make<TH2F>(
        ("BigEta" + histName).c_str(), "; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
    m_HECellEdgeDistance[depth] = histFolder.make<TH2F>(
        ("BigEta" + edgeName).c_str(), "; #Delta#eta; #Delta#phi; Events", 50, -0.1, 0.1, 50, -0.1, 0.1);
    m_HESmallCellEdgeDistance[depth] = histFolder.make<TH2F>(
        ("SmallEta" + edgeName).c_str(), "; #Delta#eta; #Delta#phi; Events", 70, -0.07, 0.07, 70, -0.07, 0.07);
    m_HESmallCellDEtaDPhi[depth] = histFolder.make<TH2F>(
        ("SmallEta" + histName).c_str(), "; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
    m_HEBigCellDEtaDPhiPosQ[depth] = histFolder.make<TH2F>(
        ("PosQBigEta" + histName).c_str(), "; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
    m_HESmallCellDEtaDPhiPosQ[depth] = histFolder.make<TH2F>(
        ("PosQSmallEta" + histName).c_str(), "; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
    m_HEBigCellDEtaDPhiNegQ[depth] = histFolder.make<TH2F>(
        ("NegQBigEta" + histName).c_str(), "; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
    m_HESmallCellDEtaDPhiNegQ[depth] = histFolder.make<TH2F>(
        ("NegQSmallEta" + histName).c_str(), "; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  }
  m_etaPhiEdgeMin = histFolder.make<TH2F>(
      "EtaPhiEdgeMin", "; Min #Delta#eta to Cell Edge; Min #Delta#phi to Cell Edge; Events", 50, 0, 0.1, 50, 0, 0.1);
  m_smallEtaPhiEdgeMin = histFolder.make<TH2F>("SmallCellEtaPhiEdgeMin",
                                               "; Min #Delta#eta to Cell Edge; Min #Delta#phi to Cell Edge; Events",
                                               50,
                                               0,
                                               0.1,
                                               50,
                                               0,
                                               0.1);
}

void Histograms::FillTree(EventInfo info) {
  pt = info.probePt;
  hitsOverThresh = info.hitsOver0p2;
  depthsFound = info.expectedHits;
  iEta = info.probeiEta;
  eta = info.ProbeEta;
  phi = info.ProbePhi;
  charge = info.probeCharge;
  vertexZ = info.vertexZ;
  for(int i=0;i<7;i++){
    if(info.foundDepths[i])
    {
      cellDetaDepth[i] = info.cellEtaEdgeDistance[i];
      cellDphiDepth[i] = info.cellPhiEdgeDistance[i];
      hitEnergies[i] = info.hitEnergies[i];
      neighborHitEnergies[i] = info.neighborHitEnergies[i];
      neighborHitEnergiesPhi[i] = info.neighborHitEnergiesPhi[i];
      xShifts[i] = info.xShifts[i];
      yShifts[i] = info.yShifts[i];
      zShifts[i] = info.zShifts[i];
      xShiftEnergies[i] = info.xShiftEnergies[i];
      yShiftEnergies[i] = info.yShiftEnergies[i];
      zShiftEnergies[i] = info.zShiftEnergies[i];
    }
    else
    {
      cellDetaDepth[i] = -1;
      cellDphiDepth[i] = -1;
      hitEnergies[i] = -1;
      neighborHitEnergies[i] = -1;
      neighborHitEnergiesPhi[i] = -1;
      xShifts[i] = 10;
      yShifts[i] = 10;
      xShiftEnergies[i] = -1;
      yShiftEnergies[i] = -1;
    }
  }
  m_treeVariables->Fill();
}

void Histograms::FillHists(EventInfo info) {
  m_eventCount->Fill(0.5, info.eventWeight);
  m_MuonTrackMass->Fill(info.diMuonMass, info.eventWeight);
  m_probePt->Fill(info.probePt, info.eventWeight);
  m_CSCMinDR->Fill(info.minDR, info.eventWeight);
  bool depthGood[4];
  for (int i = 0; i < 4; i++) {
    m_CSCMinDepthDr[i]->Fill(info.minDrByDepth[i], info.eventWeight);
    m_CSCRecMinDepthDr[i]->Fill(info.minRecDrByDepth[i], info.eventWeight);
    m_RpcMinDepthDr[i]->Fill(info.minRpcDrByDepth[i], info.eventWeight);
    if (fabs(info.minDrByDepth[i]) > 0.05) {
      m_CSCMissingDepthDr[i]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
      depthGood[i] = false;
    } else {
      m_CSCGoodDepthDr[i]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
      depthGood[i] = true;
    }
    if (fabs(info.minRecDrByDepth[i]) > 0.05) {
      m_CSCRecMissingDepthDr[i]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
    } else {
      m_CSCRecGoodDepthDr[i]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
    }
    if (fabs(info.minRpcDrByDepth[i]) > 0.05) {
      m_RpcMissingDepthDr[i]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
    } else {
      m_RpcGoodDepthDr[i]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
    }

  }
  for (int depth = 0; depth < 7; depth++) {
    if (info.hitEnergies[depth] < 0.2) {
      m_HcalMissing[depth]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
    } else {
      m_HcalGood[depth]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
    }
  }
  for (int i = 0; i < 4; i++) {
    bool goodForChecking = true;
    for (int j = 0; j < 4; j++) {
      if (j == i)
        continue;
      if (!depthGood[j])
        goodForChecking = false;
    }
    if (goodForChecking) {
      m_CSCOtherStationsGood[i]->Fill(info.minDrByDepth[i], info.eventWeight);
      if (info.minDrByDepth[i] > 0.05) {
        m_CSCOtherGoodMissingDepth[i]->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
      }
    }
  }
  m_HCALHits_EtaPhi->Fill(info.ProbeEta, info.ProbePhi, info.eventWeight);
  m_MinDR_MuonHCAL->Fill(info.minDR_muonHCAL, info.eventWeight);
  m_MinDR_RandomHCAL->Fill(info.minDR_randomHCAL, info.eventWeight);
  m_HitDepth_MuonHCAL->Fill(info.hitDepth_muonHCAL, info.eventWeight);
  m_HitDepth_RandomHCAL->Fill(info.hitDepth_randomHCAL, info.eventWeight);
  m_centerInCone->Fill(info.centerInCone, info.eventWeight);
  m_centerOutOfCone->Fill(info.centerOutOfCone, info.eventWeight);
  m_ConeHits->Fill(info.coneHits, info.eventWeight);
  m_ConeEnergy->Fill(info.coneEnergy, info.eventWeight);
  m_HEConeEnergy->Fill(info.heConeEnergy);
  m_tagHcalIso->Fill(info.tagHcalIso);
  m_probeHcalIso->Fill(info.probeHcalIso);
  m_RandomConeHits->Fill(info.randomConeHits, info.eventWeight);
  m_RandomConeEnergy->Fill(info.randomConeEnergy, info.eventWeight);
  m_HitsOverThresh->Fill(info.hitsOverThresh, info.eventWeight);
  m_RandomHitsOverThresh->Fill(info.randomHitsOverThresh, info.eventWeight);
  for (int i = 0; i < 7; i++) {
    if (!info.foundDepths[i]) {
      continue;
    }
    m_Layer_Spectra[i]->Fill(info.hitEnergies[i], info.eventWeight);
    m_Layer_PlusNeighbor_Spectra[i]->Fill(
        info.hitEnergies[i] + info.neighborHitEnergies[i] + info.neighborHitEnergiesPhi[i], info.eventWeight);
    m_RLayer_Spectra[i]->Fill(info.randHitEnergies[i], info.eventWeight);
    m_Layer_Eta[i]->Fill(info.hitEnergies[i], info.ProbeEta, info.eventWeight);
    //      m_RLayer_Eta[i]->Fill(info.randHitEnergies[i],info.ProbeEta,info.eventWeight);
    if (i < 6 && ((i + 1) < info.cellsFound)) {
      m_DepthPairSpectra[i]->Fill(info.hitEnergies[i], info.hitEnergies[i + 1], info.eventWeight);
      m_RDepthPairSpectra[i]->Fill(info.randHitEnergies[i], info.randHitEnergies[i + 1], info.eventWeight);
    }
  }
  m_MissingHits->Fill(info.expectedHits - info.hitsOverThresh, info.eventWeight);
  m_ValidIDs->Fill(info.cellsFound, info.eventWeight);
  if (fabs(info.probeiEta) > 20) {
    m_etaPhiEdgeMin->Fill(info.etaEdgeMin, info.phiEdgeMin, info.eventWeight);
  } else {
    m_smallEtaPhiEdgeMin->Fill(info.etaEdgeMin, info.phiEdgeMin, info.eventWeight);
  }
  for (int depth = 0; depth < 7; depth++) {
    if (!info.foundDepths[depth]) {
      continue;
    }
    if (fabs(info.probeiEta) > 20) {
      m_HEBigCellDEtaDPhi[depth]->Fill(info.cellCenterDEta[depth], info.cellCenterDPhi[depth], info.eventWeight);
      m_HECellEdgeDistance[depth]->Fill(
          info.cellEtaEdgeDistance[depth], info.cellPhiEdgeDistance[depth], info.eventWeight);
      if (info.probeCharge > 0) {
        m_HEBigCellDEtaDPhiPosQ[depth]->Fill(info.cellCenterDEta[depth], info.cellCenterDPhi[depth], info.eventWeight);
      } else {
        m_HEBigCellDEtaDPhiNegQ[depth]->Fill(info.cellCenterDEta[depth], info.cellCenterDPhi[depth], info.eventWeight);
      }
    } else {
      m_HESmallCellDEtaDPhi[depth]->Fill(info.cellCenterDEta[depth], info.cellCenterDPhi[depth], info.eventWeight);
      m_HESmallCellEdgeDistance[depth]->Fill(
          info.cellEtaEdgeDistance[depth], info.cellPhiEdgeDistance[depth], info.eventWeight);

      if (info.probeCharge > 0) {
        m_HESmallCellDEtaDPhiPosQ[depth]->Fill(
            info.cellCenterDEta[depth], info.cellCenterDPhi[depth], info.eventWeight);
      } else {
        m_HESmallCellDEtaDPhiNegQ[depth]->Fill(
            info.cellCenterDEta[depth], info.cellCenterDPhi[depth], info.eventWeight);
      }
    }
  }
}

void Histograms::IncCutFlow(double weight) {
  cutProgress++;
  m_cutProgress->Fill(cutProgress, weight);
}

void Histograms::IncCutFlow(double weight, int nCuts) {
  for (int i = 0; i < nCuts; i++) {
    IncCutFlow(weight);
  }
}

void Histograms::PlotTrackDisappearance(double TrackP,
                                        double TrackEta,
                                        double TrackPhi,
                                        double minDR,
                                        double minTotalImpactParameter,
                                        double TrackP_dR,
                                        double TrackEta_dR,
                                        double TrackPhi_dR) {
  bool TrackDissapears_dR0p1;
  bool TrackDissapears_IP15;
  bool TrackDissapears_IP10;
  bool TrackDissapears_IP5;
  bool TrackDissapears_IP3;
  bool TrackDissapears_IP2;

  //  std::cout << "Plotting TrackEta: " << TrackEta << " and TrackPhi: " << TrackPhi << std::endl;
  //  std::cout << "minTotalImpactParameter: " << minTotalImpactParameter << std::endl;
  m_histogram_TrackerTrack_EtaPhi_negEta->Fill(TrackEta, TrackPhi);
  m_histogram_TrackerTrack_EtaPhi_posEta->Fill(TrackEta, TrackPhi);

  //  std::cout << "Plotting TrackEta_dR: " << TrackEta_dR << " and TrackPhi_dR: " << TrackPhi_dR << std::endl;
  //  std::cout << "minDR: " << minDR << std::endl;
  m_histogram_TrackerTrack_dRMatched_EtaPhi_negEta->Fill(TrackEta_dR, TrackPhi_dR);
  m_histogram_TrackerTrack_dRMatched_EtaPhi_posEta->Fill(TrackEta_dR, TrackPhi_dR);

  if (TrackP < 50) {
    m_histogram_TrackerTrack_EtaPhi_negEta_pT0to50->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT0to50->Fill(TrackEta, TrackPhi);
  } else if (TrackP < 100) {
    m_histogram_TrackerTrack_EtaPhi_negEta_pT50to100->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT50to100->Fill(TrackEta, TrackPhi);
  } else if (TrackP < 150) {
    m_histogram_TrackerTrack_EtaPhi_negEta_pT100to150->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT100to150->Fill(TrackEta, TrackPhi);
  } else {
    m_histogram_TrackerTrack_EtaPhi_negEta_pT150toInf->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT150toInf->Fill(TrackEta, TrackPhi);
  }

  if (minDR > 0.1) {
    TrackDissapears_dR0p1 = true;
  } else {
    TrackDissapears_dR0p1 = false;
  }

  if (TrackDissapears_dR0p1 == false) {
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_dR0p1->Fill(TrackEta_dR, TrackPhi_dR);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_dR0p1->Fill(TrackEta_dR, TrackPhi_dR);
  }

  if (minTotalImpactParameter > 15) {
    TrackDissapears_IP15 = true;
  } else {
    TrackDissapears_IP15 = false;
  }

  if (minTotalImpactParameter > 5) {
    TrackDissapears_IP5 = true;
  } else {
    TrackDissapears_IP5 = false;
  }

  if (TrackDissapears_IP15 == false) {
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15->Fill(TrackEta, TrackPhi);
    if (TrackP < 50) {
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50->Fill(TrackEta, TrackPhi);
    } else if (TrackP < 100) {
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100->Fill(TrackEta, TrackPhi);
    } else if (TrackP < 150) {
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150->Fill(TrackEta, TrackPhi);
    } else {
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf->Fill(TrackEta, TrackPhi);
    }
  }

  if (minTotalImpactParameter > 10) {
    TrackDissapears_IP10 = true;
  } else {
    TrackDissapears_IP10 = false;
  }

  if (TrackDissapears_IP10 == false) {
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP10->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP10->Fill(TrackEta, TrackPhi);
  }

  if (TrackDissapears_IP5 == false) {
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP5->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP5->Fill(TrackEta, TrackPhi);
  }

  if (minTotalImpactParameter > 3) {
    TrackDissapears_IP3 = true;
  } else {
    TrackDissapears_IP3 = false;
  }

  if (TrackDissapears_IP3 == false) {
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP3->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP3->Fill(TrackEta, TrackPhi);
  }

  if (minTotalImpactParameter > 2) {
    TrackDissapears_IP2 = true;
  } else {
    TrackDissapears_IP2 = false;
  }

  if (TrackDissapears_IP2 == false) {
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP2->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP2->Fill(TrackEta, TrackPhi);
  }
}

void Histograms::PlotCSCHits(const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label) {
  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  if ((TheCSCSegments.isValid()) && (m_histogram_CSCHits_EtaPhi->GetEntries() == 0)) {
    for (CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end();
         iSegment++) {
      //CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
      //if(iDetId.station() != 1) continue;
      DetId TheDetUnitId(iSegment->cscDetId());
      const GeomDetUnit* TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
      m_histogram_CSCHits_EtaPhi->Fill(TheUnit->toGlobal(iSegment->localPosition()).eta(),
                                       TheUnit->toGlobal(iSegment->localPosition()).phi());
    }
  }
}

void Histograms::PlotHCALHits(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label) {
  edm::Handle<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);

  if ((hcalRecHits.isValid())) {
    const HBHERecHitCollection* hbhe = hcalRecHits.product();
    for (HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++) {
      HcalDetId id(hbherechit->detid());
      //std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      //Global3DPoint hbhe_position = hbhe_cell->getPosition();
      //const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
      //m_histogram_CSCHits_EtaPhi->Fill(TheUnit->toGlobal(iSegment->localPosition()).eta(),TheUnit->toGlobal(iSegment->localPosition()).phi());
    }
  }
}

void Histograms::Normalize() {}
