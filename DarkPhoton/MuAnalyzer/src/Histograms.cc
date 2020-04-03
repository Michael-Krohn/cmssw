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
#include <iostream>

Histograms::Histograms(){}

void Histograms::book(edm::Service<TFileService> fs){

  m_eventCount = fs->make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
  m_cutProgress = fs->make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 10, -.5, 9.5);

  m_MuonTrackMass = fs->make<TH1F>("MuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50  , 150  );

  m_GENMuonTrackdR = fs->make<TH1F>("DeltaR_GENMuonTrack", ";DeltaR(GENMuon,Track) ;Events", 30,  0.  , 0.2  );
  m_nMuonTrackCand = fs->make<TH1F>("nMuonTrackCand", "; # of MuonTrack Candidates ;Events", 10, -0.5  , 9.5  );
  m_nTracksPairedPerMuon = fs->make<TH1F>("nTracksPairedPerMuon", "; # of MuonTrack Candidates ;Events", 10, -0.5  , 9.5  );
  m_nMuonsPairedPerEvent = fs->make<TH1F>("nMuonsPairedPerEvent", "; # of MuonTrack Candidates ;Events", 10, -0.5  , 9.5  );
  m_nTracksNoMuon = fs->make<TH1F>("nTrackNoMuon", "; # of TrackerTrackMatched Candidates ;Events", 8, -0.5  , 7.5  );
  
  m_histogram_TrackerTack_P = fs->make<TH1F>("TrackerTrack_P", "", 100, 0, 1000);
  m_histogram_MuonTrack_P = fs->make<TH1F>("MuonTrack_P", "", 140, 0, 700);
    
  m_histogram_differenceVertex_xy = fs->make<TH2F>("differenceVertex_xy", "", 300, -1., 1., 300, -1., 1.);
  m_histogram_misshitVertex_xy = fs->make<TH2F>("misshitVertex_xy", "", 300, -1., 1., 300, -1., 1.);
  m_histogram_differenceVertex_z = fs->make<TH1F>("differenceVertex_z", "", 60, -20., 20.);
  m_histogram_misshitVertex_z = fs->make<TH1F>("misshitVertex_z", "", 60, -20., 20.);
  m_adjacentfailvertex_z = fs->make<TH1F>("adjacentfailvertex_z", "", 60, -20., 20.);

  m_MinTransverseImpactParameter = fs->make<TH1F>("MinTransverseImpactParameter", "", 300, -50, 50);
  m_MinLongiudinalImpactParameter = fs->make<TH1F>("MinLongiudinalImpactParameter", "", 300, -50, 50);
  m_MinTotalImpactParameter = fs->make<TH1F>("MinTotalImpactParameter", "", 300, 0, 1000);
  m_MinDR = fs->make<TH1F>("MinDR", "", 140, 0, 7);

  m_MinDR_Muon = fs->make<TH1F>("MinDR_Muon", "", 140, 0, 7);
  m_MinTotalImpactParameterMuon = fs->make<TH1F>("MinTotalImpactParameterMuon", "", 300, 0, 300);

  m_histogram_TrackerTrack_EtaPhi_negEta = fs->make<TH2F>("TrackerTrack_EtaPhi_negEta", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT0to50 = fs->make<TH2F>("TrackerTrack_EtaPhi_negEta_pT0to50", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT50to100 = fs->make<TH2F>("TrackerTrack_EtaPhi_negEta_pT50to100", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT100to150 = fs->make<TH2F>("TrackerTrack_EtaPhi_negEta_pT100to150", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT150toInf = fs->make<TH2F>("TrackerTrack_EtaPhi_negEta_pT150toInf", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta = fs->make<TH2F>("TrackerTrack_EtaPhi_posEta", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT0to50 = fs->make<TH2F>("TrackerTrack_EtaPhi_posEta_pT0to50", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT50to100 = fs->make<TH2F>("TrackerTrack_EtaPhi_posEta_pT50to100", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT100to150 = fs->make<TH2F>("TrackerTrack_EtaPhi_posEta_pT100to150", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT150toInf = fs->make<TH2F>("TrackerTrack_EtaPhi_posEta_pT150toInf", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_CSCHits_EtaPhi = fs->make<TH2F>("CSCHits_EtaPhi", "; #eta; #phi; Events", 100, -2.4, 2.4, 72, -ROOT::Math::Pi(), ROOT::Math::Pi());
  m_histogram_HCALHits_EtaPhi = fs->make<TH2F>("HCALHits_EtaPhi", "; i#eta; i#phi; Events", 57, -28.5, 28.5, 72, 0.5, 72.5);
  m_histogram_BlankHCALHits_EtaPhi = fs->make<TH2F>("BlankHCALHits_EtaPhi", "; i#eta; i#phi; Events", 57, -28.5, 28.5, 72, 0.5, 72.5);

  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_IP15_posEta", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
		      
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP10 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP10", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP10 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP10", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP5 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP5", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP5 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP5", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP3 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP3", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP3 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP3", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP2 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP2", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP2 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP2", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_TrackerTrack_dRMatched_EtaPhi_negEta = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_dR0p1", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_dRMatched_EtaPhi_posEta = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_dR0p1", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_TrackerTrackMatched_EtaPhi_posEta_dR0p1 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_dR0p1", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_dR0p1 = fs->make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_dR0p1", "", 100, -2.4, -1.63, 100, -3.2, 3.2);

  m_MinDR_MuonHCAL = fs->make<TH1F>("MinDR_MuonHCAL", "; Min #Delta R; Events", 140, 0, 0.4);
  m_MinDR_RandomHCAL = fs->make<TH1F>("MinDR_RandomHCAL", "; Min #Delta R; Events", 140, 0, 0.4);
  m_HitEnergy_MinDR_MuonHCAL = fs->make<TH1F>("HitEnergy_MinDR_MuonHCAL", "; Energy; Events", 200, 0, 20);
  m_HitEnergy_RandomHCAL = fs->make<TH1F>("HitEnergy_MinDR_RandomHCAL", "; Energy; Events", 200, 0, 20);
  m_HitDepth_MuonHCAL = fs->make<TH1F>("HitDepth_MuonHCAL", "; Depth; Events", 40,0,20);
  m_HitDepth_RandomHCAL = fs->make<TH1F>("HitDepth_RandomHCAL", "; Depth; Events", 40,0,20);
  m_ConeHits = fs->make<TH1F>("ConeHits", "; Hits in Cone; Events", 50,-0.5,49.5);
  m_ConeEnergy = fs->make<TH1F>("ConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,200);
  m_MissThreshConeEnergy = fs->make<TH1F>("MissConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,200);
  m_RandomConeHits = fs->make<TH1F>("RandomConeHits", "; Hits in Cone; Events", 50,-0.5,49.5);
  m_RandomConeEnergy = fs->make<TH1F>("RandomConeEnergy", "; Energy in Cone (GeV); Events", 200,0,200);
  m_HitsOverThresh = fs->make<TH1F>("HitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);
  m_RandomHitsOverThresh = fs->make<TH1F>("RandomHitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);

  for(int i=0; i<7; i++)
  {
     std::string name = "Layer"+std::to_string(i+1)+"Spectra";
     std::string rname = "R"+name;
     std::string pairname = "DepthPairSpectra_"+std::to_string(i+1)+std::to_string(i+2);
     std::string rpairname = "R"+pairname;
     std::string pairaxis = "; Depth " + std::to_string(i+1)+" Hit Energy (GeV); Depth "+std::to_string(i+2)+" Hit Energy (GeV); Events";
     m_Layer_Spectra[i] = fs->make<TH1F>(name.c_str(), "; Hit Energy (GeV); Events", 400,0,10);
     m_RLayer_Spectra[i] = fs->make<TH1F>(rname.c_str(), "; Hit Energy (GeV); Events", 400,0,10);
     if(i<6)
     {
        m_DepthPairSpectra[i] = fs->make<TH2F>(pairname.c_str(),pairaxis.c_str(),100,0,10,100,0,10);
        m_RDepthPairSpectra[i] = fs->make<TH2F>(rpairname.c_str(),pairaxis.c_str(),100,0,10,100,0,10);
     }
  }
  for(int i=0; i<7; i++)
  {
     std::string name = "Layer"+std::to_string(i+1)+"Eta";
     std::string rname = "R"+name;
     m_Layer_Eta[i] = fs->make<TH2F>(name.c_str(), "; Hit Energy (GeV); #eta ; Events", 150,0,10, 50,-2.6,2.6);
     m_RLayer_Eta[i] = fs->make<TH2F>(rname.c_str(), "; Hit Energy (GeV); #eta ; Events", 150,0,10, 50,-2.6,2.6);
  }
  m_ValidIDs = fs->make<TH1F>("ValidCellIds", "; Valid Cell Ids; Events", 101,-0.5,100.5);
  m_Missing_ValidIDs = fs->make<TH1F>("MissingHitValidCellIds", "; Valid Cell Ids; Events", 101, -0.5, 100.5);
  m_MissingHits = fs->make<TH1F>("MissingHits", "; Hit Depth; Events", 16,-8.5,7.5);
  m_RMissingHits = fs->make<TH1F>("RMissingHits", "; Hit Depth; Events", 16,-8.5,7.5);
  m_MissingHitsMap = fs->make<TH2F>("MissingHitsMap", "; i#eta; i#phi; Events", 58, -28.5, 28.5, 73, 0.5, 72.5); 
  m_RMissingHitsMap = fs->make<TH2F>("RMissingHitsMap", ";i#eta; i#phi; Events", 58, -28.5, 28.5, 73, 0.5, 72.5);
  m_MissingHitsEnergy = fs->make<TH1F>("MissingHitsEnergy", "; Hit Energy (GeV); Events", 50,0,0.2);
  m_MissingHitsDR = fs->make<TH1F>("MinDR_MissingHits", "; Min #Delta R; Events", 70, 0, 0.2);
  m_BlankHitsDR = fs->make<TH1F>("MinDR_BlankHits", "; #Delta R; Events", 70, 0, 0.2);
  m_TrackHCALDR_GoodHits = fs->make<TH1F>("MinDRTrack_GoodHits", "; #Delta R; Events", 70, 0, 0.2);
  m_TrackHCALDR_MissHit = fs->make<TH1F>("MinDRTrack_MissHit", "; #Delta R; Events", 70, 0, 0.2);
  m_BlankDepth = fs->make<TH1F>("Blank_Depth", "; Hit Depth; Events", 7,0.5,7.5);
  m_4BlankDepth = fs->make<TH1F>("TwoMiss_Depth", "; Hit Depth; Events", 7,0.5,7.5);
  m_BlankCellDetaDphiPosEta = fs->make<TH2F>("LostHitsDetaDphiPosEta", ";#eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BlankCellSmallDetaDphiPosEta = fs->make<TH2F>("SmallLostHitsDetaDphiPosEta", ";#eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BlankCellDetaDphiNegEta = fs->make<TH2F>("LostHitsDetaDphiNegEta", ";#eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BlankCellSmallDetaDphiNegEta = fs->make<TH2F>("SmallLostHitsDetaDphiNegEta", ";#eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_MuonEtaDist = fs->make<TH1F>("MuonEtaDist",";#eta; Events", 100,-10,10);
  m_TrackPt = fs->make<TH1F>("SelectedTrackPt",";#Pt (GeV); Events", 200, 0, 100);
  m_MissHitTrackPt = fs->make<TH1F>("MissHitTrackPt",";#Pt (GeV); Events", 200, 0, 100);
  m_DBremLocation = fs->make<TH2F>("DarkBremLocation",";Z;#rho;Events",150,0,600,150,0,600);
  m_DeflectingAngle = fs->make<TH1F>("DeflectingAngle",";Muon Deflection Angle;Events",100,0,3.16);
  m_initialMuE = fs->make<TH1F>("MuInitialE",";Energy (GeV);Events",100,0,200);
  m_finalMuE = fs->make<TH1F>("MuFinalE",";Energy (GeV);Events",100,0,200);
  m_dphoEnergy = fs->make<TH1F>("DarkPhotonEnergy",";Energy (GeV);Events",100,0,200);
  m_NThreshCut = fs->make<TH1F>("NPassAdjacentCut","; Fail cut in bin zero, pass in bin 1; Events", 2,-0.5,1.5);
}

/*void Histograms::IncCutFlow()
{
   cutProgress++;
   m_cutProgress->fill(cutProgress);
}*/

void Histograms::PlotTrackDisappearance(double TrackP, double TrackEta, double TrackPhi, double minDR, double minTotalImpactParameter, double TrackP_dR, double TrackEta_dR, double TrackPhi_dR){

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


  if(TrackP < 50){
    m_histogram_TrackerTrack_EtaPhi_negEta_pT0to50->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT0to50->Fill(TrackEta, TrackPhi);
  }else if(TrackP < 100){
    m_histogram_TrackerTrack_EtaPhi_negEta_pT50to100->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT50to100->Fill(TrackEta, TrackPhi);
  }else if(TrackP < 150){
    m_histogram_TrackerTrack_EtaPhi_negEta_pT100to150->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT100to150->Fill(TrackEta, TrackPhi);
  }else{
    m_histogram_TrackerTrack_EtaPhi_negEta_pT150toInf->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrack_EtaPhi_posEta_pT150toInf->Fill(TrackEta, TrackPhi);
  }

  if(minDR > 0.1){
    TrackDissapears_dR0p1 = true;
  }else{
    TrackDissapears_dR0p1 = false;
  }

  if(TrackDissapears_dR0p1 == false){
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_dR0p1->Fill(TrackEta_dR, TrackPhi_dR);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_dR0p1->Fill(TrackEta_dR, TrackPhi_dR);
  }

  if(minTotalImpactParameter > 15){
    TrackDissapears_IP15 = true;
  }else{
    TrackDissapears_IP15 = false;
  }

  if(minTotalImpactParameter > 5){
    TrackDissapears_IP5 = true;
  }else{
    TrackDissapears_IP5 = false;
  }

  if(TrackDissapears_IP15 == false){
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15->Fill(TrackEta, TrackPhi);
    if(TrackP < 50){
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50->Fill(TrackEta, TrackPhi);
    }else if(TrackP < 100){
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100->Fill(TrackEta, TrackPhi);
    }else if(TrackP < 150){
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150->Fill(TrackEta, TrackPhi);
    }else{
      m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf->Fill(TrackEta, TrackPhi);
      m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf->Fill(TrackEta, TrackPhi);
    }
  }

  if(minTotalImpactParameter > 10){
    TrackDissapears_IP10 = true;
  }else{
    TrackDissapears_IP10 = false;
  }
				  
  if(TrackDissapears_IP10== false){
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP10->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP10->Fill(TrackEta, TrackPhi);
  }
			    
  if(TrackDissapears_IP5== false){
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP5->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP5->Fill(TrackEta, TrackPhi);
  }

  if(minTotalImpactParameter > 3){
    TrackDissapears_IP3 = true;
  }else{
    TrackDissapears_IP3 = false;
  }

  if(TrackDissapears_IP3== false){
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP3->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP3->Fill(TrackEta, TrackPhi);
  }

  if(minTotalImpactParameter > 2){
    TrackDissapears_IP2 = true;
  }else{
    TrackDissapears_IP2 = false;
  }

  if(TrackDissapears_IP2== false){
    m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP2->Fill(TrackEta, TrackPhi);
    m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP2->Fill(TrackEta, TrackPhi);
  }
			      

}

void Histograms::PlotCSCHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label)
{

  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  if ((TheCSCSegments.isValid())&&(m_histogram_CSCHits_EtaPhi->GetEntries()==0))
  {
     for(CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end(); iSegment++)
     {
        //CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
	//if(iDetId.station() != 1) continue;
        DetId TheDetUnitId(iSegment->cscDetId());
	const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
	m_histogram_CSCHits_EtaPhi->Fill(TheUnit->toGlobal(iSegment->localPosition()).eta(),TheUnit->toGlobal(iSegment->localPosition()).phi());
     }
  }
	
}

void Histograms::PlotHCALHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label)
{

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  const CaloGeometry* caloGeom = TheCALOGeometry.product();

  if ((hcalRecHits.isValid()))
  {
     const HBHERecHitCollection *hbhe = hcalRecHits.product();
     for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
     {
        
	HcalDetId id(hbherechit->detid());
        //std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
	//Global3DPoint hbhe_position = hbhe_cell->getPosition();
	m_histogram_HCALHits_EtaPhi->Fill(id.ieta(),id.iphi());
        
	//const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
	//m_histogram_CSCHits_EtaPhi->Fill(TheUnit->toGlobal(iSegment->localPosition()).eta(),TheUnit->toGlobal(iSegment->localPosition()).phi());
     }
  }
	
}

void Histograms::Normalize()
{
 /* if(m_histogram_CSCHits_EtaPhi->Integral()!=0)
  {
     double norm = 1/m_histogram_CSCHits_EtaPhi->Integral();
     m_histogram_CSCHits_EtaPhi->Scale(norm);
  }*/
}
