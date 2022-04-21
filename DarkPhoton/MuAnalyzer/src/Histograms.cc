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

void Histograms::book(TFileDirectory histFolder){

  m_eventCount = histFolder.make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
  m_cutProgress = histFolder.make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 10, -.5, 9.5);
  m_MuonTrackMass = histFolder.make<TH1F>("MuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50  , 150  );

  m_failHEEnergy = histFolder.make<TH1F>("failHEEnergy",";HE Energy (GeV);Events",200,0,200);
  m_histogram_TrackerTack_P = histFolder.make<TH1F>("TrackerTrack_P", "", 100, 0, 1000);
  m_CSCMinDR = histFolder.make<TH1F>("CSCMinDR", "", 140, 0, 7);
  m_histogram_TrackerTrack_EtaPhi_negEta = histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT0to50 = histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT0to50", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT50to100 = histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT50to100", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT100to150 = histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT100to150", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_negEta_pT150toInf = histFolder.make<TH2F>("TrackerTrack_EtaPhi_negEta_pT150toInf", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta = histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT0to50 = histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT0to50", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT50to100 = histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT50to100", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT100to150 = histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT100to150", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_EtaPhi_posEta_pT150toInf = histFolder.make<TH2F>("TrackerTrack_EtaPhi_posEta_pT150toInf", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_CSCHits_EtaPhi = histFolder.make<TH2F>("CSCHits_EtaPhi", "; #eta; #phi; Events", 100, -2.4, 2.4, 72, -ROOT::Math::Pi(), ROOT::Math::Pi());
  m_HCALHits_EtaPhi = histFolder.make<TH2F>("HCALHits_EtaPhi", "; i#eta; i#phi; Events", 57, -28.5, 28.5, 72, 0.5, 72.5);

  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT0to50", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT50to100", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT100to150", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP15_pT150toInf", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_IP15_posEta", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT0to50", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT50to100", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT100to150", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP15_pT150toInf", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
		      
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP10 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP10", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP10 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP10", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP5 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP5", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP5 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP5", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP3 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP3", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP3 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP3", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_IP2 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_IP2", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_posEta_IP2 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_IP2", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_TrackerTrack_dRMatched_EtaPhi_negEta = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_dR0p1", "", 100, -2.4, -1.63, 100, -3.2, 3.2);
  m_histogram_TrackerTrack_dRMatched_EtaPhi_posEta = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_dR0p1", "", 100, 1.63, 2.4, 100, -3.2, 3.2);

  m_histogram_TrackerTrackMatched_EtaPhi_posEta_dR0p1 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_posEta_dR0p1", "", 100, 1.63, 2.4, 100, -3.2, 3.2);
  m_histogram_TrackerTrackMatched_EtaPhi_negEta_dR0p1 = histFolder.make<TH2F>("TrackerTrackMatched_EtaPhi_negEta_dR0p1", "", 100, -2.4, -1.63, 100, -3.2, 3.2);

  m_centerInCone = histFolder.make<TH1F>("CenterInCone",";Projected cells in 0.2 DR cone",10,-0.5,9.5);
  m_centerOutOfCone = histFolder.make<TH1F>("CenterOutOfCone",";Projected cells out of 0.2 DR cone",10,-0.5,9.5);
  m_MinDR_MuonHCAL = histFolder.make<TH1F>("MinDR_MuonHCAL", "; Min #Delta R; Events", 140, 0, 0.4);
  m_MinDR_RandomHCAL = histFolder.make<TH1F>("MinDR_RandomHCAL", "; Min #Delta R; Events", 140, 0, 0.4);
  m_HitDepth_MuonHCAL = histFolder.make<TH1F>("HitDepth_MuonHCAL", "; Depth; Events", 40,0,20);
  m_HitDepth_RandomHCAL = histFolder.make<TH1F>("HitDepth_RandomHCAL", "; Depth; Events", 40,0,20);
  m_ConeHits = histFolder.make<TH1F>("ConeHits", "; Hits in Cone; Events", 50,-0.5,49.5);
  m_ConeEnergy = histFolder.make<TH1F>("ConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_RandomConeHits = histFolder.make<TH1F>("RandomConeHits", "; Hits in Cone; Events", 50,-0.5,49.5);
  m_RandomConeEnergy = histFolder.make<TH1F>("RandomConeEnergy", "; Energy in Cone (GeV); Events", 200,0,20);
  m_HitsOverThresh = histFolder.make<TH1F>("HitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);
  m_HEConeEnergy = histFolder.make<TH1F>("HEConeEnergy",";HCAL Energy in 0.4 #Delta R Cone (GeV); Events",100,0,100);
  m_tagHcalIso = histFolder.make<TH1F>("tagHcalIso",";HCAL Energy in 0.4 #Delta R Cone (GeV); Events",100,0,100);
  m_probeHcalIso = histFolder.make<TH1F>("probeHcalIso",";HCAL Energy in 0.4 #Delta R Cone (GeV); Events",100,0,100);
  m_RandomHitsOverThresh = histFolder.make<TH1F>("RandomHitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);

  for(int i=0; i<7; i++)
  {
     std::string name = "Layer"+std::to_string(i+1)+"Spectra";
     std::string rname = "R"+name;
     std::string pairname = "DepthPairSpectra_"+std::to_string(i+1)+std::to_string(i+2);
     std::string rpairname = "R"+pairname;
     std::string pairaxis = "; Depth " + std::to_string(i+1)+" Hit Energy (GeV); Depth "+std::to_string(i+2)+" Hit Energy (GeV); Events";
     m_Layer_Spectra[i] = histFolder.make<TH1F>(name.c_str(), "; Hit Energy (GeV); Events", 400,0,10);
     m_RLayer_Spectra[i] = histFolder.make<TH1F>(rname.c_str(), "; Hit Energy (GeV); Events", 400,0,10);
     if(i<6)
     {
        m_DepthPairSpectra[i] = histFolder.make<TH2F>(pairname.c_str(),pairaxis.c_str(),100,0,10,100,0,10);
        m_RDepthPairSpectra[i] = histFolder.make<TH2F>(rpairname.c_str(),pairaxis.c_str(),100,0,10,100,0,10);
     }
  }
  for(int i=0; i<7; i++)
  {
     std::string name = "Layer"+std::to_string(i+1)+"Eta";
     std::string rname = "R"+name;
     m_Layer_Eta[i] = histFolder.make<TH2F>(name.c_str(), "; Hit Energy (GeV); #eta ; Events", 150,0,10, 50,-2.6,2.6);
     m_RLayer_Eta[i] = histFolder.make<TH2F>(rname.c_str(), "; Hit Energy (GeV); #eta ; Events", 150,0,10, 50,-2.6,2.6);
  }
  m_ValidIDs = histFolder.make<TH1F>("ValidCellIds", "; Valid Cell Ids; Events", 101,-0.5,100.5);
  m_HEBigCellDEtaDPhi = histFolder.make<TH2F>("BigEtaCellCenterDistance","; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_HESmallCellDEtaDPhi = histFolder.make<TH2F>("SmallEtaCellCenterDistance","; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
}

void Histograms::FillHists(EventInfo info)
{
   m_eventCount->Fill(0.5,info.eventWeight);
   m_MuonTrackMass->Fill(info.diMuonMass,info.eventWeight);
   m_CSCMinDR->Fill(info.minDR,info.eventWeight);
   m_HCALHits_EtaPhi->Fill(info.ProbeEta,info.ProbePhi,info.eventWeight);
   m_MinDR_MuonHCAL->Fill(info.minDR_muonHCAL,info.eventWeight);
   m_MinDR_RandomHCAL->Fill(info.minDR_randomHCAL,info.eventWeight);
   m_HitDepth_MuonHCAL->Fill(info.hitDepth_muonHCAL,info.eventWeight);
   m_HitDepth_RandomHCAL->Fill(info.hitDepth_randomHCAL,info.eventWeight);
   m_centerInCone->Fill(info.centerInCone,info.eventWeight);
   m_centerOutOfCone->Fill(info.centerOutOfCone,info.eventWeight);
   m_ConeHits->Fill(info.coneHits,info.eventWeight);
   m_ConeEnergy->Fill(info.coneEnergy,info.eventWeight);
   m_HEConeEnergy->Fill(info.heConeEnergy);
   m_tagHcalIso->Fill(info.tagHcalIso);
   m_probeHcalIso->Fill(info.probeHcalIso);
   m_RandomConeHits->Fill(info.randomConeHits,info.eventWeight);
   m_RandomConeEnergy->Fill(info.randomConeEnergy,info.eventWeight);
   m_HitsOverThresh->Fill(info.hitsOverThresh,info.eventWeight);
   m_RandomHitsOverThresh->Fill(info.randomHitsOverThresh,info.eventWeight);
   for(int i=0;i<7;i++)
   {
      if(i>=info.cellsFound){continue;}
      m_Layer_Spectra[i]->Fill(info.hitEnergies[i],info.eventWeight);
      m_RLayer_Spectra[i]->Fill(info.randHitEnergies[i],info.eventWeight);
      m_Layer_Eta[i]->Fill(info.hitEnergies[i],info.ProbeEta,info.eventWeight);
//      m_RLayer_Eta[i]->Fill(info.randHitEnergies[i],info.ProbeEta,info.eventWeight);
      if(i<6&&((i+1)<info.cellsFound))
      {
         m_DepthPairSpectra[i]->Fill(info.hitEnergies[i],info.hitEnergies[i+1],info.eventWeight);
         m_RDepthPairSpectra[i]->Fill(info.randHitEnergies[i],info.randHitEnergies[i+1],info.eventWeight);
      }
   }
   m_ValidIDs->Fill(info.cellsFound,info.eventWeight);
   if(fabs(info.probeiEta)>20){m_HEBigCellDEtaDPhi->Fill(info.cellCenterDEta,info.cellCenterDPhi,info.eventWeight);}
   else{m_HESmallCellDEtaDPhi->Fill(info.cellCenterDEta,info.cellCenterDPhi,info.eventWeight);}
 
}

void Histograms::IncCutFlow(double weight)
{
   cutProgress++;
   m_cutProgress->Fill(cutProgress,weight);
}

void Histograms::IncCutFlow(double weight,int nCuts)
{
   for(int i=0;i<nCuts;i++)
   {
      IncCutFlow(weight);
   }
}

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

  if ((hcalRecHits.isValid()))
  {
     const HBHERecHitCollection *hbhe = hcalRecHits.product();
     for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
     {
        
	HcalDetId id(hbherechit->detid());
        //std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
	//Global3DPoint hbhe_position = hbhe_cell->getPosition();
	//const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
	//m_histogram_CSCHits_EtaPhi->Fill(TheUnit->toGlobal(iSegment->localPosition()).eta(),TheUnit->toGlobal(iSegment->localPosition()).phi());
     }
  }
	
}

void Histograms::Normalize()
{
}
