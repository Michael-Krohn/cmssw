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
    
  m_histogram_differenceVertex_xy = fs->make<TH2F>("differenceVertex_xy", "", 300, -3., 3., 300, -3., 3.);
  m_histogram_differenceVertex_z = fs->make<TH1F>("differenceVertex_z", "", 60, -10., 10.);

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

  m_histogram_CSCHits_EtaPhi = fs->make<TH2F>("CSCHits_EtaPhi", "", 100, -2.4, 2.4, 100, -3.2, 3.2);

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

  if (TheCSCSegments.isValid())
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

void Histograms::Normalize()
{
 /* if(m_histogram_CSCHits_EtaPhi->Integral()!=0)
  {
     double norm = 1/m_histogram_CSCHits_EtaPhi->Integral();
     m_histogram_CSCHits_EtaPhi->Scale(norm);
  }*/
}
