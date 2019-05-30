#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "TH1F.h"

#include <iostream>

HCAL::HCAL(){}

void HCAL::CheckHCAL(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label){

  std::cout << "Inside CheckHCAL" << std::endl;

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  const CaloGeometry* caloGeom = TheCALOGeometry.product(); 

  if(!hcalRecHits.isValid()){
//    std::cout << "Could not find HCAL RecHits" << std::endl;
  }else{
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++){
      HcalDetId id(hbherechit->detid());

      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
 //     Global3DPoint hbhe_position = hbhe_cell->getPosition();

/*      std::cout << "hbherechit->energy(): " << hbherechit->energy() << std::endl;
       std::cout << "hbherechit->time(): " << hbherechit->time() << std::endl;

       std::cout << "id.depth(): " << id.depth() << std::endl;
       std::cout << "id.subdet(): " << id.subdet() <<  std::endl;
       std::cout << "hbhe_position.eta(): " << hbhe_position.eta() << std::endl;
       std::cout << "hbhe_position.phi(): " << hbhe_position.phi() << std::endl;
*/
    }
  }
}

double HCAL::MuonMindR(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, double MuonEta, double MuonPhi){

  double minHCALdR = 1000;
  std::cout << "inside MuonMindR" << std::endl;

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

//  std::cout << "check 3" << std::endl;

  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  const CaloGeometry* caloGeom = TheCALOGeometry.product();
//  std::cout << "check 4" << std::endl;

  if(!hcalRecHits.isValid()){
    std::cout << "Could not find HCAL RecHits" << std::endl;
  }else{
//    std::cout << "check 5" << std::endl;
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++){
       HcalDetId id(hbherechit->detid());
//       std::cout << "check 6" << std::endl;

       std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
//       const CaloCellGeometry *hbhe_cell = caloGeom->getGeometry(hbherechit->id());
//       Global3DPoint hbhe_position = caloGeom->getGeometry(hbherechit->id())->getPosition();
       Global3DPoint hbhe_position = hbhe_cell->getPosition();
       double dPhi = fabs(MuonPhi - hbhe_position.phi());
       if(dPhi > ROOT::Math::Pi()) dPhi -= 2*ROOT::Math::Pi();

       if(minHCALdR > sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)))){
	 if(hbherechit->energy()!=0)
	 {
	    minHCALdR = sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)));
	    MuonHitEnergy = hbherechit->energy();
	    MuonHitDepth = id.depth();
	 }
       }
//       hbhe_cell->reset();
    }
  }
  return minHCALdR;
}

double* HCAL::HitsPlots(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, double MuonEta, double MuonPhi, GlobalPoint MuonGlobalPoint,double ConeSize, Histograms myHistograms){ //TH1F* Spectra[7], TH2F* layer_eta[7], TH1F* missinghits, TH1F* missinghitseta){

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
  double RandPhi = MuonPhi + 1.5; 
  if(RandPhi > ROOT::Math::Pi()) RandPhi-=2*ROOT::Math::Pi();
  if(MuonEta<0&&RandPhi<-0.9&&RandPhi>-1.6){RandPhi = RandPhi-3.0;}

  static double Hits[4];
  Hits[0] = 0;
  Hits[1] = 0;
  Hits[2] = 0;
  Hits[3] = 0;

  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  const CaloGeometry* caloGeom = TheCALOGeometry.product();
  const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);
  int MuoniEta,MuoniPhi;
  HcalDetId ClosestCell = (HcalDetId)HEGeom->getClosestCell(MuonGlobalPoint); 
  MuoniEta = ClosestCell.ieta();
  MuoniPhi = ClosestCell.iphi();
  int RandiPhi = MuoniPhi + 20;
  if (RandiPhi>72) RandiPhi -= 72;
  printf("Cell-Muon delta phi: %f\n",caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta);
  printf("Muon iEta: %d, Muon iPhi: %d",MuoniEta,MuoniPhi);
  if(!hcalRecHits.isValid())
  {
    std::cout << "Could not find HCAL RecHits" << std::endl;
  }
  else
  {
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    std::deque < std::tuple <int, int, double> >  MuonHits[7];
    
    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
    {  
       HcalDetId id(hbherechit->detid());
       std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
       Global3DPoint hbhe_position = hbhe_cell->getPosition();
       int HitiEta = id.ieta();
       int HitiPhi = id.iphi();
       
       //double dPhi = fabs(MuonPhi - hbhe_position.phi());
       if((abs(hbhe_position.eta())<1.653)||(abs(hbhe_position.eta())>2.4)){continue;}
       //if(dPhi > ROOT::Math::Pi()) dPhi -= 2*ROOT::Math::Pi();

       //if((ConeSize > sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0))))&&(hbherechit->energy()!=0))
       if(HitiEta==MuoniEta&&HitiPhi==MuoniPhi&&hbherechit->energy()!=0)
       {	 
	 Hits[0]++;
	 Hits[1] += hbherechit->energy();
         if(id.depth()<8&&hbherechit->energy()!=0) {myHistograms.m_Layer_Spectra[id.depth()-1]->Fill(hbherechit->energy());}
	 myHistograms.m_Layer_Eta[id.depth()-1]->Fill(hbherechit->energy(),hbhe_position.eta());
         //if(hbherechit->energy()>Hit_Thresholds[id.depth()-1]){MuonHits[id.depth()-1].push_back(std::make_tuple(hbhe_position.eta(),hbhe_position.phi(),hbherechit->energy()));}
         MuonHits[id.depth()-1].push_back(std::make_tuple(HitiEta,HitiPhi,hbherechit->energy()));
       }
       
       //double RdPhi = fabs(RandPhi - hbhe_position.phi());
       //if(RdPhi > ROOT::Math::Pi()) RdPhi -= 2*ROOT::Math::Pi();

       //if((ConeSize > sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(RdPhi, 2.0))))&&(hbherechit->energy()!=0))
       if(HitiEta==MuoniEta&&HitiPhi==RandiPhi&&hbherechit->energy()!=0)
       {
	 Hits[2]++;
	 Hits[3] += hbherechit->energy();
         if(id.depth()<8) {myHistograms.m_RLayer_Spectra[id.depth()-1]->Fill(hbherechit->energy());}
	 myHistograms.m_RLayer_Eta[id.depth()-1]->Fill(hbherechit->energy(),hbhe_position.eta());
         if(hbherechit->energy()>Hit_Thresholds[id.depth()-1]){MuonHits[id.depth()-1].push_back(std::make_tuple(HitiEta,HitiPhi,hbherechit->energy()));}
       }

//       hbhe_cell->reset();
    }
    for(int j=1; j<5; j++)
    {
       for(auto muonhit = MuonHits[j].begin(); muonhit != MuonHits[j].end(); muonhit++)
       {
          if(std::get<2>(*muonhit)>Hit_Thresholds[j])
	  {
    
             bool missinghit = true;
  	     bool deeperhit = false;
	     double midE = 0;
	     double deepE = 0;
             for(auto deepermuonhit = MuonHits[j+2].begin(); deepermuonhit != MuonHits[j+2].end(); deepermuonhit++)
  	     {
  	        //double dphi = fabs(std::get<1>(*muonhit)-std::get<1>(*deepermuonhit));
  	        //if(dphi > ROOT::Math::Pi()) dphi -= 2*ROOT::Math::Pi();
  	        //if(sqrt( pow(std::get<0>(*muonhit)-std::get<0>(*deepermuonhit),2.0) + pow(dphi, 2.0))<ConeSize)
  	        if(std::get<1>(*muonhit)==std::get<1>(*deepermuonhit))
		{
		   if(std::get<2>(*deepermuonhit)>Hit_Thresholds[j+2])
		   {
  	              deeperhit = true;
		      deepE = std::get<2>(*deepermuonhit);
  	              for(auto middlemuonhit = MuonHits[j+1].begin(); middlemuonhit != MuonHits[j+1].end(); middlemuonhit++)
  	              {
  	                 //double dphi = fabs(std::get<1>(*muonhit)-std::get<1>(*middlemuonhit));
  	                 //if(dphi > ROOT::Math::Pi()) dphi -= 2*ROOT::Math::Pi();
                         //if(sqrt( pow(std::get<0>(*muonhit)-std::get<0>(*middlemuonhit),2.0) + pow(dphi, 2.0))<ConeSize)
  	                 if(std::get<1>(*muonhit)==std::get<1>(*middlemuonhit))
			 {
  	                    if(std::get<2>(*middlemuonhit)>midE) midE = std::get<2>(*middlemuonhit);
  	                    continue;
  	                 }
  	              }
		   }
		   if(midE>Hit_Thresholds[j+1]){missinghit=false;}
  	        }
  	     }
  	     if(deeperhit)
  	     {
  	        //double Muondphi = fabs(std::get<1>(*muonhit)-MuonPhi);
  	        //if(Muondphi > ROOT::Math::Pi()) Muondphi -= 2*ROOT::Math::Pi();
                //double MuonDR = sqrt( pow(std::get<0>(*muonhit)-MuonEta,2.0) + pow(Muondphi, 2.0));
  	        //if(MuonDR<ConeSize)
		double MuonDR = sqrt( pow(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,2.0) + pow(caloGeom->getGeometry(ClosestCell)->phiPos()-MuonPhi,2.0));
  	        if(std::get<1>(*muonhit)==MuoniPhi)
		{
                   if(missinghit) 
                   {	
  	              myHistograms.m_MissingHits->Fill(j+1.5);
  	              myHistograms.m_MissingHitsMap->Fill(std::get<0>(*muonhit),std::get<1>(*muonhit));
  	              if(midE!=0){myHistograms.m_MissingHitsEnergy->Fill(midE);}
                      myHistograms.m_MissingHitsDR->Fill(MuonDR);
		      myHistograms.m_FirstMissing_Spectra->Fill(std::get<2>(*muonhit));
		      myHistograms.m_DeepestMissing_Spectra->Fill(deepE);
  	           }
  	           else 
  	           {
  	              myHistograms.m_MissingHits->Fill(-j-2.5);
                      myHistograms.m_FirstFound_Spectra->Fill(std::get<2>(*muonhit));
		      myHistograms.m_DeepestFound_Spectra->Fill(deepE);
  	           }
  	        }
  	        else
  	        {
                   if(missinghit) 
                   {	
  	              myHistograms.m_RMissingHits->Fill(j+1.5);
  	              myHistograms.m_RMissingHitsMap->Fill(std::get<0>(*muonhit),std::get<1>(*muonhit));
  	           }
  	           else {myHistograms.m_RMissingHits->Fill(-j-2.5);}
                }
  	     }
          }
       }
    }
  }
  return Hits;
}


