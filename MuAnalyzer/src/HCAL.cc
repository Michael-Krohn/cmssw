#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "TH1F.h"

#include <algorithm>
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

double HCAL::MuonMindR(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, GlobalPoint MuonGlobalPoint){

  double MuonEta = MuonGlobalPoint.eta();
  double MuonPhi = MuonGlobalPoint.phi();
  double minHCALdR = 1000;
  std::cout << "inside MuonMindR" << std::endl;

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);


  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  const CaloGeometry* caloGeom = TheCALOGeometry.product();

  if(!hcalRecHits.isValid()){
    std::cout << "Could not find HCAL RecHits" << std::endl;
  }else{
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++){
       HcalDetId id(hbherechit->detid());
       std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
       Global3DPoint hbhe_position = hbhe_cell->getPosition();
       double dPhi = fabs(MuonPhi - hbhe_position.phi());
       if(dPhi > ROOT::Math::Pi()) dPhi -= 2*ROOT::Math::Pi();
       if(fabs(id.ieta())<18){continue;}

       if(minHCALdR > sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)))){
	 if(hbherechit->energy()!=0)
	 {
	    minHCALdR = sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)));
	    MuonMinDr = minHCALdR;
	    MuonHitEnergy = hbherechit->energy();
	    MuonHitDepth = id.depth();
	    minHCALDetId = id;
	 }
       }
//       hbhe_cell->reset();
    }
  }
  return minHCALdR;
}

void HCAL::GetConeIDs(const HcalTopology* theHBHETopology, HcalDetId *MuonAlignedCells, HcalDetId ClosestCell, const int Ndepths, const int CellsPerDepth){ 

  int startdepth = ClosestCell.depth();
  HcalDetId IteratingId = ClosestCell;
  for(int i=startdepth;i>0;i--)
  {
     MuonAlignedCells[(i-1)*CellsPerDepth]=IteratingId;
     theHBHETopology->decrementDepth(IteratingId);
  }
  IteratingId = ClosestCell;
  for(int i=startdepth+1;i<=Ndepths;i++)
  {
     if(!theHBHETopology->validHcal(MuonAlignedCells[CellsPerDepth*(i-2)])){continue;}
     theHBHETopology->incrementDepth(IteratingId);
     MuonAlignedCells[(i-1)*CellsPerDepth]=IteratingId;
  }
  for(int i=1;i<=Ndepths;i++)
  {
     if(!theHBHETopology->validHcal(MuonAlignedCells[CellsPerDepth*(i-1)])){continue;}
     HcalDetId incIEta[2];
     HcalDetId decIEta[2];
     theHBHETopology->incIEta(MuonAlignedCells[CellsPerDepth*(i-1)],incIEta);
     theHBHETopology->decIEta(MuonAlignedCells[CellsPerDepth*(i-1)],decIEta);
     HcalDetId incIPhi;
     HcalDetId decIPhi;
     if(theHBHETopology->incIPhi(MuonAlignedCells[CellsPerDepth*(i-1)],incIPhi)){MuonAlignedCells[CellsPerDepth*(i-1)+5]=incIPhi;} 
     if(theHBHETopology->decIPhi(MuonAlignedCells[CellsPerDepth*(i-1)],decIPhi)){MuonAlignedCells[CellsPerDepth*(i-1)+6]=decIPhi;}
     MuonAlignedCells[CellsPerDepth*(i-1)+1]=incIEta[0];
     MuonAlignedCells[CellsPerDepth*(i-1)+2]=incIEta[1];
     MuonAlignedCells[CellsPerDepth*(i-1)+3]=decIEta[0];
     MuonAlignedCells[CellsPerDepth*(i-1)+4]=decIEta[1];
  }
  return;
}

void HCAL::GetCornerIDs(const HcalTopology* theHBHETopology, HcalDetId *CornerCells, HcalDetId ClosestCell, const int Ndepths){ 
  const int CellsPerDepth = 8;
  int startdepth = ClosestCell.depth();
  HcalDetId IteratingId = ClosestCell;
  HcalDetId CenterIds[Ndepths];
  for(int i=startdepth;i>0;i--)
  {
     CenterIds[i-1]=IteratingId;
     theHBHETopology->decrementDepth(IteratingId);
  }
  IteratingId = ClosestCell;
  for(int i=startdepth+1;i<=Ndepths;i++)
  {
     if(!theHBHETopology->validHcal(CenterIds[i-2])){continue;}
     theHBHETopology->incrementDepth(IteratingId);
     CenterIds[i-1]=IteratingId;
  }
  for(int i=1;i<=Ndepths;i++)
  {
     if(!theHBHETopology->validHcal(CenterIds[i-1])){continue;}
     HcalDetId incIPhi;
     HcalDetId decIPhi;
     theHBHETopology->incIPhi(CenterIds[i-1],incIPhi);
     theHBHETopology->decIPhi(CenterIds[i-1],decIPhi);
     HcalDetId incIEtaincIPhi[2];
     HcalDetId decIEtadecIPhi[2];
     HcalDetId incIEtadecIPhi[2];
     HcalDetId decIEtaincIPhi[2];
     theHBHETopology->incIEta(incIPhi,incIEtaincIPhi);
     theHBHETopology->decIEta(incIPhi,decIEtaincIPhi);
     theHBHETopology->incIEta(decIPhi,incIEtadecIPhi);
     theHBHETopology->decIEta(decIPhi,decIEtadecIPhi);
     CornerCells[CellsPerDepth*(i-1)]=incIEtaincIPhi[0];
     CornerCells[CellsPerDepth*(i-1)+1]=incIEtaincIPhi[1];
     CornerCells[CellsPerDepth*(i-1)+2]=decIEtaincIPhi[0];
     CornerCells[CellsPerDepth*(i-1)+3]=decIEtaincIPhi[1];
     CornerCells[CellsPerDepth*(i-1)+4]=incIEtadecIPhi[0];
     CornerCells[CellsPerDepth*(i-1)+5]=incIEtadecIPhi[1];
     CornerCells[CellsPerDepth*(i-1)+6]=decIEtadecIPhi[0];
     CornerCells[CellsPerDepth*(i-1)+7]=decIEtadecIPhi[1];
  }
  return;
}

void HCAL::HitsPlots(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, GlobalPoint MuonGlobalPoint, GlobalPoint RandGlobalPoint, bool GoodRand, Histograms myHistograms, double CSCMuonMinDr){ 

  double MuonEta = MuonGlobalPoint.eta();
  double MuonPhi = MuonGlobalPoint.phi();
  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
  static double Hits[4];
  Hits[0] = 0;
  Hits[1] = 0;
  Hits[2] = 0;
  Hits[3] = 0;

  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  
  edm::ESHandle<HcalTopology> htopo;
  iSetup.get<HcalRecNumberingRecord>().get(htopo);
  const HcalTopology* theHBHETopology = htopo.product();
  
  const CaloGeometry* caloGeom = TheCALOGeometry.product();
  const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);
  int MuoniEta,MuoniPhi;
  HcalDetId ClosestCell = (HcalDetId)HEGeom->getClosestCell(MuonGlobalPoint); 
  MuoniEta = ClosestCell.ieta();
  MuoniPhi = ClosestCell.iphi();
  if(fabs(MuoniEta)<20){return;}
  if(MuoniEta<-16&&MuoniPhi>53&&MuoniPhi<63){return;}

  const int Ndepths = 7;
  const int CellsPerDepth = 7;
  HcalDetId MuonAlignedCells[CellsPerDepth*Ndepths];
  GetConeIDs(theHBHETopology,MuonAlignedCells,ClosestCell,Ndepths,CellsPerDepth);
  HcalDetId CornerAlignedCells[8*Ndepths];
  GetCornerIDs(theHBHETopology,CornerAlignedCells,ClosestCell,Ndepths);
  int ValidIdCount = 0;
  for(int i=0;i<CellsPerDepth*Ndepths;i++){if(theHBHETopology->validHcal(MuonAlignedCells[i])){ValidIdCount++;}}
  myHistograms.m_ValidIDs->Fill(ValidIdCount);

  HcalDetId RandClosestCell = (HcalDetId)HEGeom->getClosestCell(RandGlobalPoint);
  //int RandiPhi = MuoniPhi + 20;
  //if (RandiPhi>72) RandiPhi -= 72;
  //if (MuoniEta<-16&&RandiPhi>53&&RandiPhi<63) RandiPhi = MuoniPhi-20;
  //HcalDetId RandClosestCell(HcalEndcap,MuoniEta,RandiPhi,1);
  int RandiPhi = RandClosestCell.iphi();
  HcalDetId RandAlignedCells[CellsPerDepth*Ndepths];
  HcalDetId RandCornerAlignedCells[8*Ndepths];
  GetCornerIDs(theHBHETopology,RandCornerAlignedCells,RandClosestCell,Ndepths);

  GetConeIDs(theHBHETopology,RandAlignedCells,RandClosestCell,Ndepths,CellsPerDepth);

  if(!hcalRecHits.isValid())
  {
    std::cout << "Could not find HCAL RecHits" << std::endl;
  }
  else
  {
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    std::deque < std::tuple <int, int, double> >  MuonHits[7];
    double layerenergies[7],rlayerenergies[7];
    for(int i=0;i<7;i++)
    {
       layerenergies[i]=0;
       rlayerenergies[i]=0;
    }

    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
    {  
       HcalDetId id(hbherechit->detid());
       std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
       Global3DPoint hbhe_position = hbhe_cell->getPosition();
       
       HcalDetId *match = std::find(std::begin(MuonAlignedCells), std::end(MuonAlignedCells), id); 
       HcalDetId *cornermatch = std::find(std::begin(CornerAlignedCells), std::end(CornerAlignedCells), id);       
       int HitiEta = id.ieta();
       //int HitiPhi = id.iphi();
        
       if(fabs(HitiEta)<19){continue;}
       if(hbherechit->energy()!=0&&((match!=std::end(MuonAlignedCells))||(cornermatch!=std::end(CornerAlignedCells))))
       //if(hbherechit->energy()!=0&&(match!=std::end(MuonAlignedCells)))
       {	 
	 //Hits[0]++;
	 Hits[1] += hbherechit->energy();
         if(id.depth()<8&&hbherechit->energy()!=0) 
	 {
	    //if(hbherechit->energy()>layerenergies[id.depth()-1]){layerenergies[id.depth()-1]=hbherechit->energy();}
	    layerenergies[id.depth()-1]+=hbherechit->energy();
	 }
	 myHistograms.m_Layer_Eta[id.depth()-1]->Fill(hbherechit->energy(),hbhe_position.eta());
	 myHistograms.m_HitDepth_MuonHCAL->Fill(id.depth());
       }
       
       HcalDetId *randmatch = std::find(std::begin(RandAlignedCells), std::end(RandAlignedCells), id); 
       HcalDetId *randcornermatch = std::find(std::begin(RandCornerAlignedCells), std::end(RandCornerAlignedCells), id);       
       if(hbherechit->energy()!=0&&((randmatch!=std::end(RandAlignedCells))||(randcornermatch!=std::end(RandCornerAlignedCells)))&&GoodRand)
       {
	 Hits[2]++;
	 Hits[3] += hbherechit->energy();
	 myHistograms.m_HitDepth_RandomHCAL->Fill(id.depth());
         if(id.depth()<8) 
	 {
	    //if(hbherechit->energy()>rlayerenergies[id.depth()-1]){rlayerenergies[id.depth()-1]=hbherechit->energy();}
	    rlayerenergies[id.depth()-1]+=hbherechit->energy();
	 }
         if(hbherechit->energy()<0){printf("Negative energy. Energy is: %f\n",hbherechit->energy());}
	 myHistograms.m_RLayer_Eta[id.depth()-1]->Fill(hbherechit->energy(),hbhe_position.eta());
       }

//       hbhe_cell->reset();
    }
    int hitsoverthresh=0;
    for(int i=0;i<7;i++) 
    {
       myHistograms.m_Layer_Spectra[i]->Fill(layerenergies[i]);
       if(layerenergies[i]!=0)
       {
          if(i<6){Hits[0]++;}
          MuonHits[i].push_back(std::make_tuple(MuoniEta,MuoniPhi,layerenergies[i]));
       }
       if(GoodRand)
       {
          myHistograms.m_RLayer_Spectra[i]->Fill(rlayerenergies[i]);
          if(rlayerenergies[i]>Hit_Thresholds[i]){hitsoverthresh++;}
          if(rlayerenergies[i]!=0)
          {
             MuonHits[i].push_back(std::make_tuple(MuoniEta,RandiPhi,rlayerenergies[i]));
          }
       }
    }   
    for(int i=0;i<7;i++)
    {
       if(layerenergies[i]==0)
       {
          if(i<6)
	  {
	     if(Hits[0]==4){myHistograms.m_4BlankDepth->Fill(i+1);}
             if(Hits[0]==5){myHistograms.m_BlankDepth->Fill(i+1);}
          }
       }
    }
    if(GoodRand){myHistograms.m_HitsOverThreshold->Fill(hitsoverthresh);}

    myHistograms.m_ConeHits->Fill(Hits[0]);
    if(Hits[0]==0)
    {
       //myHistograms.m_histogram_BlankHCALHits_EtaPhi->Fill(MuoniEta,MuoniPhi);
       double bdphi = fabs(caloGeom->getGeometry(ClosestCell)->phiPos()-MuonPhi);
       if(bdphi>ROOT::Math::Pi()) bdphi -= 2*ROOT::Math::Pi();
       myHistograms.m_BlankHitsDR->Fill(sqrt( pow(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,2.0) + pow(bdphi,2.0)));
       myHistograms.m_TrackHCALDR_BlankHits->Fill(CSCMuonMinDr);
    }
    if(Hits[0]==4){myHistograms.m_histogram_BlankHCALHits_EtaPhi->Fill(MuoniEta,MuoniPhi);}
    if(Hits[0]==6){myHistograms.m_TrackHCALDR_SixHit->Fill(MuonMinDr);}
    myHistograms.m_ConeEnergy->Fill(Hits[1]);
    if(GoodRand)
    {
       myHistograms.m_RandomConeHits->Fill(Hits[2]);
       myHistograms.m_RandomConeEnergy->Fill(Hits[3]);
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
                double mdphi = fabs(caloGeom->getGeometry(ClosestCell)->phiPos()-MuonPhi);
                if(mdphi>ROOT::Math::Pi()) mdphi -= 2*ROOT::Math::Pi();
		double MuonDR = sqrt( pow(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,2.0) + pow(mdphi,2.0));
  	        //if(MuonDR<0.3)
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
  	        else if(GoodRand)
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
  return;
}


