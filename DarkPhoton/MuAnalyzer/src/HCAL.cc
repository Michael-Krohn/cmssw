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

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "DarkPhoton/MuAnalyzer/interface/MCHistograms.h"
#include "TH1F.h"

#include "Math/VectorUtil.h"
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
       if(fabs(id.ieta())<19){continue;}

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

void HCAL::GetCenterCells(const HcalTopology* theHBHETopology, HcalDetId *TrackAlignedCells, HcalDetId ClosestCell, const int Ndepths, const int CellsPerDepth)
{
   int startdepth = ClosestCell.depth();
   HcalDetId IteratingId = ClosestCell;
   for(int i=startdepth;i>0;i--)
   {
      TrackAlignedCells[(i-1)*CellsPerDepth]=IteratingId;
      theHBHETopology->decrementDepth(IteratingId);
   }
   IteratingId = ClosestCell;
   for(int i=startdepth+1;i<=Ndepths;i++)
   {
      if(!theHBHETopology->validHcal(TrackAlignedCells[CellsPerDepth*(i-2)])){continue;}
      theHBHETopology->incrementDepth(IteratingId);
      TrackAlignedCells[(i-1)*CellsPerDepth]=IteratingId;
   }

   return;
}

void HCAL::GetConeIDs(const HcalTopology* theHBHETopology, HcalDetId *TrackAlignedCells, HcalDetId ClosestCell, const int Ndepths, const int CellsPerDepth){ 

  GetCenterCells(theHBHETopology, TrackAlignedCells, ClosestCell, Ndepths, CellsPerDepth);

  for(int i=1;i<=Ndepths;i++)
  {
     if(!theHBHETopology->validHcal(TrackAlignedCells[CellsPerDepth*(i-1)])){continue;}
     HcalDetId incIEta[2];
     HcalDetId decIEta[2];
     theHBHETopology->incIEta(TrackAlignedCells[CellsPerDepth*(i-1)],incIEta);
     theHBHETopology->decIEta(TrackAlignedCells[CellsPerDepth*(i-1)],decIEta);
     HcalDetId incIPhi;
     HcalDetId decIPhi;
     if(theHBHETopology->incIPhi(TrackAlignedCells[CellsPerDepth*(i-1)],incIPhi)){TrackAlignedCells[CellsPerDepth*(i-1)+5]=incIPhi;} 
     if(theHBHETopology->decIPhi(TrackAlignedCells[CellsPerDepth*(i-1)],decIPhi)){TrackAlignedCells[CellsPerDepth*(i-1)+6]=decIPhi;}
     TrackAlignedCells[CellsPerDepth*(i-1)+1]=incIEta[0];
     TrackAlignedCells[CellsPerDepth*(i-1)+2]=incIEta[1];
     TrackAlignedCells[CellsPerDepth*(i-1)+3]=decIEta[0];
     TrackAlignedCells[CellsPerDepth*(i-1)+4]=decIEta[1];
  }
  return;
}

void HCAL::GetAdjacentCells(const HcalTopology* theHBHETopology, HcalDetId *TrackAlignedCells, HcalDetId ClosestCell, const int Ndepths, int ieta, double deta, double dphi, const int CellsPerDepth){ 

  GetCenterCells(theHBHETopology, TrackAlignedCells, ClosestCell, Ndepths, CellsPerDepth);

  for(int i=1;i<=Ndepths;i++)
  {
     if(!theHBHETopology->validHcal(TrackAlignedCells[CellsPerDepth*(i-1)])){continue;}
     HcalDetId IEta[2];
     HcalDetId Corner;
     if(deta<ROOT::Math::Pi())
     {
        if(deta>0){theHBHETopology->incIEta(TrackAlignedCells[CellsPerDepth*(i-1)],IEta);}
	else{theHBHETopology->decIEta(TrackAlignedCells[CellsPerDepth*(i-1)],IEta);}
	if(dphi<ROOT::Math::Pi())
	{
	   int high=0;
	   if(theHBHETopology->validHcal(IEta[1])){high=1;}
	   if(dphi>0){theHBHETopology->incIPhi(IEta[high],Corner);}
           else{theHBHETopology->decIPhi(IEta[0],Corner);}
	}
     }
     HcalDetId IPhi;
     if(dphi<ROOT::Math::Pi())
     {
       if(dphi>0)
       {
          if(theHBHETopology->incIPhi(TrackAlignedCells[CellsPerDepth*(i-1)],IPhi)){TrackAlignedCells[CellsPerDepth*(i-1)+3]=IPhi;} 
       }
       else
       {
          if(theHBHETopology->decIPhi(TrackAlignedCells[CellsPerDepth*(i-1)],IPhi)){TrackAlignedCells[CellsPerDepth*(i-1)+3]=IPhi;}
       }
     }
     TrackAlignedCells[CellsPerDepth*(i-1)+1]=IEta[0];
     TrackAlignedCells[CellsPerDepth*(i-1)+2]=IEta[1];
     TrackAlignedCells[CellsPerDepth*(i-1)+4]=Corner;
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

void HCAL::HitsPlots(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, GlobalPoint TrackGlobalPoint, GlobalPoint RandGlobalPoint, bool GoodRand, Histograms myHistograms, double TrackPt, GlobalPoint VertexPosition){ 

  double MuonEta = TrackGlobalPoint.eta();
  double MuonPhi = TrackGlobalPoint.phi();
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
  int TrackiEta,TrackiPhi;
  //  HcalDetId MuonClosestCell = (HcalDetId)HEGeom->getClosestCell(MuonGlobalPoint); 
  HcalDetId ClosestCell = (HcalDetId)HEGeom->getClosestCell(TrackGlobalPoint); 
  TrackiEta = ClosestCell.ieta();
  TrackiPhi = ClosestCell.iphi();
  if(fabs(TrackiEta)<19){return;}
  if(TrackiEta<-16&&TrackiPhi>52&&TrackiPhi<64){return;}
  if(fabs(TrackiEta)>28){return;}
  const int Ndepths = 7;
  const int CellsPerDepth = 5;
  HcalDetId AdjacentCells[CellsPerDepth*Ndepths];
  HcalDetId LowThreshAdjacentCells[CellsPerDepth*Ndepths];
  HcalDetId CenterCells[Ndepths];
  //GetConeIDs(theHBHETopology,TrackAlignedCells,ClosestCell,Ndepths,CellsPerDepth);
  //GetCenterCells(theHBHETopology,TrackAlignedCells,ClosestCell,Ndepths,CellsPerDepth); 
  double highetathresh, lowetathresh, phithresh;
  if(fabs(TrackiEta)>20)
  {
     highetathresh = 0.015;
     lowetathresh = 0.03;
     phithresh = 0.055;
  }
  else
  {
     lowetathresh = 0.024;
     highetathresh = 0.01;
     phithresh = 0.008;
  }
  /*lowetathresh = 0.0;
  highetathresh = 0.0;
  phithresh = 0.0;*/
  //GetCornerIDs(theHBHETopology,CornerAlignedCells,ClosestCell,Ndepths);
  
  double Tdphi = TrackGlobalPoint.phi()-caloGeom->getGeometry(ClosestCell)->phiPos();
  if(Tdphi>ROOT::Math::Pi()) Tdphi -= 2*ROOT::Math::Pi();
  if(Tdphi<-ROOT::Math::Pi()) Tdphi += 2*ROOT::Math::Pi();
  double Tdeta = TrackGlobalPoint.eta()-caloGeom->getGeometry(ClosestCell)->etaPos();
  bool etaedge = false;
  if(TrackiEta>0)
  {
     if(Tdeta>highetathresh){etaedge=true;}
     if(Tdeta<-lowetathresh){etaedge=true;}
  }
  else
  {
     if(Tdeta<-highetathresh){etaedge=true;}
     if(Tdeta>lowetathresh){etaedge=true;}
//       if(fabs(Tdeta)>lowetathresh){etaedge=true;}
  }
  bool phiedge = (fabs(Tdphi)>phithresh);
  GetAdjacentCells(theHBHETopology,LowThreshAdjacentCells,ClosestCell,Ndepths,TrackiEta,Tdeta,Tdphi,CellsPerDepth);
  if(!etaedge){Tdeta=10;}
  if(!phiedge){Tdphi=10;}
  GetAdjacentCells(theHBHETopology,AdjacentCells,ClosestCell,Ndepths,TrackiEta,Tdeta,Tdphi,CellsPerDepth);
  GetCenterCells(theHBHETopology,CenterCells,ClosestCell,Ndepths,1);
  int ValidIdCount = 0;
  for(int i=0;i<CellsPerDepth*Ndepths;i++){if(theHBHETopology->validHcal(AdjacentCells[i])){ValidIdCount++;}}

  HcalDetId RandClosestCell = (HcalDetId)HEGeom->getClosestCell(RandGlobalPoint);
  int RandiPhi = RandClosestCell.iphi();
  int RandiEta = RandClosestCell.ieta();
  HcalDetId RandAlignedCells[CellsPerDepth*Ndepths];
  HcalDetId RandCenterCells[Ndepths];
  GetCenterCells(theHBHETopology,RandCenterCells,RandClosestCell,Ndepths,1);
  double Rdphi = RandGlobalPoint.phi()-caloGeom->getGeometry(RandClosestCell)->phiPos();
  if(Rdphi>ROOT::Math::Pi()) Rdphi -= 2*ROOT::Math::Pi();
  if(Rdphi<-ROOT::Math::Pi()) Rdphi += 2*ROOT::Math::Pi();
  double Rdeta = RandGlobalPoint.eta()-caloGeom->getGeometry(RandClosestCell)->etaPos();
  bool retaedge = false;
  if(TrackiEta>0)
  {
     if(Rdeta>highetathresh){etaedge=true;}
     if(Rdeta<-lowetathresh){etaedge=true;} 
  }
  else
  {
     if(Tdeta<-highetathresh){etaedge=true;}
     if(Tdeta>lowetathresh){etaedge=true;} 
  }
  bool rphiedge = (fabs(Rdphi)>phithresh);
  
  if(!retaedge){Rdeta=10;}
  if(!rphiedge){Rdphi=10;}
  GetAdjacentCells(theHBHETopology,RandAlignedCells,RandClosestCell,Ndepths,RandiEta,Rdeta,Rdphi,CellsPerDepth);

  if(!hcalRecHits.isValid())
  {
    std::cout << "Could not find HCAL RecHits" << std::endl;
  }
  else
  {
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    std::deque < std::tuple <int, int, double> >  MuonHits[7];
    double layerenergies[7],rlayerenergies[7],lowthreshadjacent[7];    
    for(int i=0;i<7;i++)
    {
       layerenergies[i]=0;
       rlayerenergies[i]=0;
       lowthreshadjacent[i]=0;
    }

    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
    {  
       HcalDetId id(hbherechit->detid());
       std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
       Global3DPoint hbhe_position = hbhe_cell->getPosition();
       
       HcalDetId *trackmatch = std::find(std::begin(AdjacentCells), std::end(AdjacentCells), id);
       HcalDetId *centermatch = std::find(std::begin(CenterCells), std::end(CenterCells), id);
       HcalDetId *lowthreshmatch = std::find(std::begin(LowThreshAdjacentCells), std::end(LowThreshAdjacentCells), id);
       int HitiEta = id.ieta();
       if(fabs(HitiEta)<16){continue;} 
       if(hbherechit->energy()!=0)
       {
         if(trackmatch!=std::end(AdjacentCells))
         {	 
           if(centermatch==std::end(CenterCells))
    	   {
  	     if(hbherechit->energy()>Hit_Thresholds[id.depth()])
  	     {
  	        myHistograms.m_adjacentfailvertex_z->Fill(VertexPosition.z());
  	        //return;
	     }
  	   }
	   else
	   {
  	      Hits[1] += hbherechit->energy();
              if(id.depth()<8) 
  	      {
  	        //if(hbherechit->energy()>layerenergies[id.depth()-1]){layerenergies[id.depth()-1]=hbherechit->energy();}
  	        layerenergies[id.depth()-1]+=hbherechit->energy();
  	      }
  	      myHistograms.m_Layer_Eta[id.depth()-1]->Fill(hbherechit->energy(),hbhe_position.eta());
  	      myHistograms.m_HitDepth_MuonHCAL->Fill(id.depth());
	   }
         }
	 if(lowthreshmatch!=std::end(LowThreshAdjacentCells)&&centermatch==std::end(CenterCells)&&(hbherechit->energy()>lowthreshadjacent[id.depth()]))
	 {
	    lowthreshadjacent[id.depth()]=hbherechit->energy();
	 }
       }
       
       HcalDetId *randmatch = std::find(std::begin(RandAlignedCells), std::end(RandAlignedCells), id); 
       HcalDetId *randcentermatch = std::find(std::begin(RandCenterCells), std::end(RandCenterCells), id);
       if(hbherechit->energy()!=0&&randmatch!=std::end(RandAlignedCells)&&GoodRand)
       {
         if(hbherechit->energy()>Hit_Thresholds[id.depth()]&&randcentermatch==std::end(RandCenterCells)){GoodRand=false;}
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
    int rhitsoverthresh=0;
    for(int i=0;i<7;i++)
    {
       if(fabs(TrackiEta)>25||i<6)
       {
          if(layerenergies[i]<Hit_Thresholds[i])
	  {
	     myHistograms.m_NThreshCut->Fill(0);
	     if(lowthreshadjacent[i]>Hit_Thresholds[i]){return;}
	  }
       }
    }
    myHistograms.m_NThreshCut->Fill(1);
    for(int i=0;i<7;i++) 
    {
       if(fabs(TrackiEta)>25||i<6)
       {
          myHistograms.m_Layer_Spectra[i]->Fill(layerenergies[i]);
	  if(i<6){myHistograms.m_DepthPairSpectra[i]->Fill(layerenergies[i],layerenergies[i+1]);}
       }
       if(layerenergies[i]!=0)
       {
          if(i<6){Hits[0]++;}
	  if(layerenergies[i]>Hit_Thresholds[i]){hitsoverthresh++;}
          MuonHits[i].push_back(std::make_tuple(TrackiEta,TrackiPhi,layerenergies[i]));
       }
       if(GoodRand)
       {
          if(fabs(RandiEta)>25||i<6)
	  {
	     myHistograms.m_RLayer_Spectra[i]->Fill(rlayerenergies[i]);
	     if(i<6){myHistograms.m_RDepthPairSpectra[i]->Fill(rlayerenergies[i],rlayerenergies[i+1]);}
	     if(rlayerenergies[i]>Hit_Thresholds[i]){rhitsoverthresh++;}
	  }
          if(rlayerenergies[i]!=0)
          {
             MuonHits[i].push_back(std::make_tuple(TrackiEta,RandiPhi,rlayerenergies[i]));
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
    myHistograms.m_TrackPt->Fill(TrackPt);
    myHistograms.m_ConeHits->Fill(Hits[0]);
    myHistograms.m_HitsOverThresh->Fill(hitsoverthresh);
    myHistograms.m_histogram_HCALHits_EtaPhi->Fill(TrackiEta,TrackiPhi);
    if(hitsoverthresh<4)
    {
       myHistograms.m_histogram_BlankHCALHits_EtaPhi->Fill(TrackiEta,TrackiPhi);
       myHistograms.m_histogram_misshitVertex_z->Fill(VertexPosition.z());
       myHistograms.m_histogram_misshitVertex_xy->Fill(VertexPosition.x(),VertexPosition.y());
       myHistograms.m_MissThreshConeEnergy->Fill(Hits[1]);
    }
    myHistograms.m_histogram_differenceVertex_z->Fill(VertexPosition.z());
    myHistograms.m_histogram_differenceVertex_xy->Fill(VertexPosition.x(),VertexPosition.y());
    myHistograms.m_ValidIDs->Fill(ValidIdCount);
    if(Hits[0]<4)
    {
       myHistograms.m_Missing_ValidIDs->Fill(ValidIdCount);
       double bdphi = fabs(caloGeom->getGeometry(ClosestCell)->phiPos()-MuonPhi);
       if(bdphi>ROOT::Math::Pi()) bdphi -= 2*ROOT::Math::Pi();
       myHistograms.m_BlankHitsDR->Fill(sqrt( pow(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,2.0) + pow(bdphi,2.0)));
    }    
    double Tdphi = TrackGlobalPoint.phi()-caloGeom->getGeometry(ClosestCell)->phiPos();
    if(Tdphi>ROOT::Math::Pi()) Tdphi -= 2*ROOT::Math::Pi();
    if(Tdphi<-ROOT::Math::Pi()) Tdphi += 2*ROOT::Math::Pi();
    double TrackMinDr = sqrt( pow(caloGeom->getGeometry(ClosestCell)->etaPos()-TrackGlobalPoint.eta(),2.0)+pow(Tdphi,2.0));
    if(Hits[0]<6){myHistograms.m_MissHitTrackPt->Fill(TrackPt);}
    if(Hits[0]<4&&fabs(TrackiEta)>20)
    {
       myHistograms.m_TrackHCALDR_MissHit->Fill(TrackMinDr);
       if(TrackiEta>0){myHistograms.m_BlankCellDetaDphiPosEta->Fill(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,Tdphi);}
       else{myHistograms.m_BlankCellDetaDphiNegEta->Fill(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,Tdphi);}
    }
    if(Hits[0]<4&&fabs(TrackiEta)<21)
    {
       if(TrackiEta>0){myHistograms.m_BlankCellSmallDetaDphiPosEta->Fill(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,Tdphi);}
       else{myHistograms.m_BlankCellSmallDetaDphiNegEta->Fill(caloGeom->getGeometry(ClosestCell)->etaPos()-MuonEta,Tdphi);}
    }

    if(Hits[0]==6&&layerenergies[6]==0){myHistograms.m_TrackHCALDR_GoodHits->Fill(TrackMinDr);}
    myHistograms.m_ConeEnergy->Fill(Hits[1]);
    if(GoodRand)
    {
       myHistograms.m_RandomConeHits->Fill(Hits[2]);
       myHistograms.m_RandomConeEnergy->Fill(Hits[3]);
       myHistograms.m_RandomHitsOverThresh->Fill(rhitsoverthresh);
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
		if(std::get<1>(*muonhit)==TrackiPhi)
		{
                   if(missinghit) 
                   {	
  	              myHistograms.m_MissingHits->Fill(j+1.5);
  	              myHistograms.m_MissingHitsMap->Fill(std::get<0>(*muonhit),std::get<1>(*muonhit));
  	              if(midE!=0){myHistograms.m_MissingHitsEnergy->Fill(midE);}
                      myHistograms.m_MissingHitsDR->Fill(MuonDR);
  	           }
  	           else 
  	           {
  	              myHistograms.m_MissingHits->Fill(-j-2.5);
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

bool HCAL::FindMuonHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, GlobalPoint TrackGlobalPoint, MCHistograms myHistograms, double standaloneE, double weight, double vtxz, double charge, reco::TransientTrack track){ 
  double MuonEta = TrackGlobalPoint.eta();
  double MuonPhi = TrackGlobalPoint.phi();
  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
  double ReducedConeE=0;
  static double Hits[4];
  Hits[0] = 0;
  Hits[1] = 0;
  Hits[2] = 0;
  Hits[3] = 0;
  m_HitsOverThresh = 0;
  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  
  edm::ESHandle<HcalTopology> htopo;
  iSetup.get<HcalRecNumberingRecord>().get(htopo);
  const HcalTopology* theHBHETopology = htopo.product();
  
  const CaloGeometry* caloGeom = TheCALOGeometry.product();
  const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);
  int TrackiEta,TrackiPhi;
  //  HcalDetId MuonClosestCell = (HcalDetId)HEGeom->getClosestCell(MuonGlobalPoint); 
  HcalDetId ClosestCell = (HcalDetId)HEGeom->getClosestCell(TrackGlobalPoint); 
  TrackiEta = ClosestCell.ieta();
  TrackiPhi = ClosestCell.iphi();
  int BremDepth = ClosestCell.depth();
  myHistograms.m_BremDepth->Fill(BremDepth);
  myHistograms.m_WeightedBremDepth->Fill(BremDepth, weight);
  if(fabs(TrackiEta)<19){return false;}
  if(TrackiEta<-16&&TrackiPhi>52&&TrackiPhi<64){return false;}
  if(fabs(TrackiEta)>28){return false;}
  const int Ndepths = 7;
  const int CellsPerDepth = 5;
  HcalDetId AdjacentCells[CellsPerDepth*Ndepths];
  HcalDetId LowThreshAdjacentCells[CellsPerDepth*Ndepths];
  HcalDetId CenterCells[CellsPerDepth*Ndepths];
  //GetConeIDs(theHBHETopology,TrackAlignedCells,ClosestCell,Ndepths,CellsPerDepth);
  //GetCenterCells(theHBHETopology,TrackAlignedCells,ClosestCell,Ndepths,CellsPerDepth); 
  double highetathresh, lowetathresh, phithresh;
  if(fabs(TrackiEta)>20)
  {
     highetathresh = 0.015;
     lowetathresh = 0.03;
     phithresh = 0.055;
  }
  else
  {
     lowetathresh = 0.024;
     highetathresh = 0.01;
     phithresh = 0.008;
  }
  /*lowetathresh = 0.0;
  highetathresh = 0.0;
  phithresh = 0.0;*/
  //GetCornerIDs(theHBHETopology,CornerAlignedCells,ClosestCell,Ndepths);
  
  double Tdphi = TrackGlobalPoint.phi()-caloGeom->getGeometry(ClosestCell)->phiPos();
  double Tdeta = TrackGlobalPoint.eta()-caloGeom->getGeometry(ClosestCell)->etaPos();
  if(Tdphi>ROOT::Math::Pi()) Tdphi -= 2*ROOT::Math::Pi();
  if(Tdphi<-ROOT::Math::Pi()) Tdphi += 2*ROOT::Math::Pi();
  bool etaedge = false;
  if(TrackiEta>0)
  {
     if(Tdeta>highetathresh){etaedge=true;}
     if(Tdeta<-lowetathresh){etaedge=true;}
  }
  else
  {
     if(Tdeta<-highetathresh){etaedge=true;}
     if(Tdeta>lowetathresh){etaedge=true;}
//       if(fabs(Tdeta)>lowetathresh){etaedge=true;}
  }
  bool phiedge = (fabs(Tdphi)>phithresh);
  GetAdjacentCells(theHBHETopology,LowThreshAdjacentCells,ClosestCell,Ndepths,TrackiEta,Tdeta,Tdphi,CellsPerDepth);
  if(!etaedge){Tdeta=10;}
  if(!phiedge){Tdphi=10;}
  GetAdjacentCells(theHBHETopology,AdjacentCells,ClosestCell,Ndepths,TrackiEta,Tdeta,Tdphi,CellsPerDepth);
  //  GetCenterCells(theHBHETopology,CenterCells,ClosestCell,Ndepths,1);
  CellsFound = GetTransientProjectedCells(HEGeom,CenterCells,track);
  Tdphi = TrackGlobalPoint.phi()-caloGeom->getGeometry(ClosestCell)->phiPos();
  Tdeta = TrackGlobalPoint.eta()-caloGeom->getGeometry(ClosestCell)->etaPos();

  int ValidIdCount = 0;
  for(int i=0;i<CellsPerDepth*Ndepths;i++){if(theHBHETopology->validHcal(AdjacentCells[i])){ValidIdCount++;}}
  if(!hcalRecHits.isValid())
  {
    printf("Could not find HCAL RecHits.\n");
  }
  else
  {
    printf("Successfully found HCAL RecHits.\n");
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    double layerenergies[7],lowthreshadjacent[7];    
    for(int i=0;i<7;i++)
    {
       layerenergies[i]=0;
       lowthreshadjacent[i]=0;
    }

    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
    {  
       printf("Found an HCAL hit, energy is %f.\n",hbherechit->energy());
       HcalDetId id(hbherechit->detid());
       std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
//       Global3DPoint hbhe_position = hbhe_cell->getPosition();
       
       HcalDetId *trackmatch = std::find(std::begin(AdjacentCells), std::end(AdjacentCells), id);
       HcalDetId *centermatch = std::find(std::begin(CenterCells), std::end(CenterCells), id);
       HcalDetId *lowthreshmatch = std::find(std::begin(LowThreshAdjacentCells), std::end(LowThreshAdjacentCells), id);
       int HitiEta = id.ieta();
       if(fabs(HitiEta)<16){continue;} 
       if(hbherechit->energy()!=0)
       {
         if(centermatch!=std::end(CenterCells))
         {
           Hits[0]++;
           Hits[1]+=hbherechit->energy();
           if(hbherechit->energy()>0.8){ReducedConeE+= hbherechit->energy();}
         }
         if(trackmatch!=std::end(AdjacentCells))
         {	 
           if(centermatch!=std::end(CenterCells))
    	   {
              if(id.depth()<8) 
  	      {
  	        //if(hbherechit->energy()>layerenergies[id.depth()-1]){layerenergies[id.depth()-1]=hbherechit->energy();}
  	        layerenergies[id.depth()-1]+=hbherechit->energy();
  	      }
              if(id.depth()==BremDepth){myHistograms.m_BremSpectrum->Fill(hbherechit->energy());}
  	      myHistograms.m_HitDepth_MuonHCAL->Fill(id.depth());
	   }
         }
	 if(lowthreshmatch!=std::end(LowThreshAdjacentCells)&&centermatch==std::end(CenterCells)&&(hbherechit->energy()>lowthreshadjacent[id.depth()]))
	 {
	    lowthreshadjacent[id.depth()]=hbherechit->energy();
	 }
       }
//       hbhe_cell->reset();
    }
    int hitsoverthresh=0;
    if(standaloneE==0)
    { 
       double AdjacentE=0;
       for(int i=0; i<7; i++)
       {
          AdjacentE+=lowthreshadjacent[i];
       }
       myHistograms.m_NStandaloneAdjacentHitEnergies->Fill(AdjacentE);
    }
    m_failAdjacent=false;
    for(int i=0;i<7;i++)
    {
       m_hitEnergies[i]=layerenergies[i];
       if(fabs(TrackiEta)>25||i<6)
       {
          if(layerenergies[i]<Hit_Thresholds[i])
	  {
	     if(lowthreshadjacent[i]>Hit_Thresholds[i])
	     {
	        myHistograms.m_NThreshCut->Fill(0);
                myHistograms.m_AdjacentFailHitEnergy->Fill(lowthreshadjacent[i]);
                if(standaloneE==0){myHistograms.m_NStandaloneAdjacentFailHitEnergy->Fill(lowthreshadjacent[i]);}
	        m_failAdjacent=true;
	     }
	  }
       }
    }
    myHistograms.m_NThreshCut->Fill(1);
    for(int i=0;i<7;i++) 
    {
       //if(layerenergies[i]!=0)
       //{
	  if(layerenergies[i]>Hit_Thresholds[i]){hitsoverthresh++;}
	  //if(BremDepth-1<i){myHistograms.m_Depth_Spectra[i]->Fill(layerenergies[i]);}
          myHistograms.m_Depth_Spectra[i]->Fill(layerenergies[i]);
       //}
    }   
    if(m_failAdjacent)
    {
      if(charge >0)
      {
         if(MuonEta>0)
         {
           if(TrackiEta>20){myHistograms.m_BigEtaAdjFaildEtadPhiMuPlusEtaPlus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAdjFaildEtadPhiMuPlusEtaPlus->Fill(Tdeta,Tdphi);}
         }
         else
         {
           if(TrackiEta<-20){myHistograms.m_BigEtaAdjFaildEtadPhiMuPlusEtaMinus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAdjFaildEtadPhiMuPlusEtaPlus->Fill(Tdeta,Tdphi);}
         }
      }
      else
      {
         if(MuonEta>0)
         {
           if(TrackiEta>20){myHistograms.m_BigEtaAdjFaildEtadPhiMuMinusEtaPlus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAdjFaildEtadPhiMuMinusEtaPlus->Fill(Tdeta,Tdphi);}

         }
         else
         {
           if(TrackiEta<-20){myHistograms.m_BigEtaAdjFaildEtadPhiMuMinusEtaMinus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAdjFaildEtadPhiMuMinusEtaMinus->Fill(Tdeta,Tdphi);}
         }
      }
    }
    if(charge >0)
      {
         if(MuonEta>0)
         {
           if(TrackiEta>20){myHistograms.m_BigEtaAlldEtadPhiMuPlusEtaPlus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAlldEtadPhiMuPlusEtaPlus->Fill(Tdeta,Tdphi);}
         }
         else
         {
           if(TrackiEta<-20){myHistograms.m_BigEtaAlldEtadPhiMuPlusEtaMinus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAlldEtadPhiMuPlusEtaPlus->Fill(Tdeta,Tdphi);}
         }
      }
      else
      {
         if(MuonEta>0)
         {
           if(TrackiEta>20){myHistograms.m_BigEtaAlldEtadPhiMuMinusEtaPlus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAlldEtadPhiMuMinusEtaPlus->Fill(Tdeta,Tdphi);}

         }
         else
         {
           if(TrackiEta<-20){myHistograms.m_BigEtaAlldEtadPhiMuMinusEtaMinus->Fill(Tdeta,Tdphi);}
           else{myHistograms.m_SmallEtaAlldEtadPhiMuMinusEtaMinus->Fill(Tdeta,Tdphi);}
         }
      }
    myHistograms.m_ConeHits->Fill(Hits[0]);
    myHistograms.m_HitsOverThresh->Fill(hitsoverthresh);
    myHistograms.m_HitsOverThreshSplit[BremDepth-1]->Fill(hitsoverthresh);
    myHistograms.m_ConeEnergy->Fill(Hits[1]);
    myHistograms.m_ReducedHEConeEnergy->Fill(ReducedConeE);
    ConeEnergy = Hits[1];
    m_HitsOverThresh = hitsoverthresh;
    if(standaloneE==0){myHistograms.m_NMatchNStandaloneHitsOverThresh->Fill(hitsoverthresh);}
    if(standaloneE!=0)
    {
       myHistograms.m_NonMatchedHEDeposit->Fill(Hits[1]);
       myHistograms.m_NonMatchedHitsOverThresh->Fill(hitsoverthresh);
    } 
    if(standaloneE<80.){myHistograms.m_LowEHCAL->Fill(Hits[1]);}
  }
  return true;
}

int HCAL::GetProjectedCells(const CaloSubdetectorGeometry* HEGeom, HcalDetId *TrackAlignedCells, double vtxz, GlobalPoint direction)
{

   double start, step, end;
   start = 320;
   step = 5;
   end = 530;
   int j=0;
   HcalDetId lastClosestCell;
   for(int i = start; i<end; i+=step)
   {
      double testPointPerp = fabs(i*tan(direction.theta()));
      double testPointX = testPointPerp*cos(direction.phi());
      double testPointY = testPointPerp*sin(direction.phi());
      double testPointZ;
      if(direction.eta()>0){testPointZ = vtxz+i;}
      else
      {
        testPointZ = vtxz-i;	
      }
      GlobalPoint testGlobalPoint = GlobalPoint(testPointX,testPointY,testPointZ);
      HcalDetId testClosestCell = (HcalDetId)HEGeom->getClosestCell(testGlobalPoint);
      if(testClosestCell!=lastClosestCell)
      {
        lastClosestCell=testClosestCell;
        TrackAlignedCells[j] = testClosestCell;
	j = j+1;
      }
   }
   return j;
}

double HCAL::GetIsolation(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, const reco::TransientTrack track)
{
   edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
   iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
   double MatchedEnergy = 0;

   edm::ESHandle<CaloGeometry> TheCALOGeometry;
   iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);

   const CaloGeometry* caloGeom = TheCALOGeometry.product();
   const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);

/*   HcalDetId MatchedCells[20];

   GetTransientProjectedCells(HEGeom,MatchedCells,track);

   if(!hcalRecHits.isValid()) return 1000;
   const HBHERecHitCollection *hbhe = hcalRecHits.product();
   for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
   {
      HcalDetId id(hbherechit->detid());
      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());

      HcalDetId *match = std::find(std::begin(MatchedCells), std::end(MatchedCells), id);
      if(match!=std::end(MatchedCells)) MatchedEnergy+=hbherechit->energy();
   }
   */
   if(!hcalRecHits.isValid()) return 1000;
   const HBHERecHitCollection *hbhe = hcalRecHits.product();
   for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++)
   {
      HcalDetId id(hbherechit->detid());
      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      const GlobalPoint hitPos = hbhe_cell->getPosition();
      TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
      math::XYZVector idPositionRoot(hitPos.x(),hitPos.y(),hitPos.z());
      math::XYZVector trajRoot(traj.position().x(),traj.position().y(),traj.position().z());
      if(ROOT::Math::VectorUtil::DeltaR(idPositionRoot,trajRoot)<0.4)
      {
         MatchedEnergy+=hbherechit->energy();
      }
   }
 
   return MatchedEnergy;
}

int HCAL::GetTransientProjectedCells(const CaloSubdetectorGeometry* HEGeom, HcalDetId *TrackAlignedCells, reco::TransientTrack muTrack)
{

   double start, step, end;
   start = 320;
   step = 5;
   end = 530;
   int j=0;
   HcalDetId lastClosestCell;
   for(int i = start; i<end; i+=step)
   {
      double testPointPerp = fabs(i*tan(muTrack.track().theta()));
      double testPointX = testPointPerp*cos(muTrack.track().phi());
      double testPointY = testPointPerp*sin(muTrack.track().phi());
      double testPointZ; 
      if(muTrack.track().eta()>9){testPointZ = i;}
      else{testPointZ=-i;}

      GlobalPoint testGlobalPoint = GlobalPoint(testPointX,testPointY,testPointZ);
      TrajectoryStateClosestToPoint traj = muTrack.trajectoryStateClosestToPoint(testGlobalPoint);
      HcalDetId testClosestCell = (HcalDetId)HEGeom->getClosestCell(traj.position());
      if(testClosestCell!=lastClosestCell)
      {
        lastClosestCell=testClosestCell;
        TrackAlignedCells[j] = testClosestCell;
	j = j+1;
      }
   }
   return j;
}
