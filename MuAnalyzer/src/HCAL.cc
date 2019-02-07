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

#include <iostream>

HCAL::HCAL(){}

void HCAL::CheckHCAL(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label){

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

      const CaloCellGeometry *hbhe_cell = caloGeom->getGeometry(hbherechit->id());
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

       const CaloCellGeometry *hbhe_cell = caloGeom->getGeometry(hbherechit->id());
       Global3DPoint hbhe_position = hbhe_cell->getPosition();

       double dPhi = fabs(MuonPhi - hbhe_position.phi());
       if(dPhi > ROOT::Math::Pi()) dPhi -= 2*ROOT::Math::Pi();

       if(minHCALdR > sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)))){
	 minHCALdR = sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)));
	 MuonHitEnergy = hbherechit->energy();
       }

    }
  }

  return minHCALdR;
}
