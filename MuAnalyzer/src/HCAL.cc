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

  std::cout << "Inside CheckHCAL" << std::endl;

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
      std::cout << "hbherechit->id(): " << hbherechit->id() << std::endl;
//      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
//      const CaloCellGeometry *hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      std::cout << "check 1" << std::endl;
//      Global3DPoint hbhe_position = caloGeom->getGeometry(hbherechit->id())->getPosition();
      std::cout << "check 2" << std::endl;      
//      Global3DPoint hbhe_position = hbhe_cell->getPosition();

      std::cout << "hbherechit->energy(): " << hbherechit->energy() << std::endl;
       std::cout << "hbherechit->time(): " << hbherechit->time() << std::endl;

       std::cout << "id.depth(): " << id.depth() << std::endl;
       std::cout << "id.subdet(): " << id.subdet() <<  std::endl;
       std::cout << "hbhe_position.eta(): " << caloGeom->getGeometry(hbherechit->id())->getPosition().eta() << std::endl;
       std::cout << "hbhe_position.phi(): " << caloGeom->getGeometry(hbherechit->id())->getPosition().phi() << std::endl;

    }
  }
}


double HCAL::MuonMindR(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, double MuonEta, double MuonPhi){

  double minHCALdR = 1000;
  std::cout << "inside MuonMindR" << std::endl;

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

  std::cout << "check 3" << std::endl;

  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);
  const CaloGeometry* caloGeom = TheCALOGeometry.product();
  std::cout << "check 4" << std::endl;

  if(!hcalRecHits.isValid()){
    std::cout << "Could not find HCAL RecHits" << std::endl;
  }else{
    std::cout << "check 5" << std::endl;
    const HBHERecHitCollection *hbhe = hcalRecHits.product();
    for(HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++){
       HcalDetId id(hbherechit->detid());
       std::cout << "check 6" << std::endl;

       std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
//       const CaloCellGeometry *hbhe_cell = caloGeom->getGeometry(hbherechit->id());
//       Global3DPoint hbhe_position = caloGeom->getGeometry(hbherechit->id())->getPosition();
       Global3DPoint hbhe_position = hbhe_cell->getPosition();
       std::cout << "check 7" << std::endl;
       double dPhi = fabs(MuonPhi - hbhe_position.phi());
       if(dPhi > ROOT::Math::Pi()) dPhi -= ROOT::Math::Pi();

       if(minHCALdR > sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)))){
	 minHCALdR = sqrt(( pow((MuonEta - hbhe_position.eta()),2.0) + pow(dPhi, 2.0)));
	 std::cout << "check 8" << std::endl;
	 MuonHitEnergy = hbherechit->energy();
       }
//       hbhe_cell->reset();
    }
  }
  std::cout <<" check 9" << std::endl;
  return minHCALdR;
}
