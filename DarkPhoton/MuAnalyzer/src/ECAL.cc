#include "DarkPhoton/MuAnalyzer/interface/ECAL.h"
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

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "DarkPhoton/MuAnalyzer/interface/MCHistograms.h"
#include "Math/VectorUtil.h"
#include "TH1F.h"

#include <algorithm>
#include <iostream>

ECAL::ECAL(){}

double ECAL::GetIsolation(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label, edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label, const reco::TransientTrack track)
{
   edm::Handle<EcalRecHitCollection> rechitsEE;
   iEvent.getByToken(reducedEndcapRecHitCollection_Label, rechitsEE);

   edm::Handle<EcalRecHitCollection> rechitsEB;
   iEvent.getByToken(reducedBarrelRecHitCollection_Label, rechitsEB);
 
   edm::ESHandle<CaloGeometry> TheCALOGeometry;
   iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);

   const CaloGeometry* caloGeom = TheCALOGeometry.product();

   double eDR = 0;

   for(EcalRecHitCollection::const_iterator hit = rechitsEE->begin(); hit!= rechitsEE->end(); hit++)
   {
      const DetId id = (*hit).detid();
      const GlobalPoint hitPos = caloGeom->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
      TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
      math::XYZVector idPositionRoot(hitPos.x(),hitPos.y(),hitPos.z());
      math::XYZVector trajRoot(traj.position().x(),traj.position().y(),traj.position().z());
      if(ROOT::Math::VectorUtil::DeltaR(idPositionRoot,trajRoot)<0.4&&(*hit).energy()>0.3)
      {
         eDR+= (*hit).energy();
      }            
   } 
   return eDR;
}
