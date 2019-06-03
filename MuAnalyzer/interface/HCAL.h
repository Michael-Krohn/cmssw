#ifndef HCAL_H
#define HCAL_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "TH1F.h"
#include "TH2F.h"

class HCAL{
  public:
    HCAL();
    void CheckHCAL(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label);
    double MuonMindR(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, double MuonEta, double MuonPhi, GlobalPoint MuonGlobalPoint);
    void HitsPlots(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, double MuonEta, double MuonPhi, GlobalPoint MuonGlobalPoint,double ConeSize, Histograms myHistograms, double CSCMuonMinDr);
    double MuonHitEnergy;
    int    MuonHitDepth;
    double MuonMinDr;
    HcalDetId minHCALDetId;

  private:
    const double Hit_Thresholds[7] = {200,0.2,0.2,0.2,0.2,0.2,0.25};
    void GetConeIDs(const HcalTopology* theHBHETopology, HcalDetId *MuonAlignedCells, HcalDetId ClosestCell, const int Ndepths);
};

#endif
