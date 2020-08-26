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
#include "DarkPhoton/MuAnalyzer/interface/MCHistograms.h"
#include "TH1F.h"
#include "TH2F.h"

class HCAL{
  public:
    HCAL();
    void CheckHCAL(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label);
    double MuonMindR(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, GlobalPoint MuonGlobalPoint);
    void HitsPlots(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, GlobalPoint TrackGlobalPoint, GlobalPoint RandGlobalPoint, bool GoodRand, Histograms myHistograms, double TrackPt, GlobalPoint VertexPosition);
    bool FindMuonHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> HBHERecHit_Label, GlobalPoint TrackGlobalPoint, MCHistograms myHistograms, double standaloneE, double weight, double vtxz);
    double MuonHitEnergy;
    int    MuonHitDepth;
    double MuonMinDr;
    HcalDetId minHCALDetId;
    double ConeEnergy;
    int m_HitsOverThresh;
    double m_hitEnergies[7];
    bool m_failAdjacent;
    int CellsFound;

  private:
    const double Hit_Thresholds[7] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1};
    void GetConeIDs(const HcalTopology* theHBHETopology, HcalDetId *MuonAlignedCells, HcalDetId ClosestCell, const int Ndepths, const int CellsPerDepth);
    void GetCornerIDs(const HcalTopology* theHBHETopology, HcalDetId *MuonAlignedCells, HcalDetId ClosestCell, const int Ndepths);
    void GetCenterCells(const HcalTopology* theHBHETopology, HcalDetId *MuonAlignedCells, HcalDetId ClosestCell, const int Ndepths, const int CellsPerDepth);
    void GetAdjacentCells(const HcalTopology* theHBHETopology, HcalDetId *MuonAlignedCells, HcalDetId ClosestCell, const int Ndepths, int ieta, double deta, double dphi, const int CellsPerDepth);
    int GetProjectedCells(const CaloSubdetectorGeometry* HEGeom, HcalDetId *TrackAlignedCells, double vtxz, GlobalPoint direction);

};

#endif
