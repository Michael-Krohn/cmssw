#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DarkPhoton/MuAnalyzer/interface/SimHistograms.h"
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
#include "TTree.h"
#include <iostream>

SimHistograms::SimHistograms() {
}

void SimHistograms::book(TFileDirectory histFolder) {
  m_initialMuE = histFolder.make<TH1F>("MuInitialE", ";Energy (GeV);Events", 200, 0, 200);
  m_finalMuE = histFolder.make<TH1F>("MuFinalE", ";Energy (GeV);Events", 200, 0, 200);
  m_initialMuEt = histFolder.make<TH1F>("MuInitialEt", "; Et (GeV); Events", 200,0,200);
}

void SimHistograms::FillHists(SigEventInfo info) {
  m_initialMuE->Fill(info.dbremInitialE,info.eventWeight);
  m_finalMuE->Fill(info.dbremFinalE,info.eventWeight);
  m_initialMuEt->Fill(info.dbremInitialEt,info.eventWeight);
}
