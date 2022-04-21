#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DarkPhoton/MuFilter/interface/eventHistos.h"
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
#include <iostream>

eventHistos::eventHistos()
{
   cutProgress = -0.5;
}

void eventHistos::book(TFileDirectory histFolder){

  m_eventCount = histFolder.make<TH1F>("eventCount", "; ;Events", 1, 0, 1); 
  m_eventWeight = histFolder.make<TH1F>("eventWeight","; Weight; Events", 1000, 0, 1000);
  m_cutProgress = histFolder.make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 10, -.5, 9.5);
  m_MuonTrackMass = histFolder.make<TH1F>("MuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50  , 150  );
  m_TagEta = histFolder.make<TH1F>("TaggingMuonEta", "; Tagging Muon #eta; Events", 100, -2.6, 2.6);
  m_TagPt = histFolder.make<TH1F>("TaggingMuonPt","; Tagging Muon pt (GeV); Events", 200, 0, 120);
  m_TagPhi = histFolder.make<TH1F>("TaggingMuonPhi","; Tagging Muon phi; Events", 100, -3.15, 3.15);
  m_TagEtaPhi = histFolder.make<TH2F>("TaggingMuonEtaPhi","; Tagging Muon eta; Tagging Muon phi; Events", 100, -2.6,2.6, 100, -3.15, 3.15);
  m_NPassingTag = histFolder.make<TH1F>("NumberOfMuonsPassingTag","; # of Selected Muons; Events", 20, -0.5, 19.5);
  m_TagProbeVtxChi = histFolder.make<TH1F>("PairedVtxReducedChiSq", "; Paired Vertex Reduced #Chi Squared; Events", 100, 0, 10);
  m_ProbePt = histFolder.make<TH1F>("ProbeTrackPt", "; Pt of Probe Track; Events", 100, 0, 120);
  m_ProbeEta = histFolder.make<TH1F>("ProbeTrackEta", "; #eta of Probe Track; Events", 100, -2.6, 2.6);
  m_ProbePhi = histFolder.make<TH1F>("ProbeTrackPhi", "; #phi of ProbeTrack; Events", 100, -3.15, 3.15);
  m_ProbeEtaPhi = histFolder.make<TH2F>("ProbeTrackEtaPhi","; Probe Track #eta; Probe Track #phi; Events", 100, -2.6,2.6, 100, -3.15, 3.15);
  m_ProbeTrackIso = histFolder.make<TH1F>("ProbeTrackIsolation","; Track Based Isolation; Events",200,0,5);
  m_ProbeHcalIso = histFolder.make<TH1F>("ProbeHcalIsolation","; Hcal Energy in matched hits (GeV); Events", 200, 0, 100);
  m_ProbeEcalIso = histFolder.make<TH1F>("ProbeEcalIsolation","; Ecal Energy within cone (GeV); Events", 200, 0, 100);
  m_ProbeCombinedIso = histFolder.make<TH2F>("ProbeCombinedIsolation","; Track Based Isolation; Hcal energy in matched hits (GeV); Events", 200, 0, 5, 200, 0, 8);
  m_NPassingProbe = histFolder.make<TH1F>("NTracksPassingProbe","; Tracks Passing Probe Selection; Events", 20, -0.5, 19.5);
}

void eventHistos::IncCutFlow()
{
   cutProgress++;
   m_cutProgress->Fill(cutProgress);
   return;
}

void eventHistos::ResetCutFlow()
{
   cutProgress = -0.5;
}
