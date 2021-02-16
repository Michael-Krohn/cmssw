#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DarkPhoton/MuAnalyzer/interface/MuPXHistograms.h"
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

MuPXHistograms::MuPXHistograms()
{
   cutProgress = -0.5;
}

void MuPXHistograms::book(TFileDirectory histFolder){
  m_eventWeight = histFolder.make<TH1F>("eventWeight",";Weight; Events",200, -100, 100);
  m_eventCount = histFolder.make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
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
  //Multiple Paired Tracks Histograms
  m_SmallestCone = histFolder.make<TH1F>("SmallestConeSize","; Smallest #Delta R Cone that contains all paired tracks; Events", 50, 0, 4);
  m_AverageDr = histFolder.make<TH1F>("AverageDr","; Average #Delta R between all paired tracks;", 100, 0, 5);
  //Tagging muon information
  m_TagProbeDr = histFolder.make<TH1F>("TagProbeDr","; #Delta R between Tag and Selected Probe;", 100, 0, 5);
  m_TagTrackIso = histFolder.make<TH1F>("TagTrackIsolation","; Track Based Isolation; Events",200,0,5);
  m_TagHcalIso = histFolder.make<TH1F>("TagHcalIsolation","; Hcal Energy in matched hits (GeV); Events", 200, 0, 100);
  m_TagEcalIso = histFolder.make<TH1F>("TagEcalIsolation","; Ecal Energy within cone (GeV); Events", 200, 0, 100);
  //Jet Histograms
  m_NJets = histFolder.make<TH1F>("Njets","; Number of Jets in Event; Events", 200,0,200);
  m_CaloSumJetE = histFolder.make<TH2F>("CaloSumJetE",";ECAL+HCAL Energy Sum (GeV); Energy of Nearest Jet (GeV); Events", 100,0,200,100,0,200);
  m_HCALIsoJetDr = histFolder.make<TH2F>("HCALIsoJetDr",";Probe HCAL Energy within cone (GeV);#Delta R from probe to nearest PF Jet; Events", 200,0,100,100,0,5);
  m_ECALIsoJetDr = histFolder.make<TH2F>("ECALIsoJetDr",";Probe ECAL Energy within cone (GeV);#Delta R from probe to nearest PF Jet; Events", 200,0,100,100,0,5);
  m_ProbeJetDr = histFolder.make<TH1F>("ProbeJetDr",";#Delta R between selected probe and nearest jet; Events", 120,-1,5);
  m_JetPt = histFolder.make<TH1F>("JetPt",";Jet Pt (GeV); Events", 100,0,150);  
}

void MuPXHistograms::FillHists(MuPXEventInfo info)
{
     m_eventWeight->Fill(info.eventWeight);
     m_eventCount->Fill(0.5,info.eventWeight);
     m_NPassingProbe->Fill(info.nPassingProbe,info.eventWeight);
     m_NPassingTag->Fill(info.nPassingTag,info.eventWeight);
     for(int i=0; i<info.cutProgress;i++){IncCutFlow();} 
     if(!info.paired) return;
     m_ProbeTrackIso->Fill(info.probeTrackIso/info.probeTrack->pt(),info.eventWeight);
     m_ProbeEcalIso->Fill(info.probeEcalIso,info.eventWeight);
     m_ProbeHcalIso->Fill(info.probeHcalIso,info.eventWeight);
     m_TagEta->Fill(info.tagMuon->eta(),info.eventWeight);
     m_TagPhi->Fill(info.tagMuon->phi(),info.eventWeight);
     m_TagEtaPhi->Fill(info.tagMuon->eta(),info.tagMuon->phi(),info.eventWeight);
     m_TagProbeVtxChi->Fill(info.tagProbeVtxChi,info.eventWeight);
     m_TagPt->Fill(info.tagMuon->pt(),info.eventWeight);
     m_ProbePt->Fill(info.probeTrack->pt(),info.eventWeight);
     m_ProbeEta->Fill(info.probeTrack->eta(),info.eventWeight);
     m_ProbePhi->Fill(info.probeTrack->phi(),info.eventWeight);
     m_ProbeEtaPhi->Fill(info.probeTrack->eta(),info.probeTrack->phi(),info.eventWeight);
     m_MuonTrackMass->Fill(info.muonTrackMass,info.eventWeight);
     if(info.averageDr>0)
     {
        m_SmallestCone->Fill(info.smallestCone,info.eventWeight);
        m_AverageDr->Fill(info.averageDr,info.eventWeight);
     }
     m_TagProbeDr->Fill(deltaR(info.probeTrack->eta(),info.probeTrack->phi(),info.tagMuon->eta(),info.tagMuon->phi()),info.eventWeight);
     m_TagTrackIso->Fill(info.tagTrackIso/info.tagMuon->pt(),info.eventWeight);
     m_TagEcalIso->Fill(info.tagEcalIso,info.eventWeight);
     m_TagHcalIso->Fill(info.tagHcalIso,info.eventWeight);
     if(info.nearestJetDr<0.4){m_CaloSumJetE->Fill(info.probeHcalIso+info.probeEcalIso,info.nearestJetE,info.eventWeight);}
     m_NJets->Fill(info.nJets,info.eventWeight);
     m_HCALIsoJetDr->Fill(info.probeHcalIso,info.nearestJetDr,info.eventWeight);
     m_ECALIsoJetDr->Fill(info.probeEcalIso,info.nearestJetDr,info.eventWeight);
     m_ProbeJetDr->Fill(info.nearestJetDr,info.eventWeight);
}

void MuPXHistograms::IncCutFlow()
{
   cutProgress++;
   m_cutProgress->Fill(cutProgress);
   return;
}

void MuPXHistograms::ResetCutFlow()
{
   cutProgress = -0.5;
}

