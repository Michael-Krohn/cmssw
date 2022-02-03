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
#include "TNtuple.h"
#include <iostream>

MuPXHistograms::MuPXHistograms()
{
   cutProgress = -0.5;
}

void MuPXHistograms::book(TFileDirectory histFolder, bool binIso){
  m_eventWeight = histFolder.make<TH1F>("eventWeight",";Weight; Events",200, -100, 100);
  m_eventCount = histFolder.make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
  m_SigCat = histFolder.make<TH1F>("Signal Category", "; Signal Category", 2, -0.5, 1.5);
  m_cutProgress = histFolder.make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 30, -.5, 29.5);
  m_MuonTrackMass = histFolder.make<TH1F>("MuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50  , 150  );
//  if(binIso)
//  {
  for(int i=0;i<100;i++)
  {
     std::string name = "EcalMuonsTrackMass"+std::to_string(i);
     m_EcalSplitMuonTrackMass[i] = histFolder.make<TH1F>(name.c_str(), "; MuonTrackMass (GeV); Events", 100,50,150);
     std::string adjname = "adj"+name; 
     m_EcalSplitAdjMuonTrackMass[i] = histFolder.make<TH1F>(adjname.c_str(), "; MuonTrackMass (GeV); Events", 110, 50, 160);  
  }
  for(int i=0;i<100;i++)
  {
     std::string name = "TrackMuonsTrackMass"+std::to_string(i);
     m_TrackSplitMuonTrackMass[i] = histFolder.make<TH1F>(name.c_str(), "; MuonTrackMass (GeV); Events", 100,50,150);
     std::string adjname = "adj"+name;
     m_TrackSplitAdjMuonTrackMass[i] = histFolder.make<TH1F>(adjname.c_str(), "; MuonTrackMass (GeV); Events", 110, 50, 160);  
  }
 
  for(int i=0;i<20;i++)
  {
     for(int j=0;j<20;j++)
     {
        std::string name = "MuonsTrackMass"+std::to_string(i)+"_"+std::to_string(j);
        m_SplitMuonsTrackMass[i][j] = histFolder.make<TH1F>(name.c_str(), "; MuonTrackMass (GeV); Events", 100,50,150);
        std::string tuplename = "tuple"+name;
        m_TupleSplitMuonsTrackMass[i][j] = histFolder.make<TNtuple>(tuplename.c_str(), "X values and weights","x:y");
     }
  }
  m_TagEta = histFolder.make<TH1F>("TaggingMuonEta", "; Tagging Muon #eta; Events", 100, -2.6, 2.6);
  m_TagPt = histFolder.make<TH1F>("TaggingMuonPt","; Tagging Muon pt (GeV); Events", 200, 0, 120);
  m_TagPhi = histFolder.make<TH1F>("TaggingMuonPhi","; Tagging Muon phi; Events", 100, -3.15, 3.15);
  m_TagEtaPhi = histFolder.make<TH2F>("TaggingMuonEtaPhi","; Tagging Muon eta; Tagging Muon phi; Events", 100, -2.6,2.6, 100, -3.15, 3.15);
  m_NPassingTag = histFolder.make<TH1F>("NumberOfMuonsPassingTag","; # of Selected Muons; Events", 20, -0.5, 19.5);
  m_TagProbeVtxChi = histFolder.make<TH1F>("PairedVtxReducedChiSq", "; Paired Vertex Reduced #Chi Squared; Events", 100, 0, 10);
  m_ProbePt = histFolder.make<TH1F>("ProbeTrackPt", "; Pt of Probe Track; Events", 100, 0, 120);
  m_ProbeEta = histFolder.make<TH1F>("ProbeTrackEta", "; #eta of Probe Track; Events", 100, -2.6, 2.6);
  m_ProbePhi = histFolder.make<TH1F>("ProbeTrackPhi", "; #phi of ProbeTrack; Events", 100, -3.15, 3.15);
  m_ProbePhiDXY = histFolder.make<TH2F>("ProbePhiDXY", "; #phi of Probe Track; DXY of probe Track; Events", 100,-3.15,3.15,100,-2,2);
  m_PreDXYPhi = histFolder.make<TH1F>("PreDXYPhi", "; #phi of ProbeTrack; Events", 100, -3.15, 3.15);
  m_PostDXYPhi = histFolder.make<TH1F>("PostDXYPhi", "; #phi of ProbeTrack; Events", 100, -3.15, 3.15);
  m_PreDXY = histFolder.make<TH1F>("PreDXY","; dxy of Track;Events",200,-2,2);
  m_PostDXY = histFolder.make<TH1F>("PostDXY","; dxy of Track;Events",200,-2,2);
  m_ProbeDCA = histFolder.make<TH1F>("ProbeDCA","; Probe Distance of Closest Approach;Events",200,0,0.05);
  m_ProbeLongDCA = histFolder.make<TH1F>("ProbeLongDCA","; Probe Longitudinal Distance of Closest Approach;Events",200,0,0.3);
  m_ProbeDCAPhi = histFolder.make<TH2F>("ProbeDCAPhi","; Probe Distance of Closest Approach;Probe Phi;Events",200,0,0.2,100,-3.15,3.15);

  m_ProbeEtaPhi = histFolder.make<TH2F>("ProbeTrackEtaPhi","; Probe Track #eta; Probe Track #phi; Events", 100, -2.6,2.6, 100, -3.15, 3.15);
  m_ProbeTrackIso = histFolder.make<TH1F>("ProbeTrackIsolation","; Track Based Isolation; Events",200,0,5);
  m_ProbeHcalIso = histFolder.make<TH1F>("ProbeHcalIsolation","; Hcal Energy in matched hits (GeV); Events", 200, 0, 100);
  m_ProbeEcalIso = histFolder.make<TH1F>("ProbeEcalIsolation","; Ecal Energy within cone (GeV); Events", 200, 0, 100);
  m_ProbeCombinedIso = histFolder.make<TH2F>("ProbeCombinedIsolation","; Track Based Isolation; Hcal energy in matched hits (GeV); Events", 200, 0, 5, 200, 0, 8);
  m_NPassingProbe = histFolder.make<TH1F>("NTracksPassingProbe","; Tracks Passing Probe Selection; Events", 20, -0.5, 19.5);
  m_AdjustedMuTrackMass = histFolder.make<TH1F>("AdjustedMuonTrackMass","; Adjusted Invariant Mass of Tag and Probe (GeV); Events", 110, 50, 160); 
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
  m_ProbeCaloJetDr = histFolder.make<TH1F>("ProbeCaloJetDr",";#Delta R from Probe to Nearest Calo Jet; Events", 100,0,6);
  m_ProbeJetE = histFolder.make<TH1F>("ProbeCaloJetE",";Calo Jet Pt (GeV); Events", 100, 0, 200);
  m_TagCaloJetE = histFolder.make<TH1F>("TagCaloGetE","; Tag Muon Calo Jet E (GeV); Events", 100, 0, 200);
  m_CaloJetEcalE = histFolder.make<TH2F>("CaloEmEcalIso",";Calo Jet Ecal E; Ecal Iso; Events", 100, 0, 200, 100,0, 200);
  m_CaloJetEcalEnergy = histFolder.make<TH1F>("CaloEm",";Calo Jet Ecal E; Events",100,0,200);
  m_CaloJetHcalE = histFolder.make<TH2F>("CaloHadHcalIso",";Calo Jet Hcal E; Hcal Iso; Events", 100, 0, 200, 100,0, 200);
  m_CaloJetEcalHcal = histFolder.make<TH2F>("CaloEMHad",";Calo Jet Ecal E; Calo Jet Hcal E; Events", 100, 0, 200, 100,0, 200);
  m_CaloJetHcalEnergy = histFolder.make<TH1F>("CaloHad",";Calo Jet Hcal E; Events", 100,0,200);
  ///vtx stuff
  m_ProbeVtx = histFolder.make<TH1F>("ProbeVtxIndex",";Probe Vertex Index; Events",10,-1.5,8.5);
  m_TagVtx = histFolder.make<TH1F>("TagVtxIndex",";Tag Vertex Index; Events",10,-1.5,8.5);
  m_NVertices = histFolder.make<TH1F>("NPV",";Number of PV in Event; Events",100,-0.5,99.5);
  m_TagProbeVtx = histFolder.make<TH2F>("TagProbeVtxIndex",";Tag Vtx Index;Probe Vtx Index; Events",10,-1.5,8.5,10,-1.5,8.5);
  m_TagProbeVtxd0 = histFolder.make<TH1F>("TagProbeVtxD0",";Vertex D0; Events",300,0,30);
  m_TagProbeVtxdz = histFolder.make<TH1F>("TagProbeVtxdz",";Vertex dz; Events",300,0,30);
  //Pileup
  m_PileupWeights = histFolder.make<TH1F>("PileupWeights",";N Pileup Interactions; Weight", 100,0,4);
  m_PUmean = histFolder.make<TH1F>("NPUMean",";MC N Pileup Mean; Events", 100,0,100);
  m_GenProbeDr = histFolder.make<TH1F>("GenProbeDr",";#Delta R Between Gen Probe and Selected Track; Events",100,0,5);
  m_GenProbeEta = histFolder.make<TH1F>("GenProbeEta",";#eta of Gen Probe Muon; Events",200,-5.2,5.2);
  m_GenTagDr = histFolder.make<TH1F>("GenTagDr",";#Delta R Between Gen Tag and Selected Muon; Events", 100,0,5);
  m_InTrackerGenProbeDr = histFolder.make<TH1F>("InTrackerGenProbeDr",";#Delta R Between Gen Probe and Selected Track; Events",100,0,5);
  m_GenProbePt = histFolder.make<TH1F>("GenProbePt",";Pt of Gen Probes That are not Selected; Events", 100,0,200);
  m_PassLooseMatch = histFolder.make<TH1F>("PassLooseMatch",";Pass Loose Matching; Events",2,0,2);
  m_nonMuPdgId = histFolder.make<TH1F>("Non Muon Pdg Id",";Pdg Id of All Particles within #Delta R of 0.2 to Probe; Events", 600,-300,300);
  m_highPtNonMuPdgId = histFolder.make<TH1F>("High Pt Non Muon Pdg Id",";Pdg Id of all Particles >20 GeV Pt near probe; Events", 600,-300,300);
  m_GenPdgId = histFolder.make<TH1F>("Gen Mother PDG Id",";PDG Id of Gen Matched Muon; Events",200,-100,100);
  m_nLooseMuons = histFolder.make<TH1F>("N loose muons",";Number of Gen Muons; Events",10,0,10);
  m_TagPdgId = histFolder.make<TH1F>("Tag Mother PDG Id",";PDG Id of Tag Mother; Events",200,-100,100);
  m_MotherMatch = histFolder.make<TH1F>("MatchTagProbeMother",";Match PDG Ids; Events",2,0,2);
  TFileDirectory looseMatch = histFolder.mkdir("looseMatch");
  for(int i=0;i<10;i++)
  {
     for(int j=0;j<10;j++)
     {
        std::string name = "LooseMatch"+std::to_string(i)+"_"+std::to_string(j);
        m_LooseMatch[i][j] = looseMatch.make<TH1F>(name.c_str(), "; Pass Loose Matching (GeV); Events", 2,0,2);
     }
  }
  m_staMinDr = histFolder.make<TH1F>("staMinDr",";Min #Delta R from Probe to Standalone Muon; Events", 100,0,5);
  m_standaloneE = histFolder.make<TH1F>("standaloneE",";Energy of nearest standalone muon to probe (GeV);Events",100,0,100);
  m_standaloneDE = histFolder.make<TH1F>("StandaloneMuondEoverE","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1);
  m_nEle = histFolder.make<TH1F>("nElectrons",";Number of Veto electrons > 10 GeV; Events",20,0,20);
  m_nExtraEle = histFolder.make<TH1F>("nExtraElectrons",";Number of Extra Electrons in Event; Events",20,0,20);
  m_nExtraMuon = histFolder.make<TH1F>("nExtraMuons",";Number of Extra Muons in Event; Events",20,0,20);
  m_ExtraMuonProbeDr = histFolder.make<TH1F>("ExtraMuonProbeDr",";#Delta R from Probe to Extra Muon; Events",50,0,5);
  m_ExtraMuonTagDr = histFolder.make<TH1F>("ExtraMuonTagDr",";#Delta R from Tag to Extra Muon; Events",50,0,5);
  m_ExtraMuonPt = histFolder.make<TH1F>("ExtraMuonPt",";Pt of Extra Muon; Events",50,0,100);
}

void MuPXHistograms::FillHists(MuPXEventInfo info, bool binIso, bool makeTuples)
{
     m_PileupWeights->Fill(info.pileupWeight);
     m_PUmean->Fill(info.nPUmean);
     m_eventWeight->Fill(info.eventWeight);
     m_eventCount->Fill(0.5,info.eventWeight);
     m_NPassingProbe->Fill(info.nPassingProbe,info.eventWeight);
     m_NPassingTag->Fill(info.nPassingTag,info.eventWeight);
     //for(int i=0; i<info.cutProgress;i++){IncCutFlow();} 
     if(!info.paired) return;  
     m_ExtraMuonProbeDr->Fill(info.extraMuProbeDr,info.eventWeight);
     m_ExtraMuonTagDr->Fill(info.extraMuTagDr,info.eventWeight);
     m_ExtraMuonPt->Fill(info.extraMuPt,info.eventWeight);
     m_ProbeDCA->Fill(info.probeDCA,info.eventWeight);
     m_ProbeLongDCA->Fill(info.probeLongDCA,info.eventWeight);
     m_ProbeDCAPhi->Fill(info.probeDCA,info.probeTrack->phi(),info.eventWeight);
     m_SigCat->Fill(info.sigCat,info.eventWeight);
     m_staMinDr->Fill(info.staMinDr,info.eventWeight);
     m_standaloneE->Fill(info.standaloneE,info.eventWeight);
     m_standaloneDE->Fill(info.standaloneDEoverE,info.eventWeight);
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
     m_AdjustedMuTrackMass->Fill(info.adjMuonTrackMass,info.eventWeight);
     m_nEle->Fill(info.nEle,info.eventWeight);
     m_nExtraEle->Fill(info.nExtraEle,info.eventWeight);
     m_nExtraMuon->Fill(info.nExtraMuon,info.eventWeight);
     int ecalbin = (int)(info.probeEcalIso/5.);
     if(ecalbin>19){ecalbin=19;}
     int tkbin = (int)(info.probeTrackIso/info.probeTrack->pt()*40.);
     if(tkbin>19){tkbin=19;}

     if(binIso)
     {
        if(info.probeEcalIso<100)
        {
          m_EcalSplitAdjMuonTrackMass[(int)info.probeEcalIso]->Fill(info.adjMuonTrackMass,info.eventWeight);
        }
        else
        {
          m_EcalSplitAdjMuonTrackMass[99]->Fill(info.adjMuonTrackMass,info.eventWeight);
        }
        if(info.probeTrackIso/info.probeTrack->pt()<0.5)
        {
          m_TrackSplitAdjMuonTrackMass[(int)(info.probeTrackIso/info.probeTrack->pt()/0.5*100)]->Fill(info.adjMuonTrackMass,info.eventWeight);
        }
        else
        {
          m_TrackSplitAdjMuonTrackMass[99]->Fill(info.adjMuonTrackMass,info.eventWeight);
        }
        m_SplitMuonsTrackMass[ecalbin][tkbin]->Fill(info.adjMuonTrackMass,info.eventWeight); 
     }
     if(makeTuples){m_TupleSplitMuonsTrackMass[ecalbin][tkbin]->Fill(info.adjMuonTrackMass,info.eventWeight);}
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
     //Calo Jets
     m_ProbeJetE->Fill(info.probeCaloJetE,info.eventWeight); 
     m_TagCaloJetE->Fill(info.tagCaloJetE, info.eventWeight);
     m_ProbeCaloJetDr->Fill(info.minProbeCaloJetDr, info.eventWeight);
     m_CaloJetEcalE->Fill(info.probeCaloJetEcal,info.probeEcalIso, info.eventWeight);
     m_CaloJetEcalEnergy->Fill(info.probeCaloJetEcal,info.eventWeight);
     m_CaloJetHcalE->Fill(info.probeCaloJetHcal,info.probeHcalIso, info.eventWeight);
     m_CaloJetHcalEnergy->Fill(info.probeCaloJetHcal,info.eventWeight);
     m_CaloJetEcalHcal->Fill(info.probeCaloJetEcal,info.probeCaloJetHcal,info.eventWeight);
     //Vtx Stuff
     m_ProbeVtx->Fill(info.probeVtx,info.eventWeight);
     m_TagVtx->Fill(info.tagVtx,info.eventWeight);
     m_NVertices->Fill(info.nVtx,info.eventWeight);
     m_TagProbeVtx->Fill(info.tagVtx,info.probeVtx,info.eventWeight);
     m_TagProbeVtxd0->Fill(info.vtxd0,info.eventWeight);
     m_TagProbeVtxdz->Fill(info.vtxdz,info.eventWeight);
     //Gen Particle Information
     
     if(info.farMuDr>0)
     {
        m_GenPdgId->Fill(info.genPdgId,info.eventWeight);
        m_TagPdgId->Fill(info.tagPdgId,info.eventWeight);
        for(uint8_t x=0; x<info.nearbyNonMuPdgIds.size();x++){m_nonMuPdgId->Fill(info.nearbyNonMuPdgIds.at(x),info.eventWeight);}
        for(uint8_t x=0; x<info.highPtNonMuPdgIds.size();x++){m_highPtNonMuPdgId->Fill(info.highPtNonMuPdgIds.at(x),info.eventWeight);} 
        if(info.genPdgId==info.tagPdgId){m_MotherMatch->Fill(1.5,info.eventWeight);}
        else{m_MotherMatch->Fill(0.5,info.eventWeight);}
        m_nLooseMuons->Fill(info.nLooseMuons,info.eventWeight);
        if(info.passLooseMatch)
        {
           m_PassLooseMatch->Fill(1.5,info.eventWeight);
           if(binIso){m_LooseMatch[ecalbin][tkbin]->Fill(1.5,info.eventWeight);}
        }
        else
        {
           m_PassLooseMatch->Fill(0.5,info.eventWeight);
           if(binIso){m_LooseMatch[ecalbin][tkbin]->Fill(0.5,info.eventWeight);}
        }
        m_GenProbeDr->Fill(info.farMuDr,info.eventWeight);
        m_GenProbeEta->Fill(info.farMuEta,info.eventWeight);
        m_GenTagDr->Fill(info.genTagMuDr,info.eventWeight);
        if(fabs(info.farMuEta)<2.5)
        {
           m_InTrackerGenProbeDr->Fill(info.farMuDr,info.eventWeight);
           m_GenProbePt->Fill(info.farMuE,info.eventWeight);
        }
     }
}

void MuPXHistograms::IncCutFlow(double weight)
{
   cutProgress++;
   m_cutProgress->Fill(cutProgress, weight);
   return;
}

void MuPXHistograms::ResetCutFlow()
{
   cutProgress = -0.5;
}

