#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DarkPhoton/MuAnalyzer/interface/MCHistograms.h"
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

MCHistograms::MCHistograms()
{
   cutProgress = -0.5;
}

void MCHistograms::book(TFileDirectory histFolder){
  m_genEventCount = histFolder.make<TH1F>("genEventCount", "; ;Events", 1, 0, 1);
  m_genEventCount_EnBin = histFolder.make<TH1F>("genEventCount_EnBin", "; ;Events", 1, 0, 1);
  m_eventCount = histFolder.make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
  m_eventCount_PUup = histFolder.make<TH1F>("eventCount_PUup", "; ;Events", 1, 0, 1);
  m_eventCount_PUdown = histFolder.make<TH1F>("eventCount_PUdown", "; ;Events", 1, 0, 1);
  m_eventCount_IDup = histFolder.make<TH1F>("eventCount_IDup", "; ;Events", 1, 0, 1);
  m_eventCount_IDdown = histFolder.make<TH1F>("eventCount_IDdown", "; ;Events", 1, 0, 1);
  m_eventCount_ISOup = histFolder.make<TH1F>("eventCount_ISOup", "; ;Events", 1, 0, 1);
  m_eventCount_ISOdown = histFolder.make<TH1F>("eventCount_ISOdown", "; ;Events", 1, 0, 1);
  m_eventCount_Trigup = histFolder.make<TH1F>("eventCount_Trigup", "; ;Events", 1, 0, 1);
  m_eventCount_Trigdown = histFolder.make<TH1F>("eventCount_Trigdown", "; ;Events", 1, 0, 1);
  m_eventCount_EnBin = histFolder.make<TH1F>("eventCount_EnBin","; ;Events", 1,0,1);
  m_muRocWeight = histFolder.make<TH1F>("MuRocWeight",";Muon Rochester Correction Weight; Events",100,0.9,1.1);
  m_IDSF = histFolder.make<TH1F>("IDSF",";Tight muon ID scale factor; Events",100,0.9,1.1);
  m_IDSFup = histFolder.make<TH1F>("IDSFup",";Tight muon ID scale factor up; Events",100,0.9,1.1);
  m_IDSFdown = histFolder.make<TH1F>("IDSFdown",";Tight muon ID scale factor down; Events",100,0.9,1.1);
  m_ISOSF = histFolder.make<TH1F>("ISOSF",";Tight muon ISO scale factor; Events",100,0.9,1.1);
  m_ISOSFup = histFolder.make<TH1F>("ISOSFup",";Tight muon ISO scale factor up; Events",100,0.9,1.1);
  m_ISOSFdown = histFolder.make<TH1F>("ISOSFdown",";Tight muon ISO scale factor down; Events",100,0.9,1.1);
  m_TrigSF = histFolder.make<TH1F>("TrigSF",";Tight muon Trigger scale factor; Events",100,0.9,1.1);
  m_TrigSFup = histFolder.make<TH1F>("TrigSFup",";Tight muon Trigger scale factor up; Events",100,0.9,1.1);
  m_TrigSFdown = histFolder.make<TH1F>("TrigSFdown",";Tight muon Trigger scale factor down; Events",100,0.9,1.1);

  m_cutProgress = histFolder.make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 40, -.5, 39.5);
  m_signalCategory = histFolder.make<TH1F>("signalCategory",";Signal Category; Events", 4,-0.5,3.5);
  m_DiMuonMass = histFolder.make<TH1F>("DiMuonInvariantMass", "; Dimuon Invariant Mass (GeV); Events", 100, 50, 150);
  m_TagMuonEta = histFolder.make<TH1F>("TaggingMuonEta", "; Tagging Muon #eta; Events", 100, -2.6, 2.6);
  m_TagMuonPt = histFolder.make<TH1F>("TaggingMuonPt","; Tagging Muon pt (GeV); Events", 200, 0, 120);
  m_nSelectedMuons = histFolder.make<TH1F>("NumberOfSelectedMuons","; # of Selected Muons; Events", 20, -0.5, 19.5);
  m_nMatchedMuons = histFolder.make<TH1F>("NumberOfPairedProbes","; # of Paired Probes; Events", 20, -0.5, 19.5);

  m_ConeHits = histFolder.make<TH1F>("ConeHits", "; Hits in Cone; Events", 50,-0.5,49.5);
  m_ConeEnergy = histFolder.make<TH1F>("ConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_PUup = histFolder.make<TH1F>("ConeEnergy_PUup", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_PUdown = histFolder.make<TH1F>("ConeEnergy_PUdown", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_IDup = histFolder.make<TH1F>("ConeEnergy_IDup", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_IDdown = histFolder.make<TH1F>("ConeEnergy_IDdown", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_ISOup = histFolder.make<TH1F>("ConeEnergy_ISOup", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_ISOdown = histFolder.make<TH1F>("ConeEnergy_ISOdown", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_Trigup = histFolder.make<TH1F>("ConeEnergy_Trigup", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_Trigdown = histFolder.make<TH1F>("ConeEnergy_Trigdown", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_ConeEnergy_EnBin = histFolder.make<TH1F>("ConeEnergy_EnBin","; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_NMatchNHitHCALEnergy = histFolder.make<TH1F>("NMatchNHitConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_NMatchNHitLocation = histFolder.make<TH2F>("NMatchNHitLocation",";Z;#rho;Events",150,0,600,150,0,600);
  m_ProbeTrackPt = histFolder.make<TH1F>("ProbeTrackPt",";#Pt (GeV); Events", 200, 0, 100);
  m_ProbeTrackEta = histFolder.make<TH1F>("ProbeTrackEta",";#eta; Events", 200, -4, 4);
  m_ProbeTrackEtaPhi = histFolder.make<TH2F>("ProbeTrackEtaPhi",";#eta;#phi;Events",100,-2.6,2.6,100,-3.15,3.15);
  m_DBremLocation = histFolder.make<TH2F>("DarkBremLocation",";Z;#rho;Events",150,0,600,150,0,600);
  m_DeflectingAngle = histFolder.make<TH1F>("DeflectingAngle",";Muon Deflection Angle;Events",100,0,3.16);
  m_initialMuE = histFolder.make<TH1F>("MuInitialE",";Energy (GeV);Events",800,0,200);
  m_finalMuE = histFolder.make<TH1F>("MuFinalE",";Energy (GeV);Events",800,0,200);
  m_DBremScalingDistance = histFolder.make<TH1F>("DarkBremScalingDistance",";Initial E/Sampled E; Events", 100,0,4);
  m_ProbeEcalIsolation = histFolder.make<TH1F>("ProbeEcalIsolation",";Ecal Isolation; Events", 100,0,100);
  m_IsolationVariables = histFolder.make<TH2F>("IsolationVars",";ECAL Isolation;Track Isolation; Events",6,0,30,8,0,0.2);
  m_TrackMotherdE = histFolder.make<TH1F>("MotherTrackdE",";#Delta E /E ;Events",100,0,1);
  m_TrackMotherdR = histFolder.make<TH1F>("MotherTrackdR",";#Delta R ;Events",50,0,0.5);
  m_FailLocation = histFolder.make<TH2F>("FailLocation",";i#eta; i#phi; Events", 57, -30.5, 30.5, 76, 0.5, 72.5);
  m_FractionalELost = histFolder.make<TH1F>("FractionalELost",";Fraction of Energy Remaining;Events",100,0,1);
  m_HitsOverThresh = histFolder.make<TH1F>("HitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);
  m_BremDepth = histFolder.make<TH1F>("BremDepth", "; Depth of Dark Brem; Events", 8,-0.5,7.5);
  m_CSC_dR = histFolder.make<TH1F>("dRtoNearestCSC","; dR from Projected Track to Nearest CSC hit; Events", 400,0,4);
  m_MuonMotherdE = histFolder.make<TH1F>("MuonMotherdE","; #Delta E / E; Events",100,0,1);
  m_MuonMotherdR = histFolder.make<TH1F>("MuonMotherdR","; #Delta R; Events",200,0,4);
  m_CSCHitChiSq = histFolder.make<TH1F>("CSCChiSquared","; #Chi^2; Events",200,0,10);
  m_CSCHitAngleChange = histFolder.make<TH1F>("CSCDeflection","; #Theta; Events", 200, 0, 3.15);
  m_StandaloneMuonE = histFolder.make<TH1F>("StandaloneMuonE","; Energy (GeV); Events", 200, 0, 200);
  m_GlobalMuonE = histFolder.make<TH1F>("GlobaldMuonE","; #Delta E over E; Events",200,0,200);
  m_StandalonedE = histFolder.make<TH1F>("StandaloneMuondEoverE","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_PUup = histFolder.make<TH1F>("StandaloneMuondEoverE_PUup","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_PUdown = histFolder.make<TH1F>("StandaloneMuondEoverE_PUdown","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_IDup = histFolder.make<TH1F>("StandaloneMuondEoverE_IDup","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_IDdown = histFolder.make<TH1F>("StandaloneMuondEoverE_IDdown","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_ISOup = histFolder.make<TH1F>("StandaloneMuondEoverE_ISOup","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_ISOdown = histFolder.make<TH1F>("StandaloneMuondEoverE_ISOdown","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_Trigup = histFolder.make<TH1F>("StandaloneMuondEoverE_Trigup","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_Trigdown = histFolder.make<TH1F>("StandaloneMuondEoverE_Trigdown","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1); 
  m_StandalonedE_EnBin = histFolder.make<TH1F>("StandaloneMuondEoverE_EnBin","; (Standalone E - Probe E)/Probe E; Events", 200, -1,1);
  m_MatchedLocation = histFolder.make<TH2F>("MatchedLocation",";#eta; #phi; Events", 200, -2.4,2.4, 100, -3.2, 3.2);
  m_AllVertexOffset = histFolder.make<TH1F>("AllClosestApproach","; Distance of Closest Approach; Events", 200, 0, 10);
  m_EventWeights = histFolder.make<TH1F>("EventWeights","; Weight; Number of Events", 200, 0, 10);
  m_DBremR = histFolder.make<TH1F>("Dark Brem R","; Dark Brem R; Events", 300,0, 600);
  m_HEHitByDepth = histFolder.make<TH1F>("HEHitsByDepth","; Depth; Number of Events With a Hit in this depth", 8, -0.5, 7.5);
  m_DrtoStandalone = histFolder.make<TH1F>("DrtoNearestStandalone","; #Delta R to Nearest Standalone muon; Number of Events", 140, 0, 7);
  m_ProbePFJetDr = histFolder.make<TH1F>("ProbePFJetDr",";#Delta R from Probe to Nearest PF Jet; Events",140,0,7);
  m_ProbePFJetE = histFolder.make<TH1F>("ProbePFJetE",";Energy of Nearest PF Jet to Probe (GeV); Events",100,0,100);
  m_ProbePFJetMaxE = histFolder.make<TH1F>("ProbePFJetMaxE",";Energy of Largest PF Jet Nearby Probe (GeV); Events",100,0,100);
  m_CaloJetMinDr = histFolder.make<TH1F>("MinCaloJetDr",";#Delta R from Probe to Nearest Calo Jet; Events",140,0,7);
  m_CaloJetEcalE = histFolder.make<TH1F>("CaloJetEcalE",";Calo Jet Ecal E (GeV); Events", 100,0,100);
  m_CaloJetHcalE = histFolder.make<TH1F>("CaloJetHcalE",";Calo Jet Ecal E (GeV); Events", 100,0,100);
  m_CaloJetTotalE = histFolder.make<TH1F>("CaloJetTotalE",";Calo Jet Total E (GeV); Events", 100,0,100);
  m_CaloConeEnergy = histFolder.make<TH1F>("CaloConeEnergy",";Calo Cone Total E (geV); Events", 100,0,100);
  m_CaloConeHits = histFolder.make<TH1F>("CaloConeHits",";Calo Cone Hits; Events",20,-0.5,19.5);

  m_HECellsFound = histFolder.make<TH1F>("NHECellsFound","; Cells found in HE; Number of Events", 30,0,30);
  m_HEEnergyStandaloneDEoverE = histFolder.make<TH2F>("HEEStandDEOverE","; Energy in HE Hits; (Standalone P - Probe P)/Probe P; Events",30,0,10,30,-1,1);

  m_BigEtadEtadPhi = histFolder.make<TH2F>("BigEtaCellCenterDistance","; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtadEtadPhi = histFolder.make<TH2F>("SmallEtaCellCenterDistance","; #Delta#eta; #Delta#phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
/*  m_BigEtaAdjFaildEtadPhiMuMinusEtaPlus = histFolder.make<TH2F>("BMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAdjFaildEtadPhiMuMinusEtaMinus = histFolder.make<TH2F>("BMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuPlusEtaPlus = histFolder.make<TH2F>("BAllMuPlusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuPlusEtaMinus = histFolder.make<TH2F>("BAllMuPlusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuMinusEtaPlus = histFolder.make<TH2F>("BAllMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuMinusEtaMinus = histFolder.make<TH2F>("BAllMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuPlusEtaPlus = histFolder.make<TH2F>("SMuPlusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuPlusEtaMinus = histFolder.make<TH2F>("SMuPlusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuMinusEtaPlus = histFolder.make<TH2F>("SMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuMinusEtaMinus = histFolder.make<TH2F>("SMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuPlusEtaPlus = histFolder.make<TH2F>("SAllMuPlusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuPlusEtaMinus = histFolder.make<TH2F>("SAllMuPlusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuMinusEtaPlus = histFolder.make<TH2F>("SAllMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuMinusEtaMinus = histFolder.make<TH2F>("SAllMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);*/
  m_MinCSCImpactParameter = histFolder.make<TH1F>("CSCMinImpact", "; CSC Min Impact Parameter; Events", 300, 0, 1000);
  m_PairVtxChi = histFolder.make<TH1F>("PairedVtxReducedChiSq", "; Paired Vertex Reduced #Chi Squared; Events", 100, 0, 10);
  m_TrackIsolation = histFolder.make<TH1F>("TrackIsolation","; Track Based Isolation; Events",200,0,5);
  m_TrackQuality = histFolder.make<TH1F>("TrackQuality","; Track Quality; N Tracks",3,-0.5,2.5);
  m_TrackReducedChi = histFolder.make<TH1F>("TrackReducedChi","; Track Reduced #chi^{2}; Events",100,0,30);
  m_HcalPlusEcal = histFolder.make<TH1F>("HcalPlusEcalE","; Ecal Cone + Hcal cone (GeV); Events",100,0,100);
  m_Probedxy = histFolder.make<TH1F>("ProbeTrackdxy","; Probe Track dxy; Events",200,-0.3,0.3);
  m_Probedsz = histFolder.make<TH1F>("ProbeTrackdsz","; Probe Track dsz; Events",200,-7,7);
  m_ProbeDCA = histFolder.make<TH1F>("ProbeTrackDCA","; Probe Track Transverse DCA; Events",200,0,0.05);
  m_ProbeLongDCA = histFolder.make<TH1F>("ProbeTrackLongDCA","; Probe Track Long DCA; Events",200,0,0.3);
  m_ProbeDCAphi = histFolder.make<TH2F>("ProbeTrackDCAphi","; Probe Track DCA; Probe Track Phi; Events", 200,0,0.4,100,-3.14,3.14);
  m_nEle = histFolder.make<TH1F>("nElectrons",";Number of Electrons > 10GeV;Events",20,0,20);
  m_nExtraEle = histFolder.make<TH1F>("nExtraElectrons",";Number of Extra Electrons;Events",8,-0.5,7.5);
  m_nExtraMuon = histFolder.make<TH1F>("nExtraMuons",";Number of Extra Muons;Events",8,-0.5,7.5);
  m_standaloneChi2 = histFolder.make<TH1F>("standaloneChi2",";Standalone #Chi^2;Events",200,0,100);
  m_standaloneNHits = histFolder.make<TH1F>("standaloneNHits",";Standalone N Hits;Events",100,0,20);
  m_staTransDCA = histFolder.make<TH1F>("standaloneTDCA",";Standalone transverse DCA;Events",100,0,5);
  m_staLongDCA=histFolder.make<TH1F>("standaloneLongDCA",";Standalone longitudinal DCA;Events",100,0,5);
  m_standaloneChi2dEoverE = histFolder.make<TH2F>("standalonechidEoverE",";Standalone #Chi^2;Standalone dP over P;Events",50,0,100,100,-1,1);
  m_genProbedEoverE = histFolder.make<TH1F>("GenStandaloneMuondEoverE","; (Standalone E - PostBrem E)/PostBrem E; Events", 200, -1,1); 
  m_disGenProbedEoverE = histFolder.make<TH1F>("DisGenStandaloneMuondEoverE","; (Displaced Standalone E - PostBrem E)/PostBrem E; Events", 200, -1,1); 
  m_PileupWeights = histFolder.make<TH1F>("PileupWeights",";N Pileup Interactions; Weight", 100,0,4); 
  m_PUmean = histFolder.make<TH1F>("NPUMean",";MC N Pileup Mean; Events", 100,0,100);
  m_energyBinSFs = histFolder.make<TH1F>("energyBinSF",";SF from Energy Binning; Events", 50,0.5,1.5);
  m_usedEnergyBin = histFolder.make<TH1F>("energyBin",";Energy used for DBrem scaling (GeV); Events",1000,0,1000);
}

void MCHistograms::FillHists(SigEventInfo info)
{
   double eventWeightPUup = info.eventWeight*info.PuUpDownWeights[0]/info.pileupWeight; 
   double eventWeightPUdown = info.eventWeight*info.PuUpDownWeights[1]/info.pileupWeight; 
   double eventWeightIDup = info.eventWeight*info.muIdWeights[1]/info.muIdWeights[0];
   double eventWeightIDdown = info.eventWeight*info.muIdWeights[2]/info.muIdWeights[0];
   double eventWeightISOup = info.eventWeight*info.muIsoWeights[1]/info.muIsoWeights[0];
   double eventWeightISOdown = info.eventWeight*info.muIsoWeights[2]/info.muIsoWeights[0];
   double eventWeightTrigup = info.eventWeight*info.muTrigWeights[1]/info.muTrigWeights[0];
   double eventWeightTrigdown = info.eventWeight*info.muTrigWeights[2]/info.muTrigWeights[0];
   double eventWeightEnBin = info.eventWeight*info.energyBinSF;
   m_PileupWeights->Fill(info.pileupWeight);
   m_PUmean->Fill(info.nPUmean);
   m_genEventCount->Fill(0.5,info.genEventWeight);
   m_genEventCount_EnBin->Fill(0.5,info.genEventWeight*info.energyBinSF);
   m_eventCount->Fill(0.5,info.eventWeight);
   m_eventCount_PUup->Fill(0.5,eventWeightPUup);
   m_eventCount_PUdown->Fill(0.5,eventWeightPUdown);
   m_eventCount_IDup->Fill(0.5,eventWeightIDup);
   m_eventCount_IDdown->Fill(0.5,eventWeightIDdown);
   m_eventCount_ISOup->Fill(0.5,eventWeightISOup);
   m_eventCount_ISOdown->Fill(0.5,eventWeightISOdown);
   m_eventCount_Trigup->Fill(0.5,eventWeightTrigup);
   m_eventCount_Trigdown->Fill(0.5,eventWeightTrigdown);
   m_eventCount_EnBin->Fill(0.5,eventWeightEnBin);
   m_muRocWeight->Fill(info.muRocWeight,info.eventWeight);
   m_IDSF->Fill(info.muIdWeights[0],info.eventWeight);
   m_IDSFup->Fill(info.muIdWeights[1],eventWeightIDup);
   m_IDSFdown->Fill(info.muIdWeights[2],eventWeightIDdown);
   m_ISOSF->Fill(info.muIsoWeights[0],info.eventWeight);
   m_ISOSFup->Fill(info.muIsoWeights[1],eventWeightISOup);
   m_ISOSFdown->Fill(info.muIsoWeights[2],eventWeightISOdown);
   m_TrigSF->Fill(info.muTrigWeights[0],info.eventWeight);
   m_TrigSFup->Fill(info.muTrigWeights[1],eventWeightTrigup);
   m_TrigSFdown->Fill(info.muTrigWeights[2],eventWeightTrigdown);
   m_EventWeights->Fill(info.eventWeight);
   m_signalCategory->Fill(info.sigCat,info.eventWeight);
   m_DBremLocation->Fill(info.vtxZ,info.vtxPerp,info.eventWeight);
   m_DBremR->Fill(info.bremR,info.eventWeight); 
   m_DBremScalingDistance->Fill(info.dbremInitialE/info.usedBin,info.eventWeight);
   m_nSelectedMuons->Fill(info.nTagMuons,info.eventWeight);
   m_nMatchedMuons->Fill(info.nPairedProbes,info.eventWeight);
   m_TrackIsolation->Fill(info.probeTrackIso,info.eventWeight);
   m_HcalPlusEcal->Fill(info.caloConeE+info.probeEcalIso,info.eventWeight);
   m_Probedxy->Fill(info.dxy,info.eventWeight);
   m_Probedsz->Fill(info.dsz,info.eventWeight);
   m_ProbeDCA->Fill(info.dca,info.eventWeight);
   m_ProbeLongDCA->Fill(info.dcal,info.eventWeight);
   m_ProbeDCAphi->Fill(info.dca,info.probeTrackPhi,info.eventWeight);
   if(info.isLoose){m_TrackQuality->Fill(0.1,info.eventWeight);}
   if(info.isTight){m_TrackQuality->Fill(1.1,info.eventWeight);}
   if(info.isHighPurity){m_TrackQuality->Fill(2.1,info.eventWeight);}
   m_TrackReducedChi->Fill(info.probeReducedChi,info.eventWeight);
   m_TagMuonEta->Fill(info.tagMuonEta,info.eventWeight);
   m_PairVtxChi->Fill(info.pairVtxChi,info.eventWeight);
   m_ProbeEcalIsolation->Fill(info.probeEcalIso,info.eventWeight);
   m_IsolationVariables->Fill(info.probeEcalIso, info.probeTrackIso,info.eventWeight);
   m_TagMuonPt->Fill(info.tagMuonPt,info.eventWeight); 
   m_ProbeTrackPt->Fill(info.probeTrackPt,info.eventWeight);
   m_ProbeTrackEta->Fill(info.probeTrackEta,info.eventWeight);
   m_ProbeTrackEtaPhi->Fill(info.probeTrackEta,info.probeTrackPhi,info.eventWeight);
   m_DiMuonMass->Fill(info.diMuonMass,info.eventWeight);
   m_DeflectingAngle->Fill(info.deflectingAngle,info.eventWeight);
   m_FractionalELost->Fill(info.finalMuE/info.dbremInitialE,info.eventWeight);
   m_initialMuE->Fill(info.dbremInitialE,info.eventWeight);
   m_finalMuE->Fill(info.finalMuE,info.eventWeight);
   m_MinCSCImpactParameter->Fill(info.minCSCImpactParameter,info.eventWeight);
   m_CSC_dR->Fill(info.cscMinDR,info.eventWeight);
   m_CSCHitChiSq->Fill(info.cscChiSq,info.eventWeight);
   m_nEle->Fill(info.nEle,info.eventWeight);
   m_nExtraEle->Fill(info.nExtraEle,info.eventWeight);
   m_nExtraMuon->Fill(info.nExtraMuon,info.eventWeight);
   if(info.cscHitAngleChange<10)
   {
      m_CSCHitAngleChange->Fill(info.cscHitAngleChange,info.eventWeight);
      m_MatchedLocation->Fill(info.probeTrackEta,info.probeTrackPhi,info.eventWeight);   
      m_MuonMotherdE->Fill(info.probeMuondE,info.eventWeight);
      m_MuonMotherdR->Fill(info.probeMuondR,info.eventWeight);
      m_GlobalMuonE->Fill(info.globalMuonE,info.eventWeight);
   }
   if(info.standaloneE>0)
   {
      m_StandaloneMuonE->Fill(info.standaloneE,info.eventWeight);
      m_HEEnergyStandaloneDEoverE->Fill(info.coneEnergy,info.standaloneDEoverE,info.eventWeight); 
      m_StandalonedE->Fill(info.standaloneDEoverE,info.eventWeight);
      m_StandalonedE_PUup->Fill(info.standaloneDEoverE,eventWeightPUup);
      m_StandalonedE_PUdown->Fill(info.standaloneDEoverE,eventWeightPUdown);
      m_StandalonedE_IDup->Fill(info.standaloneDEoverE,eventWeightIDup);
      m_StandalonedE_IDdown->Fill(info.standaloneDEoverE,eventWeightIDdown);
      m_StandalonedE_ISOup->Fill(info.standaloneDEoverE,eventWeightISOup);
      m_StandalonedE_ISOdown->Fill(info.standaloneDEoverE,eventWeightISOdown);
      m_StandalonedE_Trigup->Fill(info.standaloneDEoverE,eventWeightTrigup);
      m_StandalonedE_Trigdown->Fill(info.standaloneDEoverE,eventWeightTrigdown);
      m_StandalonedE_EnBin->Fill(info.standaloneDEoverE,eventWeightEnBin);
   }
   m_AllVertexOffset->Fill(info.closestApproach,info.eventWeight);  
   m_HECellsFound->Fill(info.nCellsFound,info.eventWeight);
   if(fabs(info.iEta)>20){m_BigEtadEtadPhi->Fill(info.cellCenterDEta,info.cellCenterDPhi,info.eventWeight);}
   else{m_SmallEtadEtadPhi->Fill(info.cellCenterDEta,info.cellCenterDPhi,info.eventWeight);}
   if(info.staMinDr>0.5&&info.cscHitAngleChange>10)
   {
      m_NMatchNHitHCALEnergy->Fill(info.coneEnergy,info.eventWeight);
      if(info.coneEnergy<10){m_NMatchNHitLocation->Fill(info.probeTrackEta,info.probeTrackPhi,info.eventWeight);}
   }   
   for(int i=0;i<7;i++)
   {
      if(info.hitEnergies[i]>0.1){m_HEHitByDepth->Fill(i+1,info.eventWeight);}
   }
   m_DrtoStandalone->Fill(info.staMinDr,info.eventWeight);
   m_ProbePFJetE->Fill(info.nearestJetE,info.eventWeight);
   m_ProbePFJetDr->Fill(info.nearestJetDr,info.eventWeight);
   m_ProbePFJetMaxE->Fill(info.nearJetMaxE,info.eventWeight);
   m_CaloJetMinDr->Fill(info.minCaloJetDr,info.eventWeight);
   m_CaloJetEcalE->Fill(info.caloJetEcalE,info.eventWeight);
   m_CaloJetHcalE->Fill(info.caloJetHcalE,info.eventWeight);
   m_CaloJetTotalE->Fill(info.caloJetTotalE,info.eventWeight);
   m_CaloConeEnergy->Fill(info.caloConeE,info.eventWeight);
   m_CaloConeHits->Fill(info.caloConeHits,info.eventWeight);
   m_TrackMotherdR->Fill(info.genProbeMatchdR,info.eventWeight);
   m_TrackMotherdE->Fill(info.genProbeMatchdE,info.eventWeight);
   m_BremDepth->Fill(info.bremDepth,info.eventWeight);
   m_ConeEnergy->Fill(info.coneEnergy,info.eventWeight);
   m_ConeEnergy_PUup->Fill(info.coneEnergy,eventWeightPUup);
   m_ConeEnergy_PUdown->Fill(info.coneEnergy,eventWeightPUdown);
   m_ConeEnergy_IDup->Fill(info.coneEnergy,eventWeightIDup);
   m_ConeEnergy_IDdown->Fill(info.coneEnergy,eventWeightIDdown);
   m_ConeEnergy_ISOup->Fill(info.coneEnergy,eventWeightISOup);
   m_ConeEnergy_ISOdown->Fill(info.coneEnergy,eventWeightISOdown);
   m_ConeEnergy_Trigup->Fill(info.coneEnergy,eventWeightTrigup);
   m_ConeEnergy_Trigdown->Fill(info.coneEnergy,eventWeightTrigdown);
   m_ConeEnergy_EnBin->Fill(info.coneEnergy,eventWeightEnBin);
   m_ConeHits->Fill(info.coneHits,info.eventWeight);
   m_HitsOverThresh->Fill(info.hitsOverThresh,info.eventWeight);
   m_standaloneChi2->Fill(info.staChi2,info.eventWeight);
   m_standaloneNHits->Fill(info.staNHits,info.eventWeight);
   m_staTransDCA->Fill(info.staTransDCA,info.eventWeight);
   m_staLongDCA->Fill(info.staLongDCA,info.eventWeight);
   m_standaloneChi2dEoverE->Fill(info.staChi2,info.standaloneDEoverE,info.eventWeight);   
   m_energyBinSFs->Fill(info.energyBinSF,info.eventWeight);
   m_usedEnergyBin->Fill(info.usedBin,info.eventWeight);

   if(info.finalMuE>0)
   {
      if(info.standaloneE>0)
      {
         m_genProbedEoverE->Fill((info.standaloneE-info.finalMuE)/info.finalMuE,info.eventWeight);
      }
      if(info.disStaE>0)
      {
         m_disGenProbedEoverE->Fill((info.disStaE-info.finalMuE)/info.finalMuE,info.eventWeight);
      }
   } 
   else
   {
      if(info.standaloneE>0)
      {
         m_genProbedEoverE->Fill(info.standaloneDEoverE,info.eventWeight);
      }
      if(info.disStaE>0)
      {
         m_disGenProbedEoverE->Fill(info.disStaDEoverE,info.eventWeight);
      }
   }  
}

void MCHistograms::IncCutFlow(double weight)
{
   m_cutProgress->Fill(cutProgress,weight);
   cutProgress++;
   return;
}

void MCHistograms::FailLocation(double ieta, double iphi, double weight)
{
   m_FailLocation->Fill(ieta,iphi,weight);
}

void MCHistograms::ResetCutFlow()
{
   cutProgress = -0.5;
}

