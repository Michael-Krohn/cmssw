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

void MCHistograms::book(edm::Service<TFileService> fs){

  m_eventCount = fs->make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
  m_cutProgress = fs->make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 10, -.5, 9.5);

  m_MuonTrackMass = fs->make<TH1F>("MuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50  , 150  );
  m_SimInvMass = fs->make<TH1F>("SimMuonsTrackMass", "; MuonTrackMass (GeV);Events", 100, 50  , 150  );
  m_DiMuonMass = fs->make<TH1F>("DiMuonInvariantMass", "; Dimuon Invariant Mass (GeV); Events", 100, 70, 110);
  m_TagMuonEta = fs->make<TH1F>("TaggingMuonEta", "; Tagging Muon #eta; Events", 100, -2.6, 2.6);
  m_TagMuonPt = fs->make<TH1F>("TaggingMuonPt","; Tagging Muon pt (GeV); Events", 200, 0, 120);
  m_NoBackupTagMuonEta = fs->make<TH1F>("NoBackupTaggingMuonEta", "; Tagging Muon #eta; Events", 100, -2.6, 2.6);
  m_NoBackupTagMuonPt = fs->make<TH1F>("NoBackupTaggingMuonPt","; Tagging Muon pt (GeV); Events", 200, 0, 120);
  m_nSelectedMuons = fs->make<TH1F>("NumberOfSelectedMuons","; # of Selected Muons; Events", 20, -0.5, 19.5);
  m_ConeHits = fs->make<TH1F>("ConeHits", "; Hits in Cone; Events", 50,-0.5,49.5);
  m_ConeEnergy = fs->make<TH1F>("ConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,20);
  m_MuonEtaDist = fs->make<TH1F>("MuonEtaDist",";#eta; Events", 100,-10,10);
  m_TrackPt = fs->make<TH1F>("SelectedTrackPt",";#Pt (GeV); Events", 200, 0, 100);
  m_MissHitTrackPt = fs->make<TH1F>("MissHitTrackPt",";#Pt (GeV); Events", 200, 0, 100);
  m_DBremLocation = fs->make<TH2F>("DarkBremLocation",";Z;#rho;Events",150,0,600,150,0,600);
  m_WeightedDBremLocation = fs->make<TH2F>("WeightedDarkBremLocation",";Z;#rho;Events",150,0,600,150,0,600);
  m_DeflectingAngle = fs->make<TH1F>("DeflectingAngle",";Muon Deflection Angle;Events",100,0,3.16);
  m_initialMuE = fs->make<TH1F>("MuInitialE",";Energy (GeV);Events",800,0,800);
  m_finalMuE = fs->make<TH1F>("MuFinalE",";Energy (GeV);Events",800,0,800);
  m_initialMuPt = fs->make<TH1F>("MuInitialPt",";Pt (GeV);Events",200,0,400);
  m_finalMuPt = fs->make<TH1F>("MuFinalPt",";Pt (GeV);Events",200,0,400);
  m_dphoEnergy = fs->make<TH1F>("DarkPhotonEnergy",";Energy (GeV);Events",100,0,200);
  m_TrackMotherdE = fs->make<TH1F>("MotherTrackdE",";#Delta E /E ;Events",100,0,1);
  m_TrackMotherdR = fs->make<TH1F>("MotherTrackdR",";#Delta R ;Events",50,0,0.5);
  m_FractionalELost = fs->make<TH1F>("FractionalELost",";Fraction of Energy Remaining;Events",100,0,1);
  m_CSCHits_EtaPhi = fs->make<TH2F>("CSCHits_EtaPhi", "; #eta; #phi; Events", 100, -2.4, 2.4, 72, -ROOT::Math::Pi(), ROOT::Math::Pi());
  m_HitDepth_MuonHCAL = fs->make<TH1F>("HitDepth_MuonHCAL", "; Depth; Events", 40,0,20);
  m_HitsOverThresh = fs->make<TH1F>("HitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);
  m_BremDepth = fs->make<TH1F>("BremDepth", "; Depth of Dark Brem; Events", 8,-0.5,7.5);
  m_WeightedBremDepth = fs->make<TH1F>("WeightedBremDepth", "; Depth of Dark Brem; Events", 8,-0.5,7.5);
  for(int i=0;i<7;i++)
  {
    std::string name = "Depth"+std::to_string(i+1)+"HitsOverThresh";
    std::string specname = "Depth"+std::to_string(i+1)+"Spectra";
    m_HitsOverThreshSplit[i] = fs->make<TH1F>(name.c_str(), "; Hits Over Threshold; Events", 8,-0.5,7.5);
    m_Depth_Spectra[i] = fs->make<TH1F>(specname.c_str(), "; Hit Energy (GeV); Events", 400,0,10);
  }
  m_BremSpectrum = fs->make<TH1F>("BremSpectrum","; Energy of Nearest Depth to Brem; Events", 200,0,10);
  m_CSC_dR = fs->make<TH1F>("dRtoNearestCSC","; dR from Projected Track to Nearest CSC hit; Events", 140,0,7);
  m_NThreshCut = fs->make<TH1F>("NPassAdjacentCut","; Fail cut in bin zero, pass in bin 1; Events", 2,-0.5,1.5);
  m_TrackerVGlobal = fs->make<TH1F>("TrackerVGlobalMuons","; Only global in bin one, only tracker in two, both in three, neither in four; Muons", 4,0.5,4.5);
  m_MuonMotherdE = fs->make<TH1F>("MuonMotherdE","; #Delta E / E; Events",100,0,1);
  m_MuonMotherdR = fs->make<TH1F>("MuonMotherdR","; #Delta R; Events",200,0,4);
  m_CSCHitChiSq = fs->make<TH1F>("CSCChiSquared","; #Chi^2; Events",200,0,10);
  m_CSCHitAngleChange = fs->make<TH1F>("CSCDeflection","; #Theta; Events", 200, 0, 3.15);
  m_StandaloneMuonE = fs->make<TH1F>("StandaloneMuonE","; Energy (GeV); Events", 200, 0, 200);
  m_GlobalMuonE = fs->make<TH1F>("GlobaldMuonE","; #Delta E over E; Events",200,0,200);
  m_StandalonedE = fs->make<TH1F>("StandaloneMuondEoverE","; (Standalone E - PostBrem E)/PostBrem E; Events", 200, -1,1); 
  m_NonMatchedE = fs->make<TH1F>("NoMatchDeflectedE","; Deflected Muon Energy; Events", 200, 0, 200);
  m_NonMatchedLocation = fs->make<TH2F>("NonMatchedLocation",";Z;#rho;Events",150,0,600,150,0,600);
  m_MatchedLocation = fs->make<TH2F>("MatchedLocation",";#eta; #phi; Events", 200, -2.4,2.4, 100, -3.2, 3.2);
  m_NonMatchedLocationEtaPhi = fs->make<TH2F>("NonMatchedEtaPhi",";#eta; #phi; Events", 200, -2.4, 2.4, 100, -3.2, 3.2);
  m_NonMatchedHitsOverThresh = fs->make<TH1F>("NonMatchedHitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);
  m_NonMatchedHEDeposit = fs->make<TH1F>("NonMatchedConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,200);
  m_NonMatchedCSC = fs->make<TH1F>("NonMatcheddRtoNearestCSC","; dR from Projected Track to Nearest CSC hit; Events", 140,0,7);
  m_NHitStandalonedR = fs->make<TH1F>("NHitStandalonedR","; #DeltaR; Events", 200, 0, 4);
  m_NHitStandaloneMuonE = fs->make<TH1F>("NHitStandaloneMuonE","; Energy (GeV); Events", 200, 0, 200);
  m_NHitStandalonedE = fs->make<TH1F>("NHitStandaloneMuondEoverE","; (Standalone E - PostBrem E)/PostBrem E; Events", 200, -1,1); 
  m_BigCSCLocations = fs->make<TH2F>("BigCSCEtaPhi",";#eta; #phi; Events", 200, -2.4, 2.4, 100, -3.2, 3.2);
  m_LowEHCAL = fs->make<TH1F>("SmallMuonHCALEnergies",";HCAL Energy (GeV); Events", 200,0,200);
  m_NStandalone = fs->make<TH1F>("NumberOfStandaloneMuons",";Number of Standalone Muons; Events", 5,-0.5,4.5);
  m_NMatchStandaloneDr = fs->make<TH1F>("NonMatcheddRtoNearestStandalone","; dR from Projected Track to Nearest Standalone Muon; Events", 140,0,7);
  m_NMatchStandaloneDE = fs->make<TH1F>("NonMatchedDeltaEOverENearestStandalone","; #Delta E /E; Events", 200, -1,1);
  m_NMatchStandaloneE = fs->make<TH1F>("NonMatchedStandaloneE","; Energy (GeV); Events", 200, 0, 200);
  m_StandalonePlusHCALDE = fs->make<TH1F>("StandaloneMuonPlusHCALDEOverE","; (Standalone Energy + HCAL deposits)/Truth Energy); Events", 200,-1,1);
  m_NMatchNHitLocation = fs->make<TH2F>("NoMatchNoStandaloneLocation","; #eta; #phi; Events", 200, -2.4, 2.4, 100, -3.2, 3.2);
  m_NMatchNHitEnergy = fs->make<TH1F>("NonMatchedNoStandaloneE","; Energy (GeV); Events", 200, 0, 200);
  m_NMatchNHitHCALEnergy = fs->make<TH1F>("NonMatchNStandaloneConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,200);
  m_NMatchVertexOffset = fs->make<TH1F>("NonMatchClosestApproach","; Distance of Closest Approach; Events", 200, 0, 10);
  m_AllVertexOffset = fs->make<TH1F>("AllClosestApproach","; Distance of Closest Approach; Events", 200, 0, 10);
  m_NStandaloneAdjacentHitEnergies = fs->make<TH1F>("NonMatchAdjacentEnergies","; Energy in Adjacent Hits; Events", 200, 0, 200);  
  m_AdjacentFailHitEnergy = fs->make<TH1F>("FailAdjacentHitEnergy","; Energy in adjacent hit (Gev); Events", 200, 0, 20);
  m_NStandaloneAdjacentFailHitEnergy = fs->make<TH1F>("NStandaloneFailAdjacentHitEnergy","; Energy in adjacent hit (Gev); Events", 200, 0, 20);
  m_EventWeights = fs->make<TH1F>("EventWeights","; Weight; Number of Events", 200, 0, 10);
  m_SignalSelectionCuts = fs->make<TH1F>("SignalCuts","; Event Category; Events", 3, -0.5, 2.5);
  m_WeightedSignalSelectionCuts = fs->make<TH1F>("WeightedSignalCuts","; Event Category; Events", 3, -0.5, 2.5);
  m_DBremR = fs->make<TH1F>("Dark Brem R","; Dark Brem R; Events", 300,0, 600);
  m_WeightedDBremR = fs->make<TH1F>("Weighted Dark Brem R","; Dark Brem R; Events", 300,0, 600);
  m_NMatchNStandaloneHitsOverThresh = fs->make<TH1F>("NSHitsOverThresh","; Hits Over Threshold; Events", 8, -0.5, 7.5);
  m_PassSigEta = fs->make<TH1F>("PassSigEta","; Eta; Number of Events", 100, -2.4, 2.4);
  m_PassSigHEHits = fs->make<TH1F>("PassSigHEHits","; Hits in HE; Number of Events", 8, -0.5, 7.5);
  m_PassSigHEEnergy = fs->make<TH1F>("PassSigHEEnergy","; Cone energy in HE (GeV); Number of Events", 200, 0, 100);
  m_PassSigHEHitByDepth = fs->make<TH1F>("PassSigHEHitsByDepth","; Depth; Number of Events With a Hit in this depth", 8, -0.5, 7.5);
  m_PassSigDrtoCSC = fs->make<TH1F>("PassSigDrtoCSC","; #Delta R to Nearest CSC; Number of Events", 140, 0, 0.5);
  m_PassSigDrtoStandalone = fs->make<TH1F>("PassSigDrtoNearestStandalone","; #Delta R to Nearest Standalone muon; Number of Events", 140, 0, 7);
  m_PassSigLocation = fs->make<TH2F>("PassSigLocation","; #eta; #phi; Events", 200, -2.4, 2.4, 100, -3.2, 3.2);
  m_PassSigVtxZ = fs->make<TH1F>("PassSigVtxZ","; Vertex z position; Events", 100,-20,20);
  m_NPassSigVtxZ = fs->make<TH1F>("NPassSigVtxZ","; Vertex z position; Events", 100,-20,20);
  m_NPassSigLocation = fs->make<TH2F>("NPassSigLocation","; #eta; #phi; Events", 200, -2.4, 2.4, 100, -3.2, 3.2);
  m_NPassSigEta = fs->make<TH1F>("NPassSigEta","; Eta; Number of Events", 100, -2.4, 2.4);
  m_NPassSigHEHits = fs->make<TH1F>("NPassSigHEHits","; Hits in HE; Number of Events", 8, -0.5, 7.5);
  m_NPassSigHEEnergy = fs->make<TH1F>("NPassSigHEEnergy","; Cone energy in HE (GeV); Number of Events", 200, 0, 100);
  m_NPassSigHEHitByDepth = fs->make<TH1F>("NPassSigHEHitsByDepth","; Depth; Number of Events With a Hit in this depth", 8, -0.5, 7.5);
  m_NPassSigDrtoCSC = fs->make<TH1F>("NPassSigDrtoCSC","; #Delta R to Nearest CSC; Number of Events", 140, 0, 0.5);
  m_NPassSigDrtoStandalone = fs->make<TH1F>("NPassSigDrtoNearestStandalone","; #Delta R to Nearest Standalone muon; Number of Events", 140, 0, 7);
  m_FailAdjVtxZ = fs->make<TH1F>("FailAdjVtxZ","; Vertex Z position; Events", 100, -20, 20);
  m_HECellsFound = fs->make<TH1F>("NHECellsFound","; Cells found in HE; Number of Events", 30,0,30);
  m_BigEtaAdjFaildEtadPhiMuPlusEtaPlus = fs->make<TH2F>("BMuPlusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAdjFaildEtadPhiMuPlusEtaMinus = fs->make<TH2F>("BMuPlusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAdjFaildEtadPhiMuMinusEtaPlus = fs->make<TH2F>("BMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAdjFaildEtadPhiMuMinusEtaMinus = fs->make<TH2F>("BMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuPlusEtaPlus = fs->make<TH2F>("BAllMuPlusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuPlusEtaMinus = fs->make<TH2F>("BAllMuPlusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuMinusEtaPlus = fs->make<TH2F>("BAllMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigEtaAlldEtadPhiMuMinusEtaMinus = fs->make<TH2F>("BAllMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuPlusEtaPlus = fs->make<TH2F>("SMuPlusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuPlusEtaMinus = fs->make<TH2F>("SMuPlusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuMinusEtaPlus = fs->make<TH2F>("SMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAdjFaildEtadPhiMuMinusEtaMinus = fs->make<TH2F>("SMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuPlusEtaPlus = fs->make<TH2F>("SAllMuPlusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuPlusEtaMinus = fs->make<TH2F>("SAllMuPlusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuMinusEtaPlus = fs->make<TH2F>("SAllMuMinusEtaPlus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_SmallEtaAlldEtadPhiMuMinusEtaMinus = fs->make<TH2F>("SAllMuMinusEtaMinus","; #eta; #phi; Events", 50, -0.2, 0.2, 50, -0.2, 0.2);
  m_BigCSCdrEta = fs->make<TH1F>("BigCSCdrEta","; #eta; Events", 200, -2.4,2.4);
  m_NMatchedMuons = fs->make<TH1F>("NMatchedMuons", "; Number of tagging muons matching track; Events", 10,-0.5,9.5);
  m_MinCSCImpactParameter = fs->make<TH1F>("CSCMinImpact", "; CSC Min Impact Parameter; Events", 300, 0, 1000);
  m_EventsWithDpho = fs->make<TH1F>("EventsWithDpho", "; Bin; Events", 2,0,2);
  m_PairVtxChi = fs->make<TH1F>("PairedVtxReducedChiSq", "; Paired Vertex Reduced #Chi Squared; Events", 100, 0, 10);
  m_SimVtxChi = fs->make<TH1F>("SimPairedVtxReducedChiSq", "; Paired Vertex Reduced #Chi Squared; Events", 100, 0, 10);
  m_BigConeEnergy = fs->make<TH1F>("BigConeEnergy", "; Energy of hits within #Delta R of 0.3 to Probe Muon; Events", 100, 0, 20);
  m_SimProbeTrackPt = fs->make<TH1F>("SimProbeTrackPt", "; Pt of Probe Track from Z decay; Events", 100, 0, 120);
  m_SimProbeTrackEta = fs->make<TH1F>("SimProbeTrackEta", "; #eta of Potential Probe Tracks; Events", 100, -2.6, 2.6);
  m_ReducedHEConeEnergy = fs->make<TH1F>("ReducedHEConeEnergy", "; Energy of hits >0.8 GeV Along Probe Muon Path; Events", 100, 0, 20);
  m_TrackIsolation = fs->make<TH1F>("TrackIsolation","; Track Based Isolation; Events",200,0,5);
}

void MCHistograms::IncCutFlow()
{
   cutProgress++;
   m_cutProgress->Fill(cutProgress);
   return;
}

void MCHistograms::ResetCutFlow()
{
   cutProgress = -0.5;
}

void MCHistograms::PlotCSCHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label)
{

  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  if ((TheCSCSegments.isValid())&&(m_CSCHits_EtaPhi->GetEntries()==0))
  {
     for(CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end(); iSegment++)
     {
        //CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
	//if(iDetId.station() != 1) continue;
        DetId TheDetUnitId(iSegment->cscDetId());
	const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
	m_CSCHits_EtaPhi->Fill(TheUnit->toGlobal(iSegment->localPosition()).eta(),TheUnit->toGlobal(iSegment->localPosition()).phi());
     }
  }
	
}

