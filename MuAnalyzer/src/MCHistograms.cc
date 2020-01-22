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

  m_ConeHits = fs->make<TH1F>("ConeHits", "; Hits in Cone; Events", 50,-0.5,49.5);
  m_ConeEnergy = fs->make<TH1F>("ConeEnergy", "; Energy in Matched Hits (GeV); Events", 200,0,200);
  m_MuonEtaDist = fs->make<TH1F>("MuonEtaDist",";#eta; Events", 100,-10,10);
  m_TrackPt = fs->make<TH1F>("SelectedTrackPt",";#Pt (GeV); Events", 200, 0, 100);
  m_MissHitTrackPt = fs->make<TH1F>("MissHitTrackPt",";#Pt (GeV); Events", 200, 0, 100);
  m_DBremLocation = fs->make<TH2F>("DarkBremLocation",";Z;#rho;Events",150,0,600,150,0,600);
  m_DeflectingAngle = fs->make<TH1F>("DeflectingAngle",";Muon Deflection Angle;Events",100,0,3.16);
  m_initialMuE = fs->make<TH1F>("MuInitialE",";Energy (GeV);Events",100,0,200);
  m_finalMuE = fs->make<TH1F>("MuFinalE",";Energy (GeV);Events",100,0,200);
  m_dphoEnergy = fs->make<TH1F>("DarkPhotonEnergy",";Energy (GeV);Events",100,0,200);
  m_TrackMotherdE = fs->make<TH1F>("MotherTrackdE",";#DELTA E (GeV);Events",50,0,20);
  m_TrackMotherdR = fs->make<TH1F>("MotherTrackdR",";#DELTA R ;Events",50,0,0.5);
  m_FractionalELost = fs->make<TH1F>("FractionalELost",";Fraction of Energy Remaining;Events",100,0,1);
  m_CSCHits_EtaPhi = fs->make<TH2F>("CSCHits_EtaPhi", "; #eta; #phi; Events", 100, -2.4, 2.4, 72, -ROOT::Math::Pi(), ROOT::Math::Pi());
  m_HitDepth_MuonHCAL = fs->make<TH1F>("HitDepth_MuonHCAL", "; Depth; Events", 40,0,20);
  m_HitsOverThresh = fs->make<TH1F>("HitsOverThresh", "; Hits Over Threshold; Events", 8,-0.5,7.5);
  m_BremDepth = fs->make<TH1F>("BremDepth", "; Depth of Dark Brem; Events", 8,-0.5,7.5);
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

