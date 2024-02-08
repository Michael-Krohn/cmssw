// -*- C++ -*-
//
// Package:    DarkPhoton/SimAnalyzer
// Class:      SimAnalyzer
//
/**\class SimAnalyzer SimAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/SimAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Revering (adapted from analyzer developed by Michael Krohn)
//         Created:  Friday, 2 Jan 2020 13:34:23 GMT
//
//
// system include files
#include <memory>
#include <iomanip>
#include <iostream>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/transform.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

//for Standalone Muon Tracking
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

// for vertexing
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// for Global Coordinates
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

//for HCAL info
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

//Triggers
#include "FWCore/Common/interface/TriggerNames.h"

//Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

//Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DarkPhoton/MuAnalyzer/interface/SimHistograms.h"
#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"
#include "DarkPhoton/MuAnalyzer/interface/ECAL.h"
#include "DarkPhoton/MuAnalyzer/interface/Jets.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SimAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimAnalyzer(const edm::ParameterSet&);
  ~SimAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  edm::EDGetToken m_genParticleToken;
  edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
  edm::EDGetToken m_dbremweighttoken;
  edm::EDGetToken m_dbremInitialEToken;
  edm::EDGetToken m_dbremFinalEToken;
  edm::EDGetToken m_dbremAngleToken;
  edm::EDGetToken m_dbremXToken;
  edm::EDGetToken m_dbremYToken;
  edm::EDGetToken m_dbremZToken;

  double weight_;
  SimHistograms allEvents;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SimAnalyzer::SimAnalyzer(const edm::ParameterSet& iConfig)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  //now do what ever initialization is needed
  m_genParticleToken = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
  m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  m_dbremweighttoken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremweighttag"));
  m_dbremInitialEToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremInitialETag"));
  m_dbremFinalEToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremFinalETag"));
  m_dbremAngleToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremAngleTag"));
  m_dbremXToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremXTag"));
  m_dbremYToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremYTag"));
  m_dbremZToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremZTag"));

  weight_=1.;
  allEvents.book(fs->mkdir("allEvents"));
}

SimAnalyzer::~SimAnalyzer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void SimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  SigEventInfo info;
  edm::Handle<float> dbremweight, dbremInitialE, dbremFinalE, dbremAngle, dbremX, dbremY, dbremZ;
  iEvent.getByToken(m_dbremweighttoken, dbremweight);
  edm::Handle<GenEventInfoProduct> eventInfo;
  iEvent.getByToken(m_genInfoToken, eventInfo);
  weight_ = eventInfo->weight();
  if (*dbremweight > 0) {
    iEvent.getByToken(m_dbremFinalEToken, dbremFinalE);
    iEvent.getByToken(m_dbremInitialEToken, dbremInitialE);
    iEvent.getByToken(m_dbremAngleToken, dbremAngle);
    iEvent.getByToken(m_dbremXToken, dbremX);
    iEvent.getByToken(m_dbremYToken, dbremY);
    iEvent.getByToken(m_dbremZToken, dbremZ);
    weight_ = *dbremweight * weight_;
    info.dbremWeight = *dbremweight;
    info.dbremInitialE = *dbremInitialE;
    info.dbremFinalE = *dbremFinalE;
    double hyp = sqrt(*dbremX*(*dbremX)+*dbremY*(*dbremY)+*dbremZ*(*dbremZ));
    double dt = sqrt(*dbremX*(*dbremX)+*dbremY*(*dbremY));
    double sinTheta = dt/hyp;
    double et = *dbremInitialE*sinTheta;
    info.dbremInitialEt = et; 
    if (*dbremweight > 0.01) {
      return;
    }
    if (*dbremInitialE < 20) {
      return;
    }
  }
  allEvents.FillHists(info);
}

// ------------ method called once each job just before starting event loop  ------------
void SimAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SimAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimAnalyzer);
