// -*- C++ -*-
//
// Package:    DarkPhoton/DBremMatChecker
// Class:      DBremMatChecker
// 
/**\class DBremMatChecker DBremMatChecker.cc DarkPhoton/MuAnalyzer/plugins/DBremMatChecker.cc

 Description: [one line class summary]
 Read through DBrem sim file and check which material brem originated in. Only works on specially generated signal files which have the material stored as the initial muon energy.
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Revering (adapted from analyzer developed by Michael Krohn)
//         Created:  Friday, 2 Jan 2020 13:34:23 GMT
//
// system include files
#include <memory>
#include <iomanip>
#include <iostream>

#include <TH1.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// for vertexing
#include "FWCore/Framework/interface/ESHandle.h"

//Triggers
#include "FWCore/Common/interface/TriggerNames.h"

//Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DBremMatChecker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DBremMatChecker(const edm::ParameterSet&);
      ~DBremMatChecker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    edm::EDGetToken m_genParticleToken;
    edm::EDGetToken m_dbremweighttoken;
    edm::EDGetToken m_dbremInitialEToken;
    edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
    TFileDirectory histFolder;
    TH1F* m_DBremMat;
    double weight_;
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
DBremMatChecker::DBremMatChecker(const edm::ParameterSet& iConfig)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
   //now do what ever initialization is needed
  m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  
  m_dbremweighttoken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremweighttag"));
  m_dbremInitialEToken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremInitialETag"));
  histFolder = fs->mkdir("allEvents");
  m_DBremMat = histFolder.make<TH1F>("DBremMaterial",";Material A;Weight",200,0.5,200.5);
  weight_=0;
}


DBremMatChecker::~DBremMatChecker()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void DBremMatChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
  using namespace edm;
  using namespace std;
  edm::Handle<float> dbremweight, dbremInitialE;
  iEvent.getByToken(m_dbremweighttoken,dbremweight);
  edm::Handle<GenEventInfoProduct> eventInfo;
  iEvent.getByToken(m_genInfoToken, eventInfo);
  weight_  = eventInfo->weight();
  if(*dbremweight>0)
  {
    iEvent.getByToken(m_dbremInitialEToken,dbremInitialE);
    weight_ = *dbremweight*weight_;
    if(*dbremweight>0.01){return;}
    double g_per_mol = 6.2425*pow(10,21);
    m_DBremMat->Fill(*dbremInitialE/g_per_mol,weight_);
  }
  else{return;}
}

// ------------ method called once each job just before starting event loop  ------------
void 
DBremMatChecker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DBremMatChecker::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DBremMatChecker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DBremMatChecker);
