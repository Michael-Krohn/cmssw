// -*- C++ -*-
//
// Package:    DarkPhoton/LowBiasAnalyzer
// Class:      LowBiasAnalyzer
// 
/**\class LowBiasAnalyzer LowBiasAnalyzer.cc DarkPhoton/MuAnalyzer/plugins/LowBiasAnalyzer.cc

 Description: [one line class summary]
 Read through low bias DBrem sim file and count events with brem and events without.
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

class LowBiasAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit LowBiasAnalyzer(const edm::ParameterSet&);
      ~LowBiasAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
    edm::EDGetToken m_genParticleToken;
    edm::EDGetToken m_dbremweighttoken;
    edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
    TFileDirectory histFolder;
    TH1F* m_eventCount;
    TH1F* m_eventWeightDBrem;
    TH1F* m_eventWeightNoDBrem;
    TH1F* m_hasDBrem;
    double scaleFactor;
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
LowBiasAnalyzer::LowBiasAnalyzer(const edm::ParameterSet& iConfig):
  scaleFactor(iConfig.getUntrackedParameter<double>("scaleFactor",1))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
   //now do what ever initialization is needed
  m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  
  m_dbremweighttoken = consumes<float>(iConfig.getParameter<edm::InputTag>("dbremweighttag"));
  histFolder = fs->mkdir("allEvents");
  m_eventCount = histFolder.make<TH1F>("EventCount",";Events;Count",1,0,1);
  m_eventWeightDBrem = histFolder.make<TH1F>("EventWeightDBrem",";Event Weight;Count",200,0.,2);
  m_eventWeightNoDBrem = histFolder.make<TH1F>("EventWeightNoDBrem",";Event Weight;Count",200,0.,10);
  m_hasDBrem = histFolder.make<TH1F>("HasDBrem",";Has DBrem;Weight",2,-0.5,1.5);
  m_hasDBrem->GetXaxis()->SetBinLabel(1,"No DBrem");
  m_hasDBrem->GetXaxis()->SetBinLabel(2,"Yes DBrem");
  weight_=0;
}


LowBiasAnalyzer::~LowBiasAnalyzer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void LowBiasAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
  using namespace edm;
  using namespace std;
  edm::Handle<float> dbremweight;
  iEvent.getByToken(m_dbremweighttoken,dbremweight);
  edm::Handle<GenEventInfoProduct> eventInfo;
  iEvent.getByToken(m_genInfoToken, eventInfo);
  weight_  = eventInfo->weight();
  m_eventCount->Fill(0.5,weight_);
  if(*dbremweight>0)
  {
    weight_ = *dbremweight*weight_;
    if(*dbremweight*scaleFactor>10){return;}
    m_eventWeightDBrem->Fill(weight_*scaleFactor);
    m_hasDBrem->Fill(1,weight_);
  }
  else
  {
    m_eventWeightNoDBrem->Fill(weight_);
    m_hasDBrem->Fill(0.,weight_);
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
LowBiasAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LowBiasAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LowBiasAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LowBiasAnalyzer);
