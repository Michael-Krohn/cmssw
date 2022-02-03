// -*- C++ -*-
//
// Package:    DarkPhoton/WeightsFinder
// Class:      WeightsFinder
// 
/**\class WeightsFinder WeightsFinder.cc DarkPhoton/MuAnalyzer/plugins/WeightsFinder.cc

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

//Event Weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class WeightsFinder : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit WeightsFinder(const edm::ParameterSet&);
      ~WeightsFinder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

    edm::EDGetTokenT<GenEventInfoProduct> m_genInfoToken;
    edm::EDGetToken m_PUInfoToken;
    int pos_;
    int neg_;
    edm::Service<TFileService> fs;
    TH1F* trueNumInt_;
    TH1F* totalweights_;
    TH1F* nEffEvents_;
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
WeightsFinder::WeightsFinder(const edm::ParameterSet& iConfig)
{
  m_genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  m_PUInfoToken = consumes< std::vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("edmPileupInfo"));
  usesResource("TFileService");
  pos_=0;
  neg_=0;
}


WeightsFinder::~WeightsFinder()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void WeightsFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GenEventInfoProduct> eventInfo;
  iEvent.getByToken(m_genInfoToken, eventInfo);
  edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct;
  iEvent.getByToken(m_PUInfoToken,hPileupInfoProduct);
  assert(hPileupInfoProduct.isValid());
  const std::vector<PileupSummaryInfo> *inPUInfos = hPileupInfoProduct.product();
  for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = inPUInfos->begin(); itPUInfo!=inPUInfos->end(); ++itPUInfo) 
  {
     if(itPUInfo->getBunchCrossing()==0) 
     {
        trueNumInt_->Fill(itPUInfo->getTrueNumInteractions());
        break;
     }
  }

  totalweights_->Fill(0.5,eventInfo->weight());
  if(eventInfo->weight()>0){pos_++;}
  else{neg_++;}
}

// ------------ method called once each job just before starting event loop  ------------
void 
WeightsFinder::beginJob()
{
   if (!fs){throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");}
   trueNumInt_ = fs->make<TH1F>("TrueNumInteractions","TrueNumInteractions",100,0,100);
   totalweights_ = fs->make<TH1F>("totalEventWeights","TotalEventWeights",1,0,1);
   nEffEvents_ = fs->make<TH1F>("nEffEvents","nEffEvents",4,0,4);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
WeightsFinder::endJob() 
{
   nEffEvents_->SetBinContent(1,pos_-neg_);
   nEffEvents_->SetBinContent(2,pos_);
   nEffEvents_->SetBinContent(3,neg_);
   nEffEvents_->SetBinContent(4,pos_+neg_);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
WeightsFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WeightsFinder);
