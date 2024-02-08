// -*- C++ -*-
//
// Package:    MuMu/METFilter_AOD
// Class:      METFilter_AOD
// 
/**\class METFilter_AOD METFilter_AOD_AOD.cc MuMu/METFilter_AOD/plugins/METFilter_AOD_AOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Krohn
//         Created:  Mon, 18 Jun 2018 21:22:23 GMT
//
//
// system include files
#include <memory>
#include <iomanip>
#include <iostream>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include "Math/VectorUtil.h"
// user include files
#include "DarkPhoton/MuAnalyzer/interface/metFilterTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/SortedCollection.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class METFilter_AOD : public edm::stream::EDFilter<>  {
   public:
      explicit METFilter_AOD(const edm::ParameterSet&);
      ~METFilter_AOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      
    metFilterTree failingEvents;
    edm::EDGetTokenT<edm::TriggerResults> trigResults_Label;
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
METFilter_AOD::METFilter_AOD(const edm::ParameterSet& iConfig):
  trigResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResultsTag")))
{
  edm::Service<TFileService> fs;
  failingEvents.book(fs->mkdir("failingEvents"));
   //now do what ever initialization is needed
}


METFilter_AOD::~METFilter_AOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
bool
METFilter_AOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  printf("Starting the thing.\n");  
  std::cout << "Printing things.\n" << std::endl;
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trigResults_Label, trigResults);
  printf("Got the trig results.\n");
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
  std::vector<std::string> metFilterPaths= {"Flag_goodVertices",
	                                    "Flag_globalSuperTightHalo2016Filter",
					    "Flag_HBHENoiseFilter",
					    "Flag_HBHENoiseIsoFilter",
					    "Flag_EcalDeadCellTriggerPrimitiveFilter",
					    "Flag_BadPFMuonFilter",
					    "Flag_BadPFMuonDzFilter",
					    "Flag_hfNoisyHitsFilter",
					    "Flag_eeBadScFilter",
					    "Flag_ecalBadCalibFilter"};
  for(uint8_t path=0;path<metFilterPaths.size();path++)
  {
     printf("Going over a path.\n");
     bool foundtrig = false;
     for(unsigned int itrig=0; itrig!=trigResults->size(); ++itrig)
     {
        if(trigNames.triggerName(itrig).find(metFilterPaths[path])!=std::string::npos)
        {
           foundtrig = true;
           bool passTrig=trigResults->accept(itrig);
           if(!passTrig)
           {
              std::cout<<"Failed trigger" << std::endl;
              int eventNum=iEvent.id().event();
	      int lumi=iEvent.id().luminosityBlock();
	      int run=iEvent.id().run();
	      failingEvents.FillTree(run,lumi,eventNum,metFilterPaths[path],true);
	      printf("Failed trigger %s.\n",metFilterPaths[path].c_str());
              return false;
           }
	   else{printf("Passed trigger %s.\n",metFilterPaths[path].c_str());}
        }
     }
     if(!foundtrig)
     {
        int eventNum=iEvent.id().event();
        int lumi=iEvent.id().luminosityBlock();                                            
        int run=iEvent.id().run();
        failingEvents.FillTree(run,lumi,eventNum,metFilterPaths[path],false);
	printf("Failed to find path %s.\n",metFilterPaths[path].c_str());
        return false;
     }
  }
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
METFilter_AOD::beginStream(edm::StreamID)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METFilter_AOD::endStream() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METFilter_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(METFilter_AOD);
