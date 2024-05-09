
#include "TauAnalysis/MCEmbeddingTools/plugins/EmbeddingBeamSpotOnlineProducer.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/Framework/interface/MakerMacros.h"


using namespace edm;


EmbeddingBeamSpotOnlineProducer::EmbeddingBeamSpotOnlineProducer(const ParameterSet& iconf)
{
  beamSpotInput_ = iconf.getParameter<edm::InputTag>("src");
  consumes<reco::BeamSpot>(beamSpotInput_);
  produces<reco::BeamSpot>();

} 

EmbeddingBeamSpotOnlineProducer::~EmbeddingBeamSpotOnlineProducer() {}

void
EmbeddingBeamSpotOnlineProducer::produce(Event& iEvent, const EventSetup& iSetup)
{
  std::cout << "Inside produce TauAnalysis/MCEmbeddingTools/plugins/EmbeddingBeamSpotProducer.cc" << std::endl;
  //copy beam spot from input data into HLT simulation sequence
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByLabel(beamSpotInput_,beamSpotH);
  auto result = std::make_unique<reco::BeamSpot>();
  *result = *beamSpotH;
  iEvent.put(std::move(result));
  std::cout << "Done with produce TauAnalysis/MCEmbeddingTools/plugins/EmbeddingBeamSpotProducer.cc" << std::endl;


}

DEFINE_FWK_MODULE(EmbeddingBeamSpotOnlineProducer);
