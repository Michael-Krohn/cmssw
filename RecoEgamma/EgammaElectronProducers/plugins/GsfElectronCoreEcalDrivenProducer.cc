#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronTools.h"

class GsfElectronCoreEcalDrivenProducer : public edm::global::EDProducer<> {
public:
  static void fillDescriptions(edm::ConfigurationDescriptions&);

  explicit GsfElectronCoreEcalDrivenProducer(const edm::ParameterSet& conf);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  void produceEcalDrivenCore(const reco::GsfTrackRef& gsfTrackRef,
                             reco::GsfElectronCoreCollection& electrons,
                             edm::Handle<reco::TrackCollection> const& ctfTracksHandle) const;

  const bool useGsfPfRecTracks_;

  const edm::EDGetTokenT<reco::GsfPFRecTrackCollection> gsfPfRecTracksToken_;
  const edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> ctfTracksToken_;
  const edm::EDPutTokenT<reco::GsfElectronCoreCollection> putToken_;
};
using namespace reco;

void GsfElectronCoreEcalDrivenProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("gsfPfRecTracks", {"pfTrackElec"});
  desc.add<edm::InputTag>("gsfTracks", {"electronGsfTracks"});
  desc.add<edm::InputTag>("ctfTracks", {"generalTracks"});
  desc.add<bool>("useGsfPfRecTracks", true);
  descriptions.add("ecalDrivenGsfElectronCores", desc);
}

GsfElectronCoreEcalDrivenProducer::GsfElectronCoreEcalDrivenProducer(const edm::ParameterSet& config)
    : useGsfPfRecTracks_(config.getParameter<bool>("useGsfPfRecTracks")),
      gsfPfRecTracksToken_(mayConsume<GsfPFRecTrackCollection>(config.getParameter<edm::InputTag>("gsfPfRecTracks"))),
      gsfTracksToken_(consumes<reco::GsfTrackCollection>(config.getParameter<edm::InputTag>("gsfTracks"))),
      ctfTracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("ctfTracks"))),
      putToken_(produces<reco::GsfElectronCoreCollection>()) {}

void GsfElectronCoreEcalDrivenProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const {
  auto gsfTracksHandle = event.getHandle(gsfTracksToken_);
  auto ctfTracksHandle = event.getHandle(ctfTracksToken_);

  // output
  reco::GsfElectronCoreCollection electrons;

  // loop on ecal driven tracks
  if (useGsfPfRecTracks_) {
    for (auto const& gsfPfRecTrack : event.get(gsfPfRecTracksToken_)) {
      produceEcalDrivenCore(gsfPfRecTrack.gsfTrackRef(), electrons, ctfTracksHandle);
    }
  } else {
    for (unsigned int i = 0; i < gsfTracksHandle->size(); ++i) {
      produceEcalDrivenCore(edm::Ref<GsfTrackCollection>(gsfTracksHandle, i), electrons, ctfTracksHandle);
    }
  }

  event.emplace(putToken_, electrons);
}

void GsfElectronCoreEcalDrivenProducer::produceEcalDrivenCore(
    const GsfTrackRef& gsfTrackRef,
    GsfElectronCoreCollection& electrons,
    edm::Handle<TrackCollection> const& ctfTracksHandle) const {
  electrons.emplace_back(gsfTrackRef);
  auto& eleCore = electrons.back();

  if (!eleCore.ecalDrivenSeed()) {
    electrons.pop_back();
    return;
  }

  auto ctfpair = egamma::getClosestCtfToGsf(eleCore.gsfTrack(), ctfTracksHandle);
  eleCore.setCtfTrack(ctfpair.first, ctfpair.second);

  auto scRef = gsfTrackRef->extra()->seedRef().castTo<ElectronSeedRef>()->caloCluster().castTo<SuperClusterRef>();
  if (!scRef.isNull()) {
    eleCore.setSuperCluster(scRef);
  } else {
    electrons.pop_back();
    edm::LogWarning("GsfElectronCoreEcalDrivenProducer") << "Seed CaloCluster is not a SuperCluster, unexpected...";
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GsfElectronCoreEcalDrivenProducer);
