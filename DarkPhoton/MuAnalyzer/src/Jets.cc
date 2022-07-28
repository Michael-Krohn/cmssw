#include "DarkPhoton/MuAnalyzer/interface/Jets.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <TMath.h>

Jets::Jets() {
  nearestJetDr = 1000;
  nearestJetE = -1;
  nearJetMaxE = -1;
  wpLoose=0.1208;
  wpMed=0.4168;
  wpTight=0.7665;
  nJets=0;
  looseBTagJets=0;
  medBTagJets=0;
  tightBTagJets=0;
}

void Jets::SelectJets(const edm::Event& iEvent, edm::EDGetToken m_jetsToken) {
  edm::Handle<reco::PFJetCollection> pfjetHandle;
  iEvent.getByToken(m_jetsToken, pfjetHandle);
  if (!pfjetHandle.isValid()) {
    return;
  }
  const reco::PFJetCollection pfjets = *(pfjetHandle.product());
  for (reco::PFJetCollection::const_iterator it = pfjets.begin(); it != pfjets.end(); ++it) {
    double NHF, NEMF, CHF, MUF, CEMF;
    int CHM, NumNeutralParticles, NumConst;
    NHF = it->neutralHadronEnergyFraction();
    NEMF = it->neutralEmEnergyFraction();
    CHF = it->chargedHadronEnergyFraction();
    MUF = it->muonEnergyFraction();
    CEMF = it->chargedEmEnergyFraction();
    NumConst = it->chargedMultiplicity() + it->neutralMultiplicity();
    NumNeutralParticles = it->neutralMultiplicity();
    CHM = it->chargedMultiplicity();
    bool JetId;
    if (fabs(it->eta()) < 2.6) {
      JetId = (CEMF < 0.8 && CHM > 0 && CHF > 0 && NumConst > 1 && NEMF < 0.9 && MUF < 0.8 && NHF < 0.9);
    } else if (fabs(it->eta()) < 2.7) {
      JetId = (CEMF < 0.8 && CHM > 0 && NEMF < 0.99 && MUF < 0.8 && NHF < 0.9);
    } else if (fabs(it->eta()) < 3.0) {
      JetId = (NEMF > 0.01 && NEMF < 0.99 && NumNeutralParticles > 1);
    } else {
      JetId = (NEMF < 0.90 && NHF > 0.2 && NumNeutralParticles > 10);
    }
   
    if (JetId) {
      selectedJets.push_back(&(*it));
      if(it->pt()>30){nJets++;}
    }
  }
}

void Jets::FindNearestJet(const reco::Track* probe) {
  for (auto it = selectedJets.begin(); it != selectedJets.end(); ++it) {
    double jetDr = deltaR(probe->eta(), probe->phi(), (*it)->eta(), (*it)->phi());
    if (jetDr < 0.3 && (*it)->energy() > nearJetMaxE) {
      nearJetMaxE = (*it)->energy();
    }
    if (jetDr < nearestJetDr || nearestJetDr < 0) {
      nearestJetDr = jetDr;
      nearestJetE = (*it)->energy();
    }
  }
}

void Jets::SelectBTagJets(const edm::Event& iEvent, std::vector<edm::EDGetTokenT<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>, std::vector<float> > > > jetTagTokens)
{
  JetTagMap allJetBTagVals;
  for(const auto& JetTagToken : jetTagTokens)
  {
    edm::Handle<reco::JetTagCollection> bjetHandle;
    iEvent.getByToken(JetTagToken,bjetHandle);
    if (not bjetHandle.isValid()){
      edm::LogWarning("MyJets") << "B-Jet handle not valid, will skip event \n";
      return;
    }
    const reco::JetTagCollection& bTags = *(bjetHandle.product());
    for (const auto& i_jetTag : bTags){
      const auto& jetRef = i_jetTag.first;
      const auto btagVal = i_jetTag.second;
      if (not std::isfinite(btagVal)){
        continue;
      }
      if (allJetBTagVals.find(jetRef)!=allJetBTagVals.end()){
        allJetBTagVals.at(jetRef)+=btagVal;
      }
      else{
        allJetBTagVals.insert(JetTagMap::value_type(jetRef, btagVal));
      }
    }
  }

  if(allJetBTagVals.size()>0){
    for(auto jetBTagVal = allJetBTagVals.begin(); jetBTagVal!=allJetBTagVals.end(); jetBTagVal++)
    {
      auto jet = *dynamic_cast<const reco::Jet*>(jetBTagVal->first.get());
      auto btag = jetBTagVal->second;
      if(jet.pt()<30){continue;}
      if(btag>wpLoose){looseBTagJets++;}
      if(btag>wpMed){medBTagJets++;}
      if(btag>wpTight){tightBTagJets++;}
    }
  }
  return;
}
