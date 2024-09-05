#ifndef Jets_h
#define Jets_h

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

class Jets {
public:
  Jets();

  void SelectJets(const edm::Event& iEvent, edm::EDGetToken m_jetsToken);
  void FindNearestJet(const reco::Track* probe);
  void SelectBTagJets(
      const edm::Event& iEvent,
      std::vector<edm::EDGetTokenT<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>, std::vector<float> > > >
          jetTagTokens);

  std::vector<const reco::PFJet*> selectedJets;
  double nearestJetDr;
  double nearestJetE;
  double nearJetMaxE;
  struct JetRefCompare {
    inline bool operator()(const edm::RefToBase<reco::Jet>& j1, const edm::RefToBase<reco::Jet>& j2) const {
      return (j1.id() < j2.id()) || ((j1.id() == j2.id()) && (j1.key() < j2.key()));
    }
  };
  typedef std::map<edm::RefToBase<reco::Jet>, float, JetRefCompare> JetTagMap;
  double wpLoose;
  double wpMed;
  double wpTight;
  int nJets;
  int looseBTagJets;
  int medBTagJets;
  int tightBTagJets;
};

#endif
