#include "DarkPhoton/MuFilter/interface/Jets.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TMath.h>

Jets::Jets()
{
   nearestJetDr = 1000;
   nearestJetE = -1;
   nearJetMaxE = -1;
}

void Jets::SelectJets(const edm::Event& iEvent, edm::EDGetToken m_jetsToken){

  edm::Handle<reco::PFJetCollection> pfjetHandle;
  iEvent.getByToken(m_jetsToken, pfjetHandle);
  if(!pfjetHandle.isValid()){return;}
  const reco::PFJetCollection pfjets = *(pfjetHandle.product()); 
  for(reco::PFJetCollection::const_iterator it = pfjets.begin(); it!=pfjets.end();++it)
  {  
    double NHF, NEMF, CHF, MUF, CEMF;
    int CHM,NumNeutralParticles,NumConst;
    NHF = it->neutralHadronEnergyFraction();
    NEMF = it->neutralEmEnergyFraction();
    CHF = it->chargedHadronEnergyFraction();
    MUF = it->muonEnergyFraction();
    CEMF = it->chargedEmEnergyFraction();
    NumConst = it->chargedMultiplicity()+it->neutralMultiplicity();
    NumNeutralParticles =it->neutralMultiplicity();
    CHM = it->chargedMultiplicity();  
    bool JetId;
    if(fabs(it->eta())<2.6){JetId=(CEMF<0.8&&CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9);}
    else if(fabs(it->eta())<2.7){JetId= (CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 );}
    else if(fabs(it->eta())<3.0){JetId=  ( NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>1);}
    else {JetId=(NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10);}
    
    if(JetId){selectedJets.push_back(&(*it));}
  }
}

void Jets::FindNearestJet(const reco::Track* probe)
{
   for(auto it = selectedJets.begin(); it!=selectedJets.end();++it)
   {
      double jetDr=deltaR(probe->eta(),probe->phi(),(*it)->eta(),(*it)->phi());
      if(jetDr<0.3&&(*it)->energy()>nearJetMaxE){nearJetMaxE=(*it)->energy();}
      if(jetDr<nearestJetDr||nearestJetDr<0)
      {
         nearestJetDr=jetDr;
         nearestJetE=(*it)->energy();
      }
   }
}
