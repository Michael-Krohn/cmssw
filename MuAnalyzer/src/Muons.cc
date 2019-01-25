#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include <TMath.h>

Muons::Muons(){}

void Muons::SelectMuons(const edm::Event& iEvent, edm::EDGetToken m_recoMuonToken){

  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {

    //Loose ID
     if (!(iMuon->isPFMuon() && (iMuon->isGlobalMuon() || iMuon->isTrackerMuon()))) continue;

    //PFIso Loose requirement
    if ((iMuon->pfIsolationR04().sumChargedHadronPt + TMath::Max(0., iMuon->pfIsolationR04().sumNeutralHadronEt + iMuon->pfIsolationR04().sumPhotonEt - 0.5*iMuon->pfIsolationR04().sumPUPt))/iMuon->pt() > 0.25) continue;

    //pT above trigger turnon
    if (iMuon->pt() < 26 || fabs(iMuon->eta()) > 2.4) continue;
    
    selectedMuons.push_back(&(*iMuon));

  }

  if(selectedMuons.size() > 0){
    highPtSelectedMuon = selectedMuons[0];
  }

}
