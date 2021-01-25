#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include <TMath.h>

Muons::Muons()
{
  highmuonpt=0;
  secondhighmuonpt=0;
}

void Muons::SelectMuons(const edm::Event& iEvent, edm::EDGetToken m_recoMuonToken){

  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  for(std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {

    //Tight ID
//    if (!(iMuon->isPFMuon() && iMuon->isGlobalMuon() )) continue;
    if(!(iMuon->passed(reco::Muon::CutBasedIdTight))) continue;
    //PFIso Loose requirement
    if(!(iMuon->passed(reco::Muon::PFIsoTight))) continue;
    //Track iso
    if(!(iMuon->passed(reco::Muon::TkIsoTight))) continue;
//    if ((iMuon->pfIsolationR04().sumChargedHadronPt + TMath::Max(0., iMuon->pfIsolationR04().sumNeutralHadronEt + iMuon->pfIsolationR04().sumPhotonEt - 0.5*iMuon->pfIsolationR04().sumPUPt))/iMuon->pt() > 0.25) continue;

    //pT above trigger turnon
    if (iMuon->pt() < 26 || fabs(iMuon->eta()) > 2.4) continue;
     
    selectedMuons.push_back(&(*iMuon));

    if(selectedMuons.size() == 1){
      highPtSelectedMuon = selectedMuons[0];
      highmuonpt=iMuon->pt();
      secondhighPtSelectedMuon = selectedMuons[0];
      secondhighmuonpt=0;
    }
    else
    {
      if(iMuon->pt()>highmuonpt)
      {
        secondhighPtSelectedMuon=highPtSelectedMuon;
        secondhighmuonpt=iMuon->pt();
        highPtSelectedMuon=selectedMuons.back();
        highmuonpt=iMuon->pt(); 
      }
      else if(iMuon->pt()>secondhighmuonpt)
      {
        secondhighPtSelectedMuon=selectedMuons.back();
        secondhighmuonpt=iMuon->pt(); 
      }
    }
    if(selectedEndcapMuons.size()==0){highendcappt=0;}
    if (fabs(iMuon->eta())>1.653)
    {
       selectedEndcapMuons.push_back(&(*iMuon));
       if(iMuon->pt()>highendcappt)
       {
         highendcappt=iMuon->pt();
         highPtSelectedEndcapMuon=selectedMuons.back();
       }
    }
  }

}
