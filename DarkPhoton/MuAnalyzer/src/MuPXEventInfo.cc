#include "FWCore/Framework/interface/one/EDAnalyzer.h"
  
#include "DarkPhoton/MuAnalyzer/interface/MuPXEventInfo.h"
#include <iostream>

MuPXEventInfo::MuPXEventInfo()
{
    eventWeight=1;
    cutProgress=0;
    paired=false;
    muonTrackMass=-1;
    sigCat=-1;
    probeTrack=NULL;
    tagMuon=NULL;
    nPassingProbe=0;
    nPassingTag=0;
    tagProbeVtxChi=-1;
    probeTrackIso=-1;
    probeHcalIso=-1;
    probeEcalIso=-1;
    smallestCone=-1;
    averageDr=-1;
    tagTrackIso=-1;
    tagHcalIso=-1;
    tagEcalIso=-1;
    nJets=0;
    nearestJetE=-1;
    nearestJetDr=-1;
    nearJetMaxE=-1;
    probeCaloJetE=-1;
    tagCaloJetE=-1;
    probeCaloJetEcal=-1;
    probeCaloJetHcal=-1;
    farMuE=0;
    farMuDr=0;
    farMuEta=0;
    genTagMuDr=0;
    genPdgId=0;
    tagPdgId=0;
    passLooseMatch=false;
    nonMuPdgId=0;
    nLooseMuons=0;
}

