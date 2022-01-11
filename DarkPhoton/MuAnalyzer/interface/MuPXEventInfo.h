#ifndef MuPXEventInfo_h
#define MuPXEventInfo_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"  
#include "TH1F.h"

class MuPXEventInfo{
  public:
    MuPXEventInfo();
    double eventWeight;
    int cutProgress;
    int sigCat;
    bool paired;
    double muonTrackMass;
    double adjMuonTrackMass;
    const reco::Track* probeTrack;
    const reco::Muon* tagMuon;
    int nPassingProbe;
    int nPassingTag;
    double tagProbeVtxChi;
    double probeTrackIso;
    double probeHcalIso;
    double probeEcalIso;
    //Multiple Paired Track Hists
    double smallestCone;
    double averageDr;
    double tagTrackIso;
    double tagHcalIso;
    double tagEcalIso;
    double probeDCA;
    double probeLongDCA;
    int nJets;
    double nearestJetE;
    double nearestJetDr;   
    double nearJetMaxE;
    //Calo Jets
    double minProbeCaloJetDr;
    double probeCaloJetE;
    double tagCaloJetE;
    double probeCaloJetEcal;
    double probeCaloJetHcal;
    //vtx stuff
    int probeVtx;
    int tagVtx;
    int nVtx;
    double vtxd0;
    double vtxdz;
    //Pileup Weights
    double pileupWeight;
    double nPUmean;
    double farMuE;
    double farMuDr;
    double farMuEta;
    double genTagMuDr;
    double genPdgId;
    bool passLooseMatch;
    std::vector<int> nearbyNonMuPdgIds;
    std::vector<int> highPtNonMuPdgIds;
    int nLooseMuons;
    double tagPdgId;
    double staMinDr;
    double standaloneE;
    int nEle;
    int nExtraEle;
    int nExtraMuon;
};
#endif
