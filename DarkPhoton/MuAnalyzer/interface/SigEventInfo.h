#ifndef SigEventInfo_h
#define SigEventInfo_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"  
#include "TH1F.h"

class SigEventInfo{
  public:
    SigEventInfo();
    double eventWeight;
    int sigCat;
    double dbremZ;
    double dbremPerp;
    double bremR;
    double dbremInitialE;
    double dbremFinalE;
    double deflectingAngle;
    int nTagMuons;
    int nPairedProbes;
    double probeTrackIso;
    double probeEcalIso;
    bool isLoose;
    bool isTight;
    bool isHighPurity;
    double probeReducedChi;
    double tagMuonEta;
    double tagMuonPt;
    double probeTrackEta;
    double probeTrackPhi;
    double probeTrackPt;
    double pairVtxChi;
    double diMuonMass;
    //CSC parameters
    double minCSCImpactParameter;
    double cscMinDR;
    double cscChiSq;
    double cscHitAngleChange;
    //Transient Track Studies
    double closestApproach;
    double staMinDr;
    double standaloneDEoverE;
    double dxy;
    double dsz;
    //HCAL
    int nCellsFound;
    double coneEnergy;
    int coneHits;
    double hitEnergies[7];
    int hitsOverThresh;
    double vtxZ;
    double vtxPerp;
    int bremDepth;
    //Pf Jets
    double nearestJetDr;
    double nearJetMaxE;
    double nearestJetE;
    //Calo Jets
    double minCaloJetDr;
    double caloJetHcalE;
    double caloJetEcalE;
    double caloJetTotalE;
    double caloConeE;
    int caloConeHits;
    //Gen Matching
    double genProbeMatchdE;
    double genProbeMatchdR;
    double probeMuondE;
    double probeMuondR;
    double standaloneE;
    double globalMuonE;
    int particlePDG;
    double minMuDr;
};
#endif
