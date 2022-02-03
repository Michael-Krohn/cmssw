#include "FWCore/Framework/interface/one/EDAnalyzer.h"
  
#include "DarkPhoton/MuAnalyzer/interface/SigEventInfo.h"
#include <iostream>

SigEventInfo::SigEventInfo()
{
   eventWeight = 1;
   dbremZ = -1;
   dbremPerp = -1;
   sigCat = -1;
   bremDepth = -1;
   dbremInitialE = -1;
   dbremFinalE = -1;
   deflectingAngle = 20;
   nTagMuons = -1;
   nPairedProbes = -1;
   probeTrackIso = -1; 
   probeEcalIso = -1;
   isLoose=false;
   isTight=false;
   isHighPurity=false;
   probeReducedChi=-1;
   tagMuonEta=0;
   tagMuonPt=-1;
   probeTrackEta=0;
   probeTrackPhi=0;
   probeTrackPt=-1;
   pairVtxChi=-1;
   diMuonMass=-1;
   minCSCImpactParameter=-1;
   cscMinDR=-1;
   cscChiSq=-1;
   cscHitAngleChange=100;
   closestApproach=-1;
   standaloneE=-1;
   staMinDr=-1;
   nCellsFound=-1;
   coneEnergy=-1;
   for(int i=0; i<7; i++){hitEnergies[i]=-1;}
   hitsOverThresh=-1;
   vtxZ=-1;  
   bremDepth = -1;
   nearestJetDr=-1;
   nearJetMaxE=-1;
   nearestJetE=-1;
   minCaloJetDr=-1;
   caloJetHcalE=-1;
   caloJetEcalE=-1;
   caloJetTotalE=-1;
   caloConeE=-1;
   caloConeHits=-1;
   hitsOverThresh=-1;
   genProbeMatchdE=-1;
   genProbeMatchdR=-1;
   probeMuondE=-1;
   globalMuonE=-1; 
   nEle=0;
   nExtraEle=0;
   nExtraMuon=0;
   staNHits=0;
   staChi2=-1;
   staTransDCA=-1;
   staLongDCA=-1;
   finalMuE=-1.;
   disStaE=-1.;
   disStaDEoverE=10.;
   nPUmean=-1.;
   pileupWeight=-1.;
}

