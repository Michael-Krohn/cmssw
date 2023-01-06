#ifndef SimHistograms_h
#define SimHistograms_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DarkPhoton/MuAnalyzer/interface/SigEventInfo.h"

class TTree;
class TH1F;
class TH2F;

class SimHistograms {
public:
  SimHistograms();

  void book(TFileDirectory histFolder);
  void FillHists(SigEventInfo info);
public:
  //Histograms
  TH1F* m_genEventCount;
  TH1F* m_initialMuE;
  TH1F* m_finalMuE;
  TH1F* m_initialMuEt;
};
#endif
