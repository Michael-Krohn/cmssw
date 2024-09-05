#ifndef metFilterTree_h
#define metFilterTree_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "TTree.h"

class metFilterTree {
public:
  metFilterTree();

  void book(TFileDirectory histFolder);
  void FillTree(int run, int lumi, int event, std::string path, bool foundTrig);

public:
  //Event numbers
  int m_event;
  int m_lumi;
  int m_run;
  std::string m_path;
  bool m_foundTrig;
  TTree* failingEvents;
};
#endif
