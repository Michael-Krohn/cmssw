#ifndef ZptCorrection_h
#define ZptCorrection_h

#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TString.h"
#include "TAxis.h"
#include "TF1.h"

class ZptCorrection {
 public:
  ZptCorrection();
  ZptCorrection(TString datapath);
  ~ZptCorrection();
  void Setup(TString datapath);
  void Reset();
  double GetZptWeight(double pt) const;
  double GetZptWeight(double pt,double rapidity) const;
  double GetZptWeight(double pt,double rapidity,double mass) const;

 private:
  TF1* fZptWeightG=NULL;
  std::vector<TF1*> fZptWeightY;
  TAxis* fZptWeightYaxis=NULL;
  std::vector<TF1*> fZptWeightM;
  TAxis* fZptWeightMaxis=NULL; 
};

#endif
