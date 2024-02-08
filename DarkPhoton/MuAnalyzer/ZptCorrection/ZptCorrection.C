#include"ZptCorrection.h"
ZptCorrection::ZptCorrection(){
}
ZptCorrection::ZptCorrection(TString datapath){
  Setup(datapath);
}
ZptCorrection::~ZptCorrection(){
  Reset();
}
void ZptCorrection::Setup(TString datapath){
  Reset();
  std::ifstream fcheck(datapath);
  if(fcheck.good()){
    std::cout<<"[ZptCorrection::Setup] using file "+datapath<<std::endl;
  }else{
    std::cout<<"[ZptCorrection::Setup] no "+datapath<<std::endl;
    return;
  }
  TFile f(datapath);
  fZptWeightG=(TF1*)f.Get("zptweight_g");
  if(!fZptWeightG){
    std::cout<<"[ZptCorrection::Setup] no zptweight_g"<<std::endl;
    exit(ENODATA);
  }
  fZptWeightYaxis=(TAxis*)f.Get("yaxis");
  if(!fZptWeightYaxis){
    std::cout<<"[ZptCorrection::Setup] no yaxis"<<std::endl;
    exit(ENODATA);
  }
  fZptWeightY.resize(fZptWeightYaxis->GetNbins()+2,NULL);
  for(int i=1;i<fZptWeightYaxis->GetNbins()+1;i++){
    fZptWeightY[i]=(TF1*)f.Get(Form("zptweight_y%d",i));
    if(!fZptWeightY[i]){
      std::cout<<"[ZptCorrection::Setup] no zptweight_y"+TString(i)<<std::endl;
      exit(ENODATA);
    }
  }
  fZptWeightMaxis=(TAxis*)f.Get("maxis");
  if(!fZptWeightMaxis){
    std::cout<<"[ZptCorrection::Setup] no maxis"<<std::endl;
    exit(ENODATA);
  }
  fZptWeightM.resize(fZptWeightMaxis->GetNbins()+2,NULL);
  for(int i=1;i<fZptWeightMaxis->GetNbins()+1;i++){
    fZptWeightM[i]=(TF1*)f.Get(Form("zptweight_m%d",i));
    if(!fZptWeightM[i]){
      std::cout<<"[ZptCorrection::Setup] no zptweight_m"+TString(i)<<std::endl;
      exit(ENODATA);
    }
  }
}
void ZptCorrection::Reset(){
  if(fZptWeightG){
    delete fZptWeightG;
    fZptWeightG=NULL;
  }
  if(fZptWeightYaxis){
    delete fZptWeightYaxis;
    fZptWeightYaxis=NULL;
  }
  for(auto f:fZptWeightY){
    if(f) delete f;
  }
  fZptWeightY.clear();
  if(fZptWeightMaxis){
    delete fZptWeightMaxis;
    fZptWeightMaxis=NULL;
  }
  for(auto f:fZptWeightM){
    if(f) delete f;
  }
  fZptWeightM.clear();
}
double ZptCorrection::GetZptWeight(double pt) const {
  if(!fZptWeightG) return 1.;
  if(pt<0) pt=0;
  if(pt>=650) pt=649.9;
  double sf=fZptWeightG->Eval(pt);
  return sf;
}
double ZptCorrection::GetZptWeight(double pt,double rapidity) const {
  double sf=GetZptWeight(pt);

  if(!fZptWeightY.size()) return sf;
  if(std::isnan(rapidity)) return sf;
  double y=fabs(rapidity);
  if(y>=fZptWeightYaxis->GetXmax()) y=fZptWeightYaxis->GetXmax()-1e-6;  

  double ymin=fZptWeightYaxis->GetBinCenter(1);
  double ymax=fZptWeightYaxis->GetBinCenter(fZptWeightYaxis->GetNbins());
  int biny1,biny2;
  if(y<ymin){
    biny1=1;
    biny2=2;
  }else if(y>=ymax){
    biny1=fZptWeightYaxis->GetNbins()-1;
    biny2=fZptWeightYaxis->GetNbins();
  }else{
    int biny=fZptWeightYaxis->FindBin(y);
    if(y>=fZptWeightYaxis->GetBinCenter(biny)){
      biny1=biny;
      biny2=biny+1;
    }else{
      biny1=biny-1;
      biny2=biny;
    }
  }
  double y1=fZptWeightYaxis->GetBinCenter(biny1);
  double y2=fZptWeightYaxis->GetBinCenter(biny2);
  sf*=( (y2-y)*fZptWeightY[biny1]->Eval(pt) + (y-y1)*fZptWeightY[biny2]->Eval(pt) )/(y2-y1);
  return sf;
}
double ZptCorrection::GetZptWeight(double pt,double rapidity,double mass) const {
  double sf=GetZptWeight(pt,rapidity);
  if(!fZptWeightM.size()) return sf;
  if(mass==0) return sf;
  double m=mass;
  if(m<fZptWeightMaxis->GetXmin()) m=fZptWeightMaxis->GetXmin();
  if(m>=fZptWeightMaxis->GetXmax()) m=fZptWeightMaxis->GetXmax()-1e-6;

  double mmin=fZptWeightMaxis->GetBinCenter(1);
  double mmax=fZptWeightMaxis->GetBinCenter(fZptWeightMaxis->GetNbins());
  int binm1,binm2;
  if(m<mmin){
    binm1=1;
    binm2=2;
  }else if(m>=mmax){
    binm1=fZptWeightMaxis->GetNbins()-1;
    binm2=fZptWeightMaxis->GetNbins();
  }else{
    int binm=fZptWeightMaxis->FindBin(m);
    if(m>=fZptWeightMaxis->GetBinCenter(binm)){
      binm1=binm;
      binm2=binm+1;
    }else{
      binm1=binm-1;
      binm2=binm;
    }
  }
  double m1=fZptWeightMaxis->GetBinCenter(binm1);
  double m2=fZptWeightMaxis->GetBinCenter(binm2);
  sf*=( (m2-m)*fZptWeightM[binm1]->Eval(pt) + (m-m1)*fZptWeightM[binm2]->Eval(pt) )/(m2-m1);
  return sf;
}
