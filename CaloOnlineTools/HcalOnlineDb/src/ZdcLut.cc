#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/ZdcLut.h"
#include <cstdlib>
#include <map>

/**

   \class ZdcLut
   \brief Generation of ZDC Lookup tables and associate helper methods
   \brief Adopted to CMSSW HCAL LUT manager specs
   \brief by Gena Kukartsev, Brown University, Dec 08, 2009
   \author Elijah Dunn

   \brief Updated by Quan Wang @KU for QIE10 upgrade, 2022
*/

ZdcLut::ZdcLut(const HcalDbService* _conditions)
//ZdcLut::ZdcLut()
    : Empty_LUT(256,0)
{
  int lsb = 50;
  conditions = _conditions;
  // ADC to fC map for QIE10

  std::vector<float> fC_TDR{1.62607, 4.87821, 8.13035, 11.3825, 14.6346, 17.8868, 21.1389, 24.391, 27.6432, 30.8953, 34.1475, 37.3996, 40.6517, 43.9039, 47.156, 50.4081, 55.2864, 61.7906, 68.2949, 74.7992, 81.3035, 87.8077, 94.312, 100.816, 107.321, 113.825, 120.329, 126.833, 133.338, 139.842, 146.346, 152.851, 159.355, 165.859, 172.363, 178.868, 188.624, 201.633, 214.641, 227.65, 240.658, 253.667, 266.675, 279.684, 292.692, 305.701, 318.71, 331.718, 344.727, 357.735, 370.744, 383.752, 396.761, 409.769, 422.778, 435.787, 448.795, 468.308, 494.325, 520.342, 546.359, 572.376, 598.393, 624.411, 568.652, 594.369, 620.087, 645.805, 671.523, 697.24, 722.958, 748.676, 774.393, 800.111, 825.829, 851.546, 877.264, 902.982, 928.699, 954.417, 992.994, 1044.43, 1095.86, 1147.3, 1198.73, 1250.17, 1301.61, 1353.04, 1404.48, 1455.91, 1507.35, 1558.78, 1610.22, 1661.65, 1713.09, 1764.52, 1815.96, 1867.39, 1918.83, 1970.27, 2047.42, 2150.29, 2253.16, 2356.03, 2458.9, 2561.77, 2664.64, 2767.51, 2870.38, 2973.26, 3076.13, 3179, 3281.87, 3384.74, 3487.61, 3590.48, 3693.35, 3796.22, 3899.09, 4001.96, 4104.83, 4259.14, 4464.88, 4670.62, 4876.36, 5082.11, 5287.85, 5493.59, 4996.09, 5200.82, 5405.55, 5610.27, 5815, 6019.73, 6224.46, 6429.19, 6633.91, 6838.64, 7043.37, 7248.1, 7452.83, 7657.55, 7862.28, 8067.01, 8374.1, 8783.56, 9193.01, 9602.47, 10011.9, 10421.4, 10830.8, 11240.3, 11649.7, 12059.2, 12468.7, 12878.1, 13287.6, 13697, 14106.5, 14515.9, 14925.4, 15334.9, 15744.3, 16153.8, 16767.9, 17586.9, 18405.8, 19224.7, 20043.6, 20862.5, 21681.4, 22500.3, 23319.2, 24138.2, 24957.1, 25776, 26594.9, 27413.8, 28232.7, 29051.6, 29870.5, 30689.4, 31508.4, 32327.3, 33146.2, 34374.6, 36012.4, 37650.2, 39288, 40925.8, 42563.7, 44201.5, 41124.8, 42733.7, 44342.5, 45951.4, 47560.2, 49169.1, 50777.9, 52386.8, 53995.6, 55604.4, 57213.3, 58822.1, 60431, 62039.8, 63648.7, 65257.5, 67670.8, 70888.5, 74106.2, 77323.9, 80541.6, 83759.3, 86977, 90194.7, 93412.4, 96630.1, 99847.8, 103065, 106283, 109501, 112719, 115936, 119154, 122372, 125589, 128807, 133634, 140069, 146504, 152940, 159375, 165811, 172246, 178681, 185117, 191552, 197988, 204423, 210858, 217294, 223729, 230165, 236600, 243035, 249471, 255906, 262342, 271995, 284865, 297736, 310607, 323478, 336349, 349219};

  // hard coded ZDC calibration from 2018 PbPb data taking
  std::map<std::string, float> chmap{
    {"ZDCP_EM1",  0.11},
    {"ZDCP_EM2",  0.11},
    {"ZDCP_EM3",  0.08},
    {"ZDCP_EM4",  0.10},
    {"ZDCP_EM5",  0.06},
    {"ZDCP_HAD1", 1.00},
    {"ZDCP_HAD2", 0.61},
    {"ZDCP_HAD3", 0.46},
    {"ZDCP_HAD4", 0.52},

    {"ZDCM_EM1",  0.12},
    {"ZDCM_EM2",  0.12},
    {"ZDCM_EM3",  0.08},
    {"ZDCM_EM4",  0.10},
    {"ZDCM_EM5",  0.11},
    {"ZDCM_HAD1", 1.00},
    {"ZDCM_HAD2", 0.81},
    {"ZDCM_HAD3", 0.46},
    {"ZDCM_HAD4", 0.58},

    {"Pscale", 2.67},
    {"Mscale", 4.45}
  };

/*  const HcalQIECoder* channelCoder = _conditions.getHcalCoder(cell);
  const HcalQIEShape* shape = _conditions.getHcalShape(cell);
  auto adc2fC = [channelCoder, shape](unsigned int adc) {
    float fC = 0;
    for (auto capId : {0, 1, 2, 3})
      fC += channelCoder->charge(*shape, adc, capId);
    return fC / 4;
  };*/


  auto zdc_lut_calib = [&fC_TDR, &lsb](float calib, float scale) {
    std::vector<int> ret(256);
    for ( int i = 0; i < 256; i++ ) {
	//ret[i] = int(fC_TDR[i]);
//MICHAEL: Uncomment this line for true behavior. Otherwise keep above uncommented, which gives a
// 1-to-1 LUT for testing purposes.
        //ret[i] = int(fC_TDR[i] * calib /scale /lsb + 0.5);
        ret[i] = int(i);
    }
    return ret;
  };

  ZDCP_EM_LUT[0] = zdc_lut_calib(chmap["ZDCP_EM1"], chmap["Pscale"]);
  ZDCP_EM_LUT[1] = zdc_lut_calib(chmap["ZDCP_EM2"], chmap["Pscale"]);
  ZDCP_EM_LUT[2] = zdc_lut_calib(chmap["ZDCP_EM3"], chmap["Pscale"]);
  ZDCP_EM_LUT[3] = zdc_lut_calib(chmap["ZDCP_EM4"], chmap["Pscale"]);
  ZDCP_EM_LUT[4] = zdc_lut_calib(chmap["ZDCP_EM5"], chmap["Pscale"]);
  ZDCM_EM_LUT[0] = zdc_lut_calib(chmap["ZDCM_EM1"], chmap["Mscale"]);
  ZDCM_EM_LUT[1] = zdc_lut_calib(chmap["ZDCM_EM2"], chmap["Mscale"]);
  ZDCM_EM_LUT[2] = zdc_lut_calib(chmap["ZDCM_EM3"], chmap["Mscale"]);
  ZDCM_EM_LUT[3] = zdc_lut_calib(chmap["ZDCM_EM4"], chmap["Mscale"]);
  ZDCM_EM_LUT[4] = zdc_lut_calib(chmap["ZDCM_EM5"], chmap["Mscale"]);

  ZDCP_HAD_LUT[0] = zdc_lut_calib(chmap["ZDCP_HAD1"], chmap["Pscale"]);
  ZDCP_HAD_LUT[1] = zdc_lut_calib(chmap["ZDCP_HAD2"], chmap["Pscale"]);
  ZDCP_HAD_LUT[2] = zdc_lut_calib(chmap["ZDCP_HAD3"], chmap["Pscale"]);
  ZDCP_HAD_LUT[3] = zdc_lut_calib(chmap["ZDCP_HAD4"], chmap["Pscale"]);
  ZDCM_HAD_LUT[0] = zdc_lut_calib(chmap["ZDCM_HAD1"], chmap["Mscale"]);
  ZDCM_HAD_LUT[1] = zdc_lut_calib(chmap["ZDCM_HAD2"], chmap["Mscale"]);
  ZDCM_HAD_LUT[2] = zdc_lut_calib(chmap["ZDCM_HAD3"], chmap["Mscale"]);
  ZDCM_HAD_LUT[3] = zdc_lut_calib(chmap["ZDCM_HAD4"], chmap["Mscale"]);
}

ZdcLut::~ZdcLut(void) {}

// get LUT by proper ZDC channel
std::vector<int> ZdcLut::get_lut(std::string zdc_section, int zdc_side, int zdc_channel) {

  std::cout << "INSIDE get_lut" << std::endl;

  HcalZDCDetId::Section section = HcalZDCDetId::Unknown;
  if (zdc_section == "ZDC EM")
    section = HcalZDCDetId::EM;
  else if (zdc_section == "ZDC HAD")
    section = HcalZDCDetId::HAD;
  else if (zdc_section == "ZDC LUM")
    section = HcalZDCDetId::LUM;
  else if (zdc_section == "ZDC RPD")
    section = HcalZDCDetId::RPD;
  HcalZDCDetId _detid(section, (zdc_side>0), zdc_channel);
  
  if (zdc_section.find("ZDC RPD") != std::string::npos) {
    return std::vector<int>();
  }else if (zdc_section.find("ZDC EM") != std::string::npos) {
    if ( (zdc_channel>5) or (zdc_channel<1) ) return std::vector<int>();
  }else if (zdc_section.find("ZDC HAD") != std::string::npos) {
    if ( (zdc_channel>4) or (zdc_channel<1) ) return std::vector<int>();
  }

  const HcalLutMetadata* metadata = conditions->getHcalLutMetadata();
  assert(metadata != nullptr);
  float nominalgain_ = metadata->getNominalGain();

  std::cout << "nominalgain_: " << nominalgain_ << std::endl;
    
  const HcalQIECoder* channelCoder = conditions->getHcalCoder(_detid);
  const HcalQIEShape* shape = conditions->getHcalShape(_detid);

  auto adc2fC = [channelCoder, shape](unsigned int adc) {
    float fC = 0;
    for (auto capId : {0, 1, 2, 3}){
      fC += channelCoder->charge(*shape, adc, capId);
    }
    return fC / 4;
  };

  float gain = 0;
  float ped = 0;
  const HcalCalibrations& calibrations = conditions->getHcalCalibrations(_detid);
  for (auto capId : {0, 1, 2, 3}) {
    ped += calibrations.effpedestal(capId);
    gain += calibrations.LUTrespcorrgain(capId);
  }
  ped /= 4.0;
  gain /= 4.0;

  std::cout << "ped: " << ped << std::endl;
  std::cout << "gain: " << gain << std::endl;

  /*const HcalLutMetadatum* meta = metadata->getValues(_detid);
  std::cout << "check 1" << std::endl;
  float rcalib = meta->getRCalib();  
  std::cout << "rcalib: " << rcalib << std::endl;*/

  for (unsigned int adc = 0; adc < 256; ++adc) {
     std::cout << "adc: " << adc << std::endl;
     std::cout << "adc2fC(adc): " << adc2fC(adc) << std::endl;
  }
  
  if (zdc_section.find("ZDC EM") != std::string::npos) {
      if ( (zdc_channel>5) or (zdc_channel<1) ) return std::vector<int>();
      if ( zdc_side == -1 ) {
          return ZDCM_EM_LUT[zdc_channel-1];
      }
      if ( zdc_side == 1 ) {
          return ZDCP_EM_LUT[zdc_channel-1];
      }
  } else if (zdc_section.find("ZDC HAD") != std::string::npos) {
      if ( (zdc_channel>4) or (zdc_channel<1) ) return std::vector<int>();
      if ( zdc_side == -1 ) {
          return ZDCM_HAD_LUT[zdc_channel-1];
      }
      if ( zdc_side == 1 ) {
          return ZDCP_HAD_LUT[zdc_channel-1];
      }
  } else {
    std::cout << zdc_section << ": unknown ZDC section, exiting... No, continue..." << std::endl;
    //exit(-1);
  }
  return std::vector<int>();
}
