#ifndef ZdcLut_h
#define ZdcLut_h

/**

   \class ZdcLut
   \brief Generation of ZDC Lookup tables and associate helper methods
   \brief Adopted to CMSSW HCAL LUT manager specs
   \brief by Gena Kukartsev, Brown University, Dec 08, 2009
   \author Elijah Dunn

   \brief Updated by Quan Wang for QIE10
*/

#include <iostream>
#include <vector>
#include <string>

class ZdcLut {
public:
  ZdcLut(const HcalDbService* _conditions);
  ~ZdcLut();

  //get_lut returns a specific lut based on side, fiber, and fiber_channel
  //std::vector <int> get_lut(int side_num, int fiber_num, int channel_num){ return side[side_num].fiber[fiber_num].channel[channel_num].LUT; }
  //std::vector<int> get_lut(int emap_side, int emap_htr_fiber, int emap_fi_ch);

  std::vector<int> get_lut(std::string zdc_section, int zdc_side, int zdc_channel);

  //int simple_loop(void);

private:
  //variable
  std::vector<int> ZDCP_EM_LUT[5];
  std::vector<int> ZDCM_EM_LUT[5];
  std::vector<int> ZDCP_HAD_LUT[4];
  std::vector<int> ZDCM_HAD_LUT[4];
  std::vector<int> Empty_LUT;

  const HcalDbService* conditions;
};

#endif
