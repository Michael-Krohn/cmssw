# Group 0 LUT: bits 12-15 for LLP fine grain bits
Bits 12-15 of the uHTR Group 0 LUT are used to set the energy thresholds for fine grain bits. This is described in Section 3.5 of the [Specifications for HCAL uHTR Firmware](https://cms-docdb.cern.ch/cgi-bin/PublicDocDB/RetrieveFile?docid=12306&filename=LLPbits_uhtr_spec.pdf&version=22). There is not yet infrastructure in HE to support these high bits to flatten the TDC space.

In central CMSSW, from [PR #35599](https://github.com/cms-sw/cmssw/pull/35599) the high bits for Group 0 are set in both HB and HE. A single file modification is needed to only set the bits in HB, to agree with current firmware. File: HcaluLUTTPGCoder.cc. Group 0 LUTs have been successfully created and loaded at P5 with these modifications already.

Done in CMSSW_12_2_1.