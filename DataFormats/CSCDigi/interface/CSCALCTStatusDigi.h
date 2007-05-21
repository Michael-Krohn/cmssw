#ifndef CSCALCTStatusDigi_CSCALCTStatusDigi_h
#define CSCALCTStatusDigi_CSCALCTStatusDigi_h

/** \class CSCALCTStatusDigi
 *
 *  Digi for CSC ALCT info available in DDU
 *  
 *  $Date: 2007/05/18 18:52:02 $
 *  $Revision: 1.1 $
 *
 */

#include <vector>
#include <boost/cstdint.hpp>

class CSCALCTStatusDigi{

public:

  /// Constructor for all variables 
  CSCALCTStatusDigi (const uint16_t * header, const uint16_t * trailer );

  /// Default constructor.
  CSCALCTStatusDigi () {}

private:

  uint16_t header_[8];
  uint16_t trailer_[4];
};

#include<iostream>
/// needed by COBRA
inline std::ostream & operator<<(std::ostream & o, const CSCALCTStatusDigi& digi) {
  o << " "; 
  o <<"\n";
 
  return o;
}

#endif
