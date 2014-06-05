#include <iostream>

#include "FWCore/Utilities/interface/typelookup.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"

#include "CalibCalorimetry/EcalLaserAnalyzer/interface/MEEBGeom.h"
#include "CalibCalorimetry/EcalLaserAnalyzer/interface/MEEEGeom.h"
// #include "CalibCalorimetry/EcalLaserAnalyzer/interface/ME.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

EcalLaserDbService::EcalLaserDbService () 
  : 
  mAlphas_ (0),
  mAPDPNRatiosRef_ (0),
  mAPDPNRatios_ (0)
 {}



const EcalLaserAlphas* EcalLaserDbService::getAlphas () const {
  return mAlphas_;
}

const EcalLaserAPDPNRatiosRef* EcalLaserDbService::getAPDPNRatiosRef () const {
  return mAPDPNRatiosRef_;
}

const EcalLaserAPDPNRatios* EcalLaserDbService::getAPDPNRatios () const {
  return mAPDPNRatios_;
}


float EcalLaserDbService::getLaserCorrection (DetId const & xid, edm::Timestamp const & iTime) const {
  
  float correctionFactor = 1.0;

  const EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap& laserRatiosMap =  mAPDPNRatios_->getLaserMap();
  const EcalLaserAPDPNRatios::EcalLaserTimeStampMap& laserTimeMap =  mAPDPNRatios_->getTimeMap();
  const EcalLaserAPDPNRatiosRefMap& laserRefMap =  mAPDPNRatiosRef_->getMap();
  const EcalLaserAlphaMap& laserAlphaMap =  mAlphas_->getMap();

  EcalLaserAPDPNRatios::EcalLaserAPDPNpair apdpnpair;
  EcalLaserAPDPNRatios::EcalLaserTimeStamp timestamp;
  EcalLaserAPDPNref apdpnref;
  EcalLaserAlpha alpha;

  if (xid.det()==DetId::Ecal) {
    //    std::cout << " XID is in Ecal : ";
  } else {
    //    std::cout << " XID is NOT in Ecal : ";
    edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL" << endl;
    return correctionFactor;
  } 

//  int hi = -1;
//  if (xid.subdetId()==EcalBarrel) {
//    //    std::cout << "EcalBarrel" << std::endl;
//    //    std::cout << "--> rawId() = " << xid.rawId() << "   id() = " << EBDetId( xid ).hashedIndex() << std::endl;
//    hi = EBDetId( xid ).hashedIndex();
//  } else if (xid.subdetId()==EcalEndcap) {
//    //    std::cout << "EcalEndcap" << std::endl;
//    hi = EEDetId( xid ).hashedIndex() + EBDetId::MAX_HASH + 1;
//
//  } else {
//    //    std::cout << "NOT EcalBarrel or EcalEndCap" << std::endl;
//    edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL Barrel or Endcap" << endl;
//    return correctionFactor;
//  }

  int iLM;
  if (xid.subdetId()==EcalBarrel) {
    EBDetId ebid( xid.rawId() );
    iLM = MEEBGeom::lmr(ebid.ieta(), ebid.iphi());
  } else if (xid.subdetId()==EcalEndcap) {
    EEDetId eeid( xid.rawId() );
    // SuperCrystal coordinates
    MEEEGeom::SuperCrysCoord iX = (eeid.ix()-1)/5 + 1;
    MEEEGeom::SuperCrysCoord iY = (eeid.iy()-1)/5 + 1;    
    iLM = MEEEGeom::lmr(iX, iY, eeid.zside());    
  } else {
    edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL Barrel or Endcap" << endl;
    return correctionFactor;
  }
  //  std::cout << " LM num ====> " << iLM << endl;

  // get alpha, apd/pn ref, apd/pn pairs and timestamps for interpolation

  EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap::const_iterator itratio = laserRatiosMap.find(xid);
  if (itratio != laserRatiosMap.end()) {
    apdpnpair = (*itratio);
  } else {
    edm::LogError("EcalLaserDbService") << "error with laserRatiosMap!" << endl;     
    return correctionFactor;
  }

  if (iLM-1< (int)laserTimeMap.size()) {
    timestamp = laserTimeMap[iLM-1];  
  } else {
    edm::LogError("EcalLaserDbService") << "error with laserTimeMap!" << endl;     
    return correctionFactor;
  }

  EcalLaserAPDPNRatiosRefMap::const_iterator itref = laserRefMap.find(xid);
  if ( itref != laserRefMap.end() ) {
    apdpnref = (*itref);
  } else { 
    edm::LogError("EcalLaserDbService") << "error with laserRefMap!" << endl;     
    return correctionFactor;
  }

  EcalLaserAlphaMap::const_iterator italpha = laserAlphaMap.find(xid);
  if ( italpha != laserAlphaMap.end() ) {
    alpha = (*italpha);
  } else {
    edm::LogError("EcalLaserDbService") << "error with laserAlphaMap!" << endl;     
    return correctionFactor;
  }

  //    std::cout << " APDPN pair " << apdpnpair.p1 << " , " << apdpnpair.p2 << std::endl; 
  //    std::cout << " TIME pair " << timestamp.t1.value() << " , " << timestamp.t2.value() << " iLM " << iLM << std::endl; 
  //    std::cout << " LM module " << iLM << std::endl;
  //    std::cout << " APDPN ref " << apdpnref << std::endl; 
  //    std::cout << " ALPHA " << alpha << std::endl; 
  
  // should implement some default in case of error...

  // should do some quality checks first
  // ...

  // we will need to treat time differently...
  // is time in DB same format as in MC?  probably not...
  
  // interpolation

  edm::TimeValue_t t = iTime.value();
  edm::TimeValue_t t_i = 0, t_f = 0;
  float p_i = 0, p_f = 0;

  if ( t >= timestamp.t1.value() && t < timestamp.t2.value() ) {
          t_i = timestamp.t1.value();
          t_f = timestamp.t2.value();
          p_i = apdpnpair.p1;
          p_f = apdpnpair.p2;
  } else if ( t >= timestamp.t2.value() && t <= timestamp.t3.value() ) {
          t_i = timestamp.t2.value();
          t_f = timestamp.t3.value();
          p_i = apdpnpair.p2;
          p_f = apdpnpair.p3;
  } else if ( t < timestamp.t1.value() ) {
          t_i = timestamp.t1.value();
          t_f = timestamp.t2.value();
          p_i = apdpnpair.p1;
          p_f = apdpnpair.p2;
          //edm::LogWarning("EcalLaserDbService") << "The event timestamp t=" << t 
          //        << " is lower than t1=" << t_i << ". Extrapolating...";
  } else if ( t > timestamp.t3.value() ) {
          t_i = timestamp.t2.value();
          t_f = timestamp.t3.value();
          p_i = apdpnpair.p2;
          p_f = apdpnpair.p3;
          //edm::LogWarning("EcalLaserDbService") << "The event timestamp t=" << t 
          //        << " is greater than t3=" << t_f << ". Extrapolating...";
  }

  if ( apdpnref != 0 && (t_i - t_f) != 0) {
    float interpolatedLaserResponse = p_i/apdpnref + (t-t_i)*(p_f-p_i)/apdpnref/(t_f-t_i);
    if ( interpolatedLaserResponse <= 0 ) {
      edm::LogWarning("EcalLaserDbService") << "The interpolated laser correction is <= zero! (" 
                    << interpolatedLaserResponse << "). Using 1. as correction factor.";
            return correctionFactor;
    } else {
      correctionFactor = 1/pow(interpolatedLaserResponse,alpha);
    }
    
  } else {
    edm::LogError("EcalLaserDbService") 
            << "apdpnref (" << apdpnref << ") "
            << "or t_i-t_f (" << (t_i - t_f) << " is zero!";
    return correctionFactor;
  }
  
  return correctionFactor;
}


TYPELOOKUP_DATA_REG(EcalLaserDbService);
