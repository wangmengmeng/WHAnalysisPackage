#include "SHarper/HEEPAnalyzer/interface/HEEPEleTypeCodes.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

heep::ComCodes heep::EleTypeCodes::codes_ = EleTypeCodes::setCodes_();

heep::ComCodes heep::EleTypeCodes::setCodes_()
{
  heep::ComCodes codes;
  codes.setCode("barrel",BARREL);
  codes.setCode("endcap",ENDCAP);
  codes.setCode("golden",GOLDEN);
  codes.setCode("narrow",NARROW);
  codes.setCode("bigBrem",BIGBREM);
  codes.setCode("showering",SHOWERING);
  codes.setCode("crack",CRACK);
  return codes;
}

int heep::EleTypeCodes::makeTypeCode(int eleType)
{
  int typeCode=0x0;
  if(eleType<100) typeCode |=BARREL;
  else{
    typeCode |=ENDCAP;
    eleType-=100;
  }
  if(eleType/10==0) typeCode |=GOLDEN;
  else if(eleType/10==1) typeCode |= NARROW;
  else if(eleType/10==2) typeCode |= BIGBREM;
  else if(eleType/10==3) typeCode |= SHOWERING;
  else if(eleType/10==4) typeCode |= CRACK;
  else edm::LogWarning("heep::EleTypeCodes") <<"EleTypeCode::makeTypeCode Warning :ele type "<<eleType<<" unrecognised ";
  return typeCode;
}
