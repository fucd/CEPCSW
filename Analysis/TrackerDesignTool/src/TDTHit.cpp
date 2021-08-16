#include "TDTHit.h"

TDTHit::TDTHit(){
  m_startPoint     = TVector3(0,0,0);
  m_endPoint       = TVector3(0,0,0);
  m_isSens         = false;
  m_isPixel        = false;
  m_isGas          = false;
  m_XtoX0          = 0;
  m_resolutionRPhi = 0;
  m_resolutionZ    = 0;
}

TDTHit::TDTHit(const TDTHit& h){
  m_startPoint     = h.GetStartPoint();
  m_endPoint       = h.GetEndPoint();
  m_isSens         = h.isSensitive();
  m_isPixel        = h.isPixel();
  m_isGas          = h.isGas();
  m_XtoX0          = h.GetRadiation();
  m_resolutionRPhi = h.GetResolutionRPhi();
  m_resolutionZ    = h.GetResolutionZ();
}

TDTHit::TDTHit(const TVector3& start, const TVector3& end, double radiation, bool isSens, double resRPhi, double resZ, bool isPixel, bool isGas){
  m_startPoint     = start;
  m_endPoint       = end;
  m_isSens         = isSens;
  m_isPixel        = isPixel;
  m_isGas          = isGas;
  m_XtoX0          = radiation;
  m_resolutionRPhi = resRPhi;
  m_resolutionZ    = resZ;
}

TDTHit::TDTHit(double x1, double y1, double z1, double x2, double y2, double z2, double radiation, bool isSens, double resRPhi, double resZ, bool isPixel, bool isGas){
  m_startPoint     = TVector3(x1,y1,z1);
  m_endPoint       = TVector3(x2,y2,z2);
  m_isSens         = isSens;
  m_isPixel        = isPixel;
  m_isGas          = isGas;
  m_XtoX0          = radiation;
  m_resolutionRPhi = resRPhi;
  m_resolutionZ    = resZ;
}

