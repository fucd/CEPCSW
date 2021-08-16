#ifndef TDTHit_h
#define TDTHit_h

#include "TVector3.h"

class TDTHit {
 public:
  TDTHit();
  TDTHit(const TDTHit& h);
  TDTHit(const TVector3& start, const TVector3& end, double radiation, bool isSens, double resRPhi=0, double resZ=0, bool isPixel=false, bool isGas=false);  
  TDTHit(double x1, double y1, double z1, double x2, double y2, double z2, double radiation, bool isSens, double resRPhi=0, double resZ=0, bool isPixel=false, bool isGas=false);

  ~TDTHit(){};

  void SetStartPoint(const TVector3& start)        { m_startPoint = start; };
  void SetStartPoint(double x, double y, double z) { m_startPoint = TVector3(x,y,z); };
  void SetEndPoint(const TVector3& end)            { m_endPoint = end; };;
  void SetEndPoint(double x, double y, double z)   { m_endPoint = TVector3(x,y,z); };
  void SetRadiation(double radiation)              { m_XtoX0 = radiation; };  
  void SetResolutionRPhi(double res)               { m_resolutionRPhi = res; };
  void SetResolutionZ(double res)                  { m_resolutionZ = res; };
  void SetPhi(double phi)                          { m_phi = phi; };
  void SetDeltaPhi(double delta)                   { m_deltaPhi = delta; }; 
  void SetPathLength(double length)                { m_pathLength = length; };  

  TVector3 GetStartPoint() const { return m_startPoint; };
  TVector3 GetEndPoint()   const { return m_endPoint; };
  //double   GetRPos()       const { return (m_startPoint+m_endPoint).Perp()/2.; };
  //double   GetZPos()       const { return (m_startPoint.z()+m_endPoint.z())/2.; };
  double   GetRPos()       const { return m_endPoint.Perp(); };
  double   GetZPos()       const { return m_endPoint.z(); };

  double   GetRadiation()      const { return m_XtoX0; };
  double   GetResolutionRPhi() const { return m_resolutionRPhi; };
  double   GetResolutionZ()    const { return m_resolutionZ; };
  double   GetPhi()            const { return m_phi; }; 
  double   GetDeltaPhi()       const { return m_deltaPhi; }; 
  double   GetPathLength()     const { return m_pathLength; };

  bool     isSensitive() const { return m_isSens; };
  bool     isPixel()     const { return m_isPixel; };
  bool     isGas()       const { return m_isGas; };
  bool     isPassive()   const { return !m_isSens; };
  
 private:
  TVector3 m_startPoint;
  TVector3 m_endPoint;
  bool     m_isSens;
  bool     m_isPixel;
  bool     m_isGas;

  double   m_XtoX0;
  double   m_resolutionRPhi;
  double   m_resolutionZ;
  double   m_phi;
  double   m_deltaPhi;
  double   m_pathLength;
};
typedef std::unique_ptr<TDTHit> TDTHitPtr;
typedef std::vector<TDTHitPtr>  TDTHitVec;
#endif
  
