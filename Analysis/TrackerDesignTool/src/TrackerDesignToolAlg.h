#ifndef TrackerDesignToolAlg_h
#define TrackerDesignToolAlg_h

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

//#include "edm4hep/TrackCollection.h"
//#include "edm4hep/TrackerHit.h"
//#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleCollection.h"
//#include "edm4hep/SimTrackerHitCollection.h"
//#include "edm4hep/TrackerHitCollection.h"
//#include "edm4hep/MCRecoTrackerAssociationCollection.h"

#include "GaudiKernel/NTuple.h"
#include "TMatrixT.h"

#include "TDTHit.h"
#include "DD4hep/Detector.h"

class TrackerDesignToolAlg : public Algorithm {
 public :
  TrackerDesignToolAlg(const std::string& name, ISvcLocator* pSvcLocator);
  ~TrackerDesignToolAlg(){};
  
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;
  
 private :
  TDTHitVec        findHits(dd4hep::Detector* description, double pt, double theta, double phi=0);
  TMatrixT<double> computeCovarianceMatrixRPhi(TDTHitVec& hits, double pt, double theta, bool msOn=true);
  TMatrixT<double> computeCovarianceMatrixRZ(TDTHitVec& hits, double pt, double theta, bool msOn=true);

  DataHandle<edm4hep::MCParticleCollection> m_mcColHdl{"MCParticle", Gaudi::DataHandle::Reader, this};

  Gaudi::Property<float> m_magField{this, "MagneticField", 3};
  Gaudi::Property<float> m_minStep{this, "MinStep", 1e-5};
  Gaudi::Property<std::vector<std::string> > m_mapSens{this, "Sensitives", {"SiActiveLayer"} };
  Gaudi::Property<std::vector<float> >       m_mapResU{this, "ResolutionU", {0.0040} };
  Gaudi::Property<std::vector<float> >       m_mapResV{this, "ResolutionV", {0.0040} };
  Gaudi::Property<bool> _useTPC{this, "useTPC", true};
  
  NTuple::Tuple* m_tuple;

  //NTuple::Item<long>    m_nParticles;
  //NTuple::Array<double> vx;
};
#endif
