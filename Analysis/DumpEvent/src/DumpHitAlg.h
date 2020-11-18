#ifndef DumpHitAlg_h
#define DumpHitAlg_h 1

// Include files
#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"

#include "GaudiKernel/NTuple.h"

class DumpHitAlg : public Algorithm {
 public:
  // Constructor of this form must be provided
  DumpHitAlg( const std::string& name, ISvcLocator* pSvcLocator );

  // Three mandatory member functions of any algorithm
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;

 private:
  DataHandle<edm4hep::SimTrackerHitCollection> _inSimTrackerHitColHdl{"FTDCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inTrackerHitColHdl{"FTDCollection", Gaudi::DataHandle::Reader, this};

  NTuple::Tuple*       m_tuple;
  NTuple::Item<long>   m_nHits;
  NTuple::Array<long>  m_cellID;
  NTuple::Array<float> m_x;
  NTuple::Array<float> m_y;
  NTuple::Array<float> m_z;
  
  int _nEvt;
  std::string m_thisName;
};

#endif // HISTOGRAMS_HISTOALGORITHM_H
