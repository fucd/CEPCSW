#ifndef CompareTrackAlg_h
#define CompareTrackAlg_h 1

// Include files
#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"

#include "AIDA/IAxis.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IProfile1D.h"
#include "AIDA/IProfile2D.h"
#include "GaudiKernel/NTuple.h"

using namespace AIDA;
// Forward declarations
class HistogramSvc;

class CompareTrackAlg : public Algorithm {

public:
  // Constructor of this form must be provided
  CompareTrackAlg( const std::string& name, ISvcLocator* pSvcLocator );

  // Three mandatory member functions of any algorithm
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;

private:
  DataHandle<edm4hep::MCParticleCollection> _inMCColHdl{"MCParticle", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackCollection> _inTrackColHdl{"SiTracks", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackCollection> _inTrackMarlinColHdl{"SiTracks_Marlin", Gaudi::DataHandle::Reader, this};

  Gaudi::Property<double> m_field{this, "Field", 3.0};
  Gaudi::Property<double> m_cut{this, "SigmaCut", 5.0};

  IHistogram1D* m_hPtMC;
  IHistogram1D* m_hThetaMC;
  IHistogram1D* m_hPt[2];
  IHistogram1D* m_hEffPt[2];
  IHistogram1D* m_hTheta[2];
  IHistogram1D* m_hEffTheta[2];
  IHistogram1D* m_hPullD0[2];
  IHistogram1D* m_hPullPhi0[2];
  IHistogram1D* m_hPullOmega[2];
  IHistogram1D* m_hPullZ0[2];
  IHistogram1D* m_hPullTanLambda[2];
  IHistogram1D* m_hDiffPR[2];
  IHistogram1D* m_hDiffInvPt[2];
  IHistogram1D* m_hDiffD0[2];
  NTuple::Tuple*       m_tuple;
  NTuple::Array<float> m_dInvPt;
  NTuple::Array<float> m_dD0;
  
  int _nEvt;
  int m_nTrack[2];
  int m_nMatched[2];
  int m_nFilled[2];
  std::string m_thisName;
};

#endif // HISTOGRAMS_HISTOALGORITHM_H
