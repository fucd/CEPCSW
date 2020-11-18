/***********************************************************************************\
* (c) Copyright 1998-2019 CERN for the benefit of the LHCb and ATLAS collaborations *
*                                                                                   *
* This software is distributed under the terms of the Apache version 2 licence,     *
* copied verbatim in the file "LICENSE".                                            *
*                                                                                   *
* In applying this licence, CERN does not waive the privileges and immunities       *
* granted to it by virtue of its status as an Intergovernmental Organization        *
* or submit itself to any jurisdiction.                                             *
\***********************************************************************************/
// Include files
#include "DumpHitAlg.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "DataHelper/HelixClass.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include <math.h>

DECLARE_COMPONENT( DumpHitAlg )

//------------------------------------------------------------------------------
DumpHitAlg::DumpHitAlg( const std::string& name, ISvcLocator* pSvcLocator )
    : Algorithm( name, pSvcLocator ) {
  declareProperty("SimTrackerHitCollection", _inSimTrackerHitColHdl, "Handle of the Input Hit collection from CEPCSW");
  declareProperty("TrackerHitCollection", _inTrackerHitColHdl, "Handle of the Input Hit collection from CEPCSW");

  m_thisName = name;
}

//------------------------------------------------------------------------------
StatusCode DumpHitAlg::initialize(){
  info() << "Booking Ntuple" << endmsg;

  NTuplePtr nt1(ntupleSvc(), "MyTuples/Track"+m_thisName);
  if ( !nt1 ) {
    m_tuple = ntupleSvc()->book("MyTuples/Track"+m_thisName,CLID_ColumnWiseTuple,"Tracking result");
    if ( 0 != m_tuple ) {
      m_tuple->addItem ("nhit",      m_nHits, 0, 5000 ).ignore();
      m_tuple->addIndexedItem ("id",        m_nHits, m_cellID ).ignore();
      m_tuple->addIndexedItem ("x",         m_nHits, m_x ).ignore();
      m_tuple->addIndexedItem ("y",         m_nHits, m_y ).ignore();
      m_tuple->addIndexedItem ("z",         m_nHits, m_z ).ignore();
    }
    else { // did not manage to book the N tuple....
      fatal() << "Cannot book MyTuples/Track"+m_thisName <<endmsg; 
      return StatusCode::FAILURE;
    }
  }
  else{
    m_tuple = nt1;
  }

  _nEvt = 0;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode DumpHitAlg::execute(){
  const edm4hep::SimTrackerHitCollection* hitCol = nullptr;
  try {
    hitCol = _inSimTrackerHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inTrackerHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  if(hitCol){
    m_nHits = 0;
    for(auto hit : *hitCol){
      m_cellID[m_nHits] = hit.getCellID();
      const edm4hep::Vector3d& pos = hit.getPosition();
      m_x[m_nHits] = pos.x;
      m_y[m_nHits] = pos.y;
      m_z[m_nHits] = pos.z;
      m_nHits++;
    }
  }
  m_tuple->write();
  _nEvt++;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode DumpHitAlg::finalize(){
  debug() << "Finalizing..." << endmsg;

  return StatusCode::SUCCESS;
}
