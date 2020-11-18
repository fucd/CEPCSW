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
#include "CompareTrackAlg.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "DataHelper/HelixClass.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include <math.h>

DECLARE_COMPONENT( CompareTrackAlg )

//------------------------------------------------------------------------------
CompareTrackAlg::CompareTrackAlg( const std::string& name, ISvcLocator* pSvcLocator )
    : Algorithm( name, pSvcLocator ) {
  declareProperty("MCParticleCollection", _inMCColHdl, "Handle of the Input MCParticle collection");
  declareProperty("CEPCSWTrackCollection", _inTrackColHdl, "Handle of the Input Track collection from CEPCSW");
  declareProperty("MarlinTrackCollection", _inTrackMarlinColHdl, "Handle of the Input Track collection from Marlin");

  m_hPtMC = 0;
  m_hThetaMC = 0;
  for(int i=0;i<2;i++){
    m_hPt[i] = 0;
    m_hEffPt[i] = 0;
    m_hTheta[i] = 0;
    m_hEffTheta[i] = 0;
    m_hPullD0[i] = 0;
    m_hPullPhi0[i] = 0;
    m_hPullOmega[i] = 0;
    m_hPullZ0[i] = 0;
    m_hPullTanLambda[i] = 0;
  }

  m_thisName = name;
}

//------------------------------------------------------------------------------
StatusCode CompareTrackAlg::initialize(){
  info() << "Booking Histograms" << endmsg;

  // Book 1D histogram with fixed and variable binning
  m_hPtMC    = histoSvc()->book( "PtMC_"+m_thisName, "1D fix binning", 100, 0., 100. );
  m_hThetaMC = histoSvc()->book( "ThetaMC_"+m_thisName, "1D fix binning", 100, 0., CLHEP::pi );
  std::string type[2] = {"CEPCSW","Marlin"};
  for(int i=0;i<2;i++){
    std::string suffix = type[i] + m_thisName;
    debug() << "suffix = " << suffix <<endmsg;
    m_hPt[i]    = histoSvc()->book( "Pt_"+suffix, "1D fix binning", 100, 0., 100. );
    m_hTheta[i] = histoSvc()->book( "Theta_"+suffix, "1D fix binning", 100, 0., CLHEP::pi );
    m_hEffPt[i]    = histoSvc()->book( "effPt_"+suffix, "1D fix binning", 100, 0., 100. );
    m_hEffTheta[i] = histoSvc()->book( "effTheta_"+suffix, "1D fix binning", 100, 0., CLHEP::pi );
    
    m_hDiffPR[i] = histoSvc()->book( "diffPR_"+suffix, "1D fix binning", 80, -0.04, 0.04); 
    
    m_hPullD0[i]   = histoSvc()->book( "pullD0_"+suffix, "1D fix binning", 100, -5, 5);
    m_hPullPhi0[i] = histoSvc()->book( "pullPhi0_"+suffix, "1D fix binning", 100, -5, 5);
    m_hPullOmega[i] = histoSvc()->book( "pullOmega_"+suffix, "1D fix binning", 100, -5, 5);
    m_hPullZ0[i]   = histoSvc()->book( "pullZ0_"+suffix, "1D fix binning", 100, -5, 5);
    m_hPullTanLambda[i]   = histoSvc()->book( "pullTanLambda_"+suffix, "1D fix binning", 100, -5, 5);

    m_hDiffInvPt[i]    = histoSvc()->book( "InvPt_"+suffix, "1D fix binning", 2000, -0.1, 0.1 );
    m_hDiffD0[i] = histoSvc()->book( "D0_"+suffix, "1D fix binning", 1000, -0.5, 0.5 ); 

    if ( 0 == m_hPt[i] || 0 == m_hTheta[i] || 0 == m_hEffPt[i] || 0 == m_hEffTheta[i] ||
	 0 == m_hPullD0[i] || 0 == m_hPullPhi0[i] || 0 == m_hPullOmega[i] || 0 == m_hPullZ0[i] || 0 == m_hPullTanLambda[i] ||
	 0 == m_hDiffPR[i] || 0 == m_hDiffInvPt[i] ) {
      error() << "----- Cannot book or register histograms -----" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_nTrack[0] = m_nTrack[1] = 0;
  m_nMatched[0] = m_nMatched[1] = 0;
  m_nFilled[0] = m_nFilled[1] = 0;
  info() << "Finished booking Histograms" << endmsg;

  NTuplePtr nt1(ntupleSvc(), "MyTuples/Tracking"+m_thisName);
  if ( !nt1 ) {
    m_tuple = ntupleSvc()->book("MyTuples/Tracking"+m_thisName,CLID_ColumnWiseTuple,"Tracking result");
    if ( 0 != m_tuple ) {
      m_tuple->addItem ("dInvPt", (long)2, m_dInvPt ).ignore();
      m_tuple->addItem ("dD0",    (long)2, m_dD0 ).ignore();
    }
    else { // did not manage to book the N tuple....
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
StatusCode CompareTrackAlg::execute(){
  const edm4hep::MCParticleCollection* mcCol = nullptr;
  try{
    mcCol = _inMCColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inMCColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  const edm4hep::TrackCollection* trackCols[2] = {nullptr, nullptr};
  try {
    auto trackCol = _inTrackColHdl.get();
    trackCols[0] = trackCol;
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inTrackColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }
  try {
    auto trackCol = _inTrackMarlinColHdl.get();
    trackCols[1] = trackCol;
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inTrackMarlinColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  std::vector<edm4hep::ConstMCParticle> particles;
  std::map<edm4hep::ConstMCParticle, int> map_fit_flag[2];
  particles.reserve(1000);
  if(mcCol){
    for(auto particle : *mcCol){
      edm4hep::Vector3d pos = particle.getVertex();
      edm4hep::Vector3f mom = particle.getMomentum();
      if(pos[0]!=0||pos[1]!=0||pos[2]!=0) continue;
      particles.push_back(particle);
      float pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      m_hPtMC->fill(pt);
      float theta = acos(mom[2]/sqrt(pt*pt+mom[2]*mom[2]));
      m_hThetaMC->fill(theta);
      map_fit_flag[0][particles.back()] = 0;
      map_fit_flag[1][particles.back()] = 0;
    }
  }
  
  for(int i=0; i<2; i++){
    m_dInvPt[i] = 999;
    m_dD0[i]    = 999;
    if(trackCols[i]==nullptr) continue;
    int   iTrack = -1;
    for(auto track : *trackCols[i]){
      float D0, phi, omega, Z0, tanLambda, p, theta;
      float D0MC, phiMC, omegaMC, Z0MC, tanLambdaMC, pMC, thetaMC;
      float sigma_D0, sigma_phi, sigma_omega, sigma_Z0, sigma_tanLambda;
      // since possible more than one location=1 TrackState (not deleted in reconstruction), always use last one 
      for(std::vector<edm4hep::TrackState>::const_iterator it=track.trackStates_end()-1; it!=track.trackStates_begin()-1; it--){
	edm4hep::TrackState trackState = *it;
	if(trackState.location!=1)continue;
	D0 = trackState.D0;
	phi = trackState.phi;
	if(phi>CLHEP::pi) phi = phi - CLHEP::twopi;
	omega = trackState.omega;
	Z0 = trackState.Z0;
	tanLambda = trackState.tanLambda;
	sigma_D0 = std::sqrt(trackState.covMatrix[0]);
	sigma_phi = std::sqrt(trackState.covMatrix[2]);
	sigma_omega = std::sqrt(trackState.covMatrix[5]);
	sigma_Z0 = std::sqrt(trackState.covMatrix[9]);
        sigma_tanLambda = std::sqrt(trackState.covMatrix[14]);
	HelixClass helix_fit;
	helix_fit.Initialize_Canonical(phi,D0,Z0,omega,tanLambda,m_field);
	float px = helix_fit.getMomentum()[0];
	float py = helix_fit.getMomentum()[1];
	float pz = helix_fit.getMomentum()[2];
	p = sqrt(px*px+py*py+pz*pz);
	theta = acos(pz/p);
	//debug() << "TrackCollection " << i << " " << track.id() << ": "
	//	<<trackState.location << " " <<  D0 << " " << phi << " " << omega << " " << Z0 << " " << tanLambda << endmsg;
	m_nTrack[i]++;
	bool match = false;
	for(unsigned ip =0; ip<particles.size(); ip++){
	  edm4hep::ConstMCParticle particle = particles[ip];
	  debug() << ip << ":" << particle <<endmsg;
	  edm4hep::Vector3d pos = particle.getVertex();
	  edm4hep::Vector3f mom = particle.getMomentum();
	  HelixClass helix;
	  float posV[3] = {pos[0],pos[1],pos[2]};
	  float momV[3] = {mom[0],mom[1],mom[2]};
          helix.Initialize_VP(posV,momV,-1/*particle.getCharge()*/,m_field);
          phiMC = helix.getPhi0();
	  if(phiMC>CLHEP::pi) phiMC = phiMC - CLHEP::twopi;
	  float phiC = phi;
	  if(phiC-phiMC>3.1415926){
            phiC = phiC - 3.1415926*2;
          }
          else if(phiMC-phiC>3.1415926){
            phiC = phiC + 3.1415926*2;
          }
          D0MC = helix.getD0();
          omegaMC = helix.getOmega();
          Z0MC = helix.getZ0();
          tanLambdaMC = helix.getTanLambda();
	  //debug() << "MCParticle   " << ip << " " << particle.id() << ": "
	  //	  << "   " <<  D0MC << " " << phiMC << " " << omegaMC << " " << Z0MC << " " << tanLambdaMC << endmsg;
	  if(fabs(D0-D0MC)/sigma_D0>m_cut||fabs(phiC-phiMC)/sigma_phi>m_cut||fabs(omega-omegaMC)/sigma_omega>m_cut||
	     fabs(Z0-Z0MC)/sigma_Z0>m_cut||fabs(tanLambda-tanLambdaMC)/sigma_tanLambda>m_cut){
	    debug() << fabs(D0-D0MC)/sigma_D0 << " " << fabs(phiC-phiMC)/sigma_phi << " " << fabs(omega-omegaMC)/sigma_omega << " "
	    	    << fabs(Z0-Z0MC)/sigma_Z0 << " " << fabs(tanLambda-tanLambdaMC)/sigma_tanLambda << endmsg;
	    continue;
	  }
	  //debug() << "matched!" << endmsg;
	  match = true;
	  pMC = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
	  thetaMC = acos(mom[2]/pMC);
	  if(map_fit_flag[i].find(particle)!=map_fit_flag[i].end()){
	    map_fit_flag[i][particle]++;
	    float dInvPt = 1/(p*sin(theta))-1/(pMC*sin(thetaMC));
	    if(fabs(dInvPt)<fabs(m_dInvPt[i])){
	      m_dInvPt[i] = dInvPt;
	    }
	    if(fabs(D0-D0MC)<fabs(m_dD0[i])){
	      m_dD0[i] = D0-D0MC;
	    }
	  }
	}
	if(match){
	  if(phi-phiMC>3.1415926){
            phi = phi - 3.1415926*2;
          }
          else if(phiMC-phi>3.1415926){
            phi = phi + 3.1415926*2;
          }
	  m_hPullD0[i]        ->fill((D0-D0MC)/sigma_D0);
	  m_hPullPhi0[i]      ->fill((phi-phiMC)/sigma_phi);
	  m_hPullOmega[i]     ->fill((omega-omegaMC)/sigma_omega);
	  m_hPullZ0[i]        ->fill((Z0-Z0MC)/sigma_Z0);
	  m_hPullTanLambda[i] ->fill((tanLambda-tanLambdaMC)/sigma_tanLambda);
	  HelixClass helix_fit;
	  helix_fit.Initialize_Canonical(phi,D0,Z0,omega,tanLambda,m_field);
	  float px = helix_fit.getMomentum()[0];
	  float py = helix_fit.getMomentum()[1];
	  float pz = helix_fit.getMomentum()[2];
	  p = sqrt(px*px+py*py+pz*pz);
	  theta = acos(pz/p);
	  m_hDiffPR[i] ->fill((p-pMC)/pMC);//(p*sin(theta)-pMC*sin(thetaMC));//((p-pMC)/pMC);
	  m_hDiffInvPt[i]->fill(1/(p*sin(theta))-1/(pMC*sin(thetaMC)));
	  m_hDiffD0[i]->fill(D0-D0MC);
	  m_nMatched[i]++;
	}
	break;
      }
    }
    for(std::map<edm4hep::ConstMCParticle, int>::iterator it=map_fit_flag[i].begin();it!=map_fit_flag[i].end();it++){
      if(it->second==0)continue;
      edm4hep::ConstMCParticle particle = it->first;
      edm4hep::Vector3d pos = particle.getVertex();
      edm4hep::Vector3f mom = particle.getMomentum();
      float pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      if(it->second>1){
	warning() << "More than one track matched to this particle " << endmsg;
	warning() << pos << mom << endmsg;
      }
      m_hPt[i]->fill(pt);
      float theta = acos(mom[2]/sqrt(pt*pt+mom[2]*mom[2]));
      m_hTheta[i]->fill(theta);
      m_nFilled[i]++;
    }
  }
  m_tuple->write();
  _nEvt++;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode CompareTrackAlg::finalize(){
  debug() << "Finalizing..." << endmsg;

  info() << "Type  Fit  Match  Fill" << endmsg;  
  for(int i=0;i<2;i++){
    info() << i << " :   " << m_nTrack[i] << " " << m_nMatched[i] << " " << m_nFilled[i] << endmsg; 
  }
  /* //AIDA/IHistogram does not support set, move to root analysis macro
  for(int i=0;i<=100;i++){
    if(m_hPtMC->GetBinContent(i)!=0){
      m_hEffPt->SetBinContent(i,double(m_hPt->GetBinContent(i))/double(m_hPtMC->GetBinContent(i)));
      m_hEffPt->SetBinError(i,sqrt(m_hPt->GetBinContent(i)*(1.0-m_hPt->GetBinContent(i))/m_hPtMC->GetBinContent(i)));
    }
    else{
      m_hEffPt->SetBinContent(i,0);
      m_hEffPt->SetBinError(i,0);
    }
    if(m_hThetaMC->GetBinContent(i)!=0){
      m_hEffTheta->SetBinContent(i,double(m_hTheta->GetBinContent(i))/double(m_hThetaMC->GetBinContent(i)));
      m_hEffTheta->SetBinError(i,sqrt(m_hTheta->GetBinContent(i)*(1.0-m_hTheta->GetBinContent(i))/m_hThetaMC->GetBinContent(i)));
    }
    else{
      m_hEffTheta->SetBinContent(i,0);
      m_hEffTheta->SetBinError(i,0);
    }
  }
  */
  return StatusCode::SUCCESS;
}
