#include "TrackerDesignToolAlg.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Exceptions.h"
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
//#include "DDRec/Material.h"
#include "DD4hep/DD4hepUnits.h"

#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TVirtualGeoTrack.h"
#include "TVirtualMagField.h"
#include "TGeoGlobalMagField.h"
#include "TEveTrackPropagator.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveManager.h"
#include "TRint.h"
#include "TMatrixTSym.h"
//#include "CLHEP/Units/SystemOfUnits.h"
#include "TMath.h"
#include "TGeoSystemOfUnits.h"
#include <math.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
//#include <iostream>

DECLARE_COMPONENT( TrackerDesignToolAlg )

//------------------------------------------------------------------------------
TrackerDesignToolAlg::TrackerDesignToolAlg( const std::string& name, ISvcLocator* pSvcLocator )
 : Algorithm( name, pSvcLocator ) {
  declareProperty("MCParticleCollection", m_mcColHdl, "Handle of the Input MCParticle collection");
}

//------------------------------------------------------------------------------
StatusCode TrackerDesignToolAlg::initialize(){
  if(m_mapSens.size()!=m_mapResU.size()||m_mapSens.size()!=m_mapResV.size()){
    fatal() << "sensitive number not equal to set number for resolutions!" << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "Booking Ntuple" << endmsg;
  std::string tupleName = "MyTuples/TDT";
  NTuplePtr nt1(ntupleSvc(), tupleName.c_str());
  if ( !nt1 ) {
    m_tuple = ntupleSvc()->book(tupleName.c_str(),CLID_ColumnWiseTuple,"result");
    if ( 0 != m_tuple ) {
      //m_tuple->addItem        ("nmc",                 m_nParticles, 0, 1000 ).ignore();
      //m_tuple->addIndexedItem ("vx",                  m_nParticles, vx                 ).ignore();
    }
    else { // did not manage to book the N tuple....
      fatal() << "Cannot book " + tupleName <<endmsg; 
      return StatusCode::FAILURE;
    }
  }
  else{
    m_tuple = nt1;
  }

  auto geomSvc = service<IGeomSvc>("GeomSvc");
  if ( !geomSvc ) {
    info() << "Failed to find GeomSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  dd4hep::Detector* description = geomSvc->lcdd();

  try{
    //TRint *theApp = new TRint("ROOT example", 0, 0);
    //TEveManager::Create();
    //TDTHitVec hits = findHits(description, pt[i]*TGeoUnit::GeV, 75.*TMath::DegToRad(), 0*TMath::DegToRad());
    double p[9] = {1,3,6,10,20,30,50,80,100};
    double sigma[5][9],sigmaMSOff[5][9];
    double theta[2] = {85.*TMath::DegToRad(), 35.*TMath::DegToRad()};
    std::string name[5]={"#sigmaPt/Pt","#sigmaPhi","#sigmaD0","#sigmaCtg","#sigmaZ0"};
    for(int itheta = 0; itheta<1; itheta++){
      for(int i=0;i<9;i++){
	double pt = p[i]*TGeoUnit::GeV*sin(theta[itheta]);
	TDTHitVec hits = findHits(description, pt, theta[itheta], 90*TMath::DegToRad());
	
	auto covMatrixRZ = computeCovarianceMatrixRZ(hits, pt, theta[itheta]);
	sigma[3][i] = sqrt(covMatrixRZ(0,0));
	sigma[4][i] = sqrt(covMatrixRZ(1,1));
	auto covMatrixRPhi = computeCovarianceMatrixRPhi(hits, pt, theta[itheta]);
	sigma[0][i] = sqrt(covMatrixRPhi(0,0))*fabs(pt/TGeoUnit::MeV/(0.3*m_magField)*TGeoUnit::mm);
	sigma[1][i] = sqrt(covMatrixRPhi(1,1));
	sigma[2][i] = sqrt(covMatrixRPhi(2,2));
	auto covMatrixRZOff = computeCovarianceMatrixRZ(hits, pt, theta[itheta], false);
	sigmaMSOff[3][i] = sqrt(covMatrixRZOff(0,0));
	sigmaMSOff[4][i] = sqrt(covMatrixRZOff(1,1));
	auto covMatrixRPhiOff = computeCovarianceMatrixRPhi(hits, pt, theta[itheta], false);
	sigmaMSOff[0][i] = sqrt(covMatrixRPhiOff(0,0))*fabs(pt/TGeoUnit::MeV/(0.3*m_magField)*TGeoUnit::mm);
	sigmaMSOff[1][i] = sqrt(covMatrixRPhiOff(1,1))*TMath::RadToDeg();
	sigmaMSOff[2][i] = sqrt(covMatrixRPhiOff(2,2));
      }
      for(int i=0;i<5;i++){
	std::cout << name[i] << " MSON  ";
	for(int j=0;j<9;j++){
	  std::cout << sigma[i][j] << ", "; 
	}
	std::cout << std::endl;
	std::cout << name[i] << " MSOFF ";
	for(int j=0;j<9;j++){
	  std::cout << sigmaMSOff[i][j] << ", ";
	}
	std::cout << std::endl;
      }
    }
  }
  catch(std::runtime_error& e){
    
  }

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode TrackerDesignToolAlg::execute(){
  const edm4hep::MCParticleCollection* mcCols = nullptr;
  try {
    mcCols = m_mcColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << m_mcColHdl.fullKey() << " is unavailable" << endmsg;
  }
  
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode TrackerDesignToolAlg::finalize(){
  debug() << "Finalizing..." << endmsg;
  
  return StatusCode::SUCCESS;
}

TDTHitVec TrackerDesignToolAlg::findHits(dd4hep::Detector* description, double pt, double theta, double phi){
  //if(description->world().volume()->GetField()) description->world().volume()->GetField()->Print();
  //else info() << "no field set" << endmsg;
  TGeoManager* tgeoMgr = description->world().volume()->GetGeoManager();
  //TVirtualMagField* field = TGeoGlobalMagField::Instance()->GetField();
  //double pos[3]={0,0,0};
  //double B[3];
  //if(field){
  //  field->Field(pos,B);
  //  info() << B[0] << " " << B[1] << " "<< B[2] << endmsg;
  //}
  //else{
  //  info() << "no field set" << endmsg;
  //}

  //auto list = new TEveTrackList();
  //auto prop = list->GetPropagator();
  //prop->SetFitDaughters(kFALSE);
  //prop->SetMaxZ(1000);

  //prop->SetStepper(TEveTrackPropagator::kRungeKutta);
  //list->SetName("RK Propagator");


  //TEveTrack *track = 0;
  //prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 3.));
  //list->SetElementName(Form("%s, constB", list->GetElementName()));
  //auto rc = new TParticle(11,0,0,0,0,0,100,0,0,100.0000966049533374149533793973,0,0,0,0);

  //track = new TEveTrack(rc, -1, prop);
  //track->SetName(Form("Charge %d", -1));
  
  //list->SetLineColor(kMagenta);
  //track->SetLineColor(list->GetLineColor());

  //gEve->AddElement(list);
  //list->AddElement(track);

  //track->MakeTrack();

  //TGeoUniformMagField* field0 = new TGeoUniformMagField(0,0,3*TGeoUnit::tesla);   
  //TGeoGlobalMagField::Instance()->SetField(field0);
  //info() << "==========" << endmsg;
  TVector3 start(0,0,0);
  //double endpoint[3];
  TVector3 momentum(pt*cos(phi),pt*sin(phi),pt/tan(theta));
  TVector3 direction = momentum.Unit();
  //double direction[3] = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
  double L=0;

  tgeoMgr->AddTrack(0, 11);

  TGeoNode *node1 = tgeoMgr->InitTrack(start[0],start[1],start[2],direction[0],direction[1],direction[2]);
  if(!node1)
    throw std::runtime_error("No geometry node found at given location. Either there is no node placed here or position is outside of top volume.");

  double R = fabs(pt/TGeoUnit::MeV/(0.3*m_magField)*TGeoUnit::mm);
  double Phi0 = TGeoUnit::halfpi + phi;
  TVector3 point2center(R*sin(phi), -R*cos(phi), 0);
  TVector3 center = start + point2center;
  double step = 0.01*TGeoUnit::mm;
  double step_z = step*cos(theta);
  double step_rphi = step*sin(theta);
  //double p = pt/sin(theta);
  TVector3 currentPos = start;
  TVector3 previousPos = start;
  info() << "pt=" << pt << " center=" << center[0] << "," << center[1] << "," << center[2] << endmsg;
  TDTHitVec hits;
  double totalXtoX0 = 0;
  double volumeLength = 0;
  double totalLength = 0;
  double totalDeltaPhi = 0;
  for(int istep=0; istep<500000; istep++){
    TGeoNode * node2 = tgeoMgr->FindNextBoundaryAndStep( step/TGeoUnit::cm, 1) ;
    if( !node2 || tgeoMgr->IsOutside() ){
      info() << node2 << " " << tgeoMgr->IsOutside() << endmsg;
      break;
    }
    const double *position    =  tgeoMgr->GetCurrentPoint();
    if(sqrt(position[0]*position[0]+position[1]*position[1])>1900*TGeoUnit::mm||fabs(position[2])>2300*TGeoUnit::mm){
      //info() << position[0] << "," << position[1] << "," << position[2] << endmsg;
      break;
    }
    const double *previouspos =  tgeoMgr->GetLastPoint();
    double length = tgeoMgr->GetStep();
    //info() << "length=" << length << " " << previouspos[0] << "," << previouspos[1] << "," << previouspos[2]
    //	   << " -> " << position[0] << "," << position[1] << "," << position[2] << endmsg; 
    if(tgeoMgr->IsOnBoundary()){
      volumeLength += length;
      totalLength  += length;
      currentPos = TVector3(position);
      totalDeltaPhi = (start - center).Phi() - (currentPos - center).Phi();
      //currentPhi = (currentPos - center).Phi();
      direction.SetPhi(Phi0 - totalDeltaPhi -TGeoUnit::halfpi);
      //info() << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << "  " << direction[0] << "," << direction[1] << "," << direction[2] << endmsg;
      tgeoMgr->SetCurrentDirection(direction[0],direction[1],direction[2]);
    }
    else{
      if(fabs(length/step-1)>0.01) info() << length << " " << step << endmsg;
      volumeLength  += step;
      totalLength   += step;
      totalDeltaPhi += step_rphi/R;
      //TRotation rot;
      //rot.RotateZ(-totalDeltaPhi);
      //currentPhi -= step_dphi;
      TVector3 centerToPos = start-center;
      //info() << "before " << centerToPos[0] << "," << centerToPos[1] << "," << centerToPos[2] << endmsg; 
      centerToPos.RotateZ(-totalDeltaPhi);
      //info() << "after " << centerToPos[0] << "," << centerToPos[1] << "," << centerToPos[2] << endmsg;
      currentPos = center + centerToPos + TVector3(0,0,totalLength*cos(theta));// TVector3(R*cos(currentPhi),R*sin(currentPhi),totalLength*cos(theta));
      direction.SetPhi(Phi0 - totalDeltaPhi - TGeoUnit::halfpi);
      //info() << totalDeltaPhi << " " << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << "  " << direction[0] << "," << direction[1] << "," << direction[2] << endmsg;
      tgeoMgr->SetCurrentPoint(currentPos[0],currentPos[1],currentPos[2]);
      tgeoMgr->SetCurrentDirection(direction[0],direction[1],direction[2]);
      node1 = tgeoMgr->InitTrack(currentPos[0],currentPos[1],currentPos[2],direction[0],direction[1],direction[2]);
      if(node1!=node2) info() << "not OnBoundary, next node 1 != node 2" << endmsg;
      continue;
    }

    //TVirtualGeoTrack *track = tgeoMgr->GetLastTrack();
    //track->AddPoint( position[0], position[1], position[2], 0.);
    info() << "length=" << volumeLength << " " << previousPos[0] << "," << previousPos[1] << "," << previousPos[2]
           << " -> " << currentPos[0] << "," << currentPos[1] << "," << currentPos[2]
    	   << "    " << direction[0] << "," << direction[1] << "," << direction[2]
	   << " total = " << totalLength << endmsg;
    TGeoMedium* medium = node1->GetMedium();
    double XtoX0 = volumeLength/medium->GetMaterial()->GetRadLen();
    double radius = sqrt(position[0]*position[0]+position[1]*position[1]);
    if(radius<190) totalXtoX0 += XtoX0;

    double sinAlpha = radius/2./R;
    double corrRPhi = 1./sqrt(1.-sinAlpha*sinAlpha);
    //std::cout << "fucd debug: pt=" << pt << " " << XtoX0 << "*" << corrRPhi << " r=" << radius << " l=" << volumeLength
    //	      << " X0=" << medium->GetMaterial()->GetRadLen() 
    //	      << " + " << XtoX0 << " " << totalXtoX0 << std::endl;
    //XtoX0 *= corrRPhi;

    dd4hep::Volume volume = node1->GetVolume();
    bool isSens = true;//volume.isSensitive();
    std::string name = volume->GetName();
    //std::cout << "fucd debug: " << name << std::endl;
    double resU=0, resV=0;
    if(isSens){
      for(unsigned i=0; i<m_mapSens.size(); i++){
        if(name.find(m_mapSens[i])!=-1){
          resU = m_mapResU[i];
          resV = m_mapResV[i];
	}
      }
      for(unsigned i=0; i<m_mapSens.size(); i++){
	if(m_mapSens[i]==name){
	  resU = m_mapResU[i];
	  resV = m_mapResV[i];
	}
      }
      if(resU==0&&resV==0) isSens = false;
    }
    
    node1 = tgeoMgr->InitTrack(currentPos[0],currentPos[1],currentPos[2],direction[0],direction[1],direction[2]);
    if(node1!=node2){
      debug() << "node 1 != node 2" << node1->GetVolume()->GetName() << " " << node2->GetVolume()->GetName() << endmsg;
      TGeoNode * node2 = tgeoMgr->FindNextBoundaryAndStep( step/TGeoUnit::cm, 1) ;
      //node1 = node2;
      volumeLength = 0;
      previousPos = currentPos;
      //continue;
    }

    TVector3 normal(tgeoMgr->GetNormal()); 
    TDTHitPtr hit(new TDTHit(previousPos[0], previousPos[1], previousPos[2], position[0], position[1], position[2], XtoX0, isSens));
    //hit->SetDirection(direction);
    hit->SetPhi(totalDeltaPhi);
    hit->SetDeltaPhi(direction.Angle(normal));
    hit->SetPathLength(totalLength);
    if(isSens){
      hit->SetResolutionRPhi(resU*TGeoUnit::mm);
      hit->SetResolutionZ(resV*TGeoUnit::mm);
    }
    hits.push_back(std::move(hit));
    //info() << volume->GetName() << " " << resU << " " << resV << endmsg;

    volumeLength = 0;
    previousPos = currentPos;
    //node1 = tgeoMgr->InitTrack(currentPos[0],currentPos[1],currentPos[2],direction[0],direction[1],direction[2]);
    //if(node1!=node2) info() << "node 1 != node 2" << endmsg;
    node1 = node2;
  }

  for(int i=0;i<hits.size();i++){
    info() << i << ": r = " << hits.at(i)->GetRPos() << " length = " << hits.at(i)->GetPathLength() << " X0 = " << hits.at(i)->GetRadiation()
	    << " resRPhi = " << hits.at(i)->GetResolutionRPhi() << " resZ = " << hits.at(i)->GetResolutionZ() 
	    << endmsg;
  }

  tgeoMgr->ClearTracks();
  tgeoMgr->CleanGarbage();

  return hits;
}

TMatrixT<double> TrackerDesignToolAlg::computeCovarianceMatrixRPhi(TDTHitVec& hits, double pt, double theta, bool msOn){
  int nHits = hits.size();
  //std::cout << "fucd debug: nHits=" << nHits << std::endl;
  TMatrixT<double>      covMatrixRPhi(3,3);
  TMatrixTSym<double>   varMatrixRPhi(nHits,nHits);
  varMatrixRPhi.ResizeTo(nHits,nHits);
  varMatrixRPhi.Zero();

  int iStart          = nHits;
  int iEnd            = nHits-1;
  int nHitsUsed       = 0;
  int nActiveHitsUsed = 0;

  for(auto it=hits.rbegin(); it!=hits.rend(); it++){
    if((*it)->isSensitive()) nActiveHitsUsed++;
    iStart--;
  }
  if(nActiveHitsUsed<3){
    std::string message = "Variance matrix V(NxN) in R-Phi -> refPointRPos[mm]=";
    message+= " in combination with propagator direction doesn't provide sufficient number of hits for tracking!";
    info() << message << endmsg;
    return covMatrixRPhi;
  }
  nHitsUsed = nHits - iStart;

  std::vector<double> msThetaOverSinSq;
  //std::cout << "fucd debug: iStart=" << iStart << " nHitsUsed=" << nHitsUsed << " nActiveHitsUsed=" << nActiveHitsUsed << std::endl;
  for (int i=iStart; i<=iEnd; i++) {
    // MS theta
    double msTheta = 0.0;
    // Material in terms of rad. lengths
    double XtoX0 = hits.at(i)->GetRadiation();

    if(XtoX0>0){
      //msTheta = (13.6*dd4hep::MeV * 13.6*dd4hep::MeV) / (pt/dd4hep::MeV * pt/dd4hep::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));
      msTheta = (13.6 * 13.6) / (pt/TGeoUnit::MeV * pt/TGeoUnit::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));
      double A = hits.at(i)->GetRPos()/2./fabs(pt/TGeoUnit::MeV/(0.3*m_magField)*TGeoUnit::mm);     // r_i/2R
      double corrFactor = 1 + A*A*cos(theta)*cos(theta)/(1-A*A);
      //msTheta *= corrFactor;
      msTheta *= (cos(hits.at(i)->GetDeltaPhi())*cos(hits.at(i)->GetDeltaPhi()));
      //info() << i << ": " << hits.at(i)->GetRPos() << " " << XtoX0 << " " << msTheta/corrFactor << "*" << corrFactor << endmsg;
    }
    else {
      msTheta = 0;
    }
    msThetaOverSinSq.push_back(msTheta);
  }
  //std::cout << "fucd debug: " << msThetaOverSinSq.size() << std::endl;
  for (int c=iStart ; c<=iEnd; c++) {

    // Dummy value for correlations involving inactive surfaces
    if (hits.at(c)->isPassive()) {
      for (int r = 0; r <= c; r++) varMatrixRPhi(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r=iStart; r <= c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (hits.at(r)->isPassive()) varMatrixRPhi(r, c) = 0.0;

        // Correlations between two active surfaces
        else {

          double sum = 0.0;

	  if(msOn){
	    for (int i=iStart; i<r; i++) {
	      //std::cout << i-(nHits-nHitsUsed) << " = " << i << "-(" << nHits << "-" << nHitsUsed << ")" << std::endl;
	      sum += msThetaOverSinSq.at(i-(nHits-nHitsUsed)) * (hits.at(c)->GetRPos() - hits.at(i)->GetRPos()) * (hits.at(r)->GetRPos() - hits.at(i)->GetRPos());
	      double r1 = hits.at(i)->GetRPos();
	      double r2 = hits.at(c)->GetRPos();
	      double deltaPhi = hits.at(c)->GetPhi() - hits.at(i)->GetPhi();
	      //sum += msThetaOverSinSq.at(i-(nHits-nHitsUsed)) * (hits.at(c)->GetPathLength() - hits.at(i)->GetPathLength()) * (hits.at(r)->GetPathLength() - hits.at(i)->GetPathLength());
	      //info() << r1 << " " << r2 << " " << deltaPhi << " " << (r1*r1+r2*r2+2*r1*r2*cos(deltaPhi/2)) << endmsg;
	    }
	    //if (r == c)
	      //sum += msThetaOverSinSq.at(r-(nHits-nHitsUsed)) * (hits.at(r)->GetPathLength() - hits.at(r-1)->GetPathLength()) * (hits.at(r)->GetPathLength() - hits.at(r-1)->GetPathLength());;
	  }
	  if (r == c) {
            // Get hit
            //const Hit* const myHit = hits.at(r).get();

            // Get the hit's Rphi resolution
            //const double trackRadius = getRadius(myHit->getZPos());
            const double resolutionRPhi = hits.at(r)->GetResolutionRPhi(); //myHit->getResolutionRphi(trackRadius, deltaTheta);

            // Add resolutionRPhi to the correlation matrix
            sum = sum + pow(resolutionRPhi, 2.);
	    //info() << r << ": " << sum << endmsg;
          }
	  //std::cout << "fucd debug2: sum=" << sum << " pt=" << pt/TGeoUnit::MeV << std::endl;
          varMatrixRPhi(r, c) = sum;
          if (r != c) varMatrixRPhi(c, r) = sum;
        }
      }
    }
  }
  int  rActual = -1;          // Row, at which to move the active row due to a sequence of zero rows or inactive hits inbetween
  bool lookForActive = false; // Start looking for shift of active rows, after first passive row found
  //std::cout << "fucd debug:================" << std::endl;
  for (int r=0; r<nHits; r++) {
    if ((hits.at(r)->isPassive() || r<iStart) && (!lookForActive)) {
      // Next hit has to be active (set as active and considered in track fitting (see iStart))
      if ((r+1)<nHits && hits.at(r+1)->isSensitive() && (r+1)>=iStart) lookForActive = true;
      
      // Previous hit has to be passive or not being considered in track fitting (see iStart))
      if (!((r-1)>=0 && (hits.at(r-1)->isPassive() || (r-1)<iStart)) ) rActual = r;
    }
    // Shift active layer to zero-th row + i active layers, which have already been shifted by number of zero layers
    else if ((hits.at(r)->isSensitive()) && (lookForActive)) {
      for (int c=0; c<nHits; c++) {
        varMatrixRPhi(rActual, c) = varMatrixRPhi(r, c);
        varMatrixRPhi(c, rActual) = varMatrixRPhi(c, r);
      }
      
      varMatrixRPhi(rActual, rActual) = varMatrixRPhi(r, r);
      rActual++;
    }
  }
  //std::cout << "fucd debug: rActual = " << rActual << std::endl; 
  // If some rows/colums were zero -> matrix rank needs to be adjusted
  int nResized = rActual;
  if (nResized!=0) varMatrixRPhi.ResizeTo(nResized, nResized);
  // Compute 3x3 covariance matrix of the track parameters in R-Phi projection
  unsigned int offset = iStart;

  int n = varMatrixRPhi.GetNrows();

  TMatrixT<double> V(varMatrixRPhi);
  TMatrixT<double> diffsT(3, n);
  TMatrixT<double> diffs(n, 3);
  // Set up partial derivative matrices diffs and diffsT -> using Karimaki approach & parabolic aproximations to define these matrices
  for (auto i=iStart; i<=iEnd; i++) {
    if (hits.at(i)->isSensitive()) {
      diffs(i - offset, 0) = 0.5*hits.at(i)->GetRPos()*hits.at(i)->GetRPos();// computeDfOverDRho(hits.at(i)->GetRPos(),hits.at(i)->GetZPos());
      diffs(i - offset, 1) = +hits.at(i)->GetRPos();
      diffs(i - offset, 2) = 1;
    }
    else offset++;
  }
  
  //for(int i=0; i<nResized; i++){
  //for(int j=0; j<nResized; j++){
  //  std::cout << V[i][j] << " ";
  //}
  //std::cout << std::endl;
  //} 
  //std::cout << diffsT << V << diffs << std::endl;
  // Transpose
  diffsT.Transpose(diffs);
  covMatrixRPhi = diffsT * V.Invert() * diffs;
  covMatrixRPhi.Invert();

  return covMatrixRPhi;
}


TMatrixT<double> TrackerDesignToolAlg::computeCovarianceMatrixRZ(TDTHitVec& hits, double pt, double theta, bool msOn){
  int nHits = hits.size();
  //std::cout << "fucd debug: nHits=" << nHits << std::endl;
  TMatrixT<double>      covMatrixRZ(2,2);
  TMatrixTSym<double>   varMatrixRZ(nHits,nHits);
  varMatrixRZ.ResizeTo(nHits,nHits);
  varMatrixRZ.Zero();

  int iStart          = nHits;
  int iEnd            = nHits-1;
  int nHitsUsed       = 0;
  int nActiveHitsUsed = 0;

  for(auto it=hits.rbegin(); it!=hits.rend(); it++){
    if((*it)->isSensitive()) nActiveHitsUsed++;
    iStart--;
  }
  if(nActiveHitsUsed<3){
    std::string message = "Variance matrix V(NxN) in R-Phi -> refPointRPos[mm]=";
    message+= " in combination with propagator direction doesn't provide sufficient number of hits for tracking!";
    info() << message << endmsg;
    return covMatrixRZ;
  }

  nHitsUsed = nHits - iStart;
  std::vector<double> msThetaOverSinSq;
  //std::cout << "fucd debug: iStart=" << iStart << " nHitsUsed=" << nHitsUsed << " nActiveHitsUsed=" << nActiveHitsUsed << std::endl;
  for (int i=iStart; i<iEnd; i++) {
    double msTheta = 0.0;
    double XtoX0 = hits.at(i)->GetRadiation();
    if (XtoX0>0) {
      //msTheta = (13.6*dd4hep::MeV * 13.6*dd4hep::MeV) / (pt/dd4hep::MeV * pt/dd4hep::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));
      msTheta = (13.6 * 13.6) / (pt/TGeoUnit::MeV * pt/TGeoUnit::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));
      double A = hits.at(i)->GetRPos()/2./fabs(pt/TGeoUnit::MeV/(0.3*m_magField)*TGeoUnit::mm);  // r_i/2R                                                          
      double corrFactor = pow( cos(theta)*cos(theta)/sin(theta)/sqrt(1-A*A) + sin(theta) ,2);
      //msTheta *=corrFactor;
      msTheta *= pow(1./sin(theta),2);
      //std::cout << "fucd debug: " << msTheta << " = " << msTheta/corrFactor << " x " << corrFactor << std::endl;  
    }
    else {
      msTheta = 0;
    }
    msThetaOverSinSq.push_back(msTheta);
  }
  //std::cout << "fucd debug: " << msThetaOverSinSq.size() << std::endl;
  for (int c=iStart ; c<=iEnd; c++) {
    if (hits.at(c)->isPassive()) {
      for (int r=0; r<=c; r++) varMatrixRZ(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {
      for (int r=iStart; r<=c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (hits.at(r)->isPassive()) varMatrixRZ(r, c) = 0.0;
        // Correlations between two active surfaces
        else {
          double sum = 0.0;
	  if(msOn){
	    for (int i=iStart; i<r; i++) sum += msThetaOverSinSq.at(i-(nHits-nHitsUsed))
					   * (hits.at(c)->GetRPos() - hits.at(i)->GetRPos())
					   * (hits.at(r)->GetRPos() - hits.at(i)->GetRPos());
	  }
          if (r == c) {
            double prec = hits.at(r)->GetResolutionZ();//fabs(pt/TGeoUnit::MeV/(0.3*m_magField)));
            sum = sum + prec * prec;
          }
	  //std::cout << "fucd debug1: sum=" << sum << " pt=" << pt/TGeoUnit::MeV << std::endl;
          varMatrixRZ(r, c) = sum;
          if (r != c) varMatrixRZ(c, r) = sum;
#undef CORRELATIONS_OFF_DEBUG
#ifdef CORRELATIONS_OFF_DEBUG
          if (r!=c) {
            varMatrixRZ(c, r)=0;
            varMatrixRZ(r, c)=0;
          }
#endif
	}
      }
    }
  }
  // Remove zero rows and columns in covariance matrix
  int  rActual = -1;          // Row, at which to move the active row due to a sequence of zero rows or inactive hits inbetween
  bool lookForActive = false; // Start looking for shift of active rows, after first passive row found
  
  for (int r=0; r<nHits; r++) {
    // Keep actual row @ zero for first N passive (zero) layers
    if ((hits.at(r)->isPassive() || r<iStart) && (!lookForActive)) {
      // Next hit has to be active (set as active and considered in track fitting (see iStart))
      if ((r+1)<nHits && hits.at(r+1)->isSensitive() && (r+1)>=iStart) lookForActive = true;
      // Previous hit has to be passive or not being considered in track fitting (see iStart))
      if (!((r-1)>=0 && (hits.at(r-1)->isPassive() || (r-1)<iStart)) ) rActual = r;
    }
    // Shift active layer to zero-th row + i active layers, which have already been shifted by number of zero layers
    else if ((hits.at(r)->isSensitive()) && (lookForActive)) {
      for (int c=0; c<nHits; c++) {
        varMatrixRZ(rActual, c) = varMatrixRZ(r, c);
        varMatrixRZ(c, rActual) = varMatrixRZ(c, r);
      }
      
      varMatrixRZ(rActual, rActual) = varMatrixRZ(r, r);
      rActual++;
    }
  }
  //std::cout << "fucd debug: rActual=" << rActual << std::endl;
  // If some rows/colums were zero -> matrix rank needs to be adjusted
  int nResized = rActual;
  if (nResized!=0) varMatrixRZ.ResizeTo(nResized, nResized);
  unsigned int offset       = iStart;
  unsigned int varMatrixDim = varMatrixRZ.GetNrows();

  TMatrixT<double> V(varMatrixRZ);        // Local copy to be inverted
  TMatrixT<double> diffsT(2, varMatrixDim); // Derivatives of track parameters transposed (in R-Z -> 2 track parameters)
  TMatrixT<double> diffs(varMatrixDim, 2);  // Derivatives of track parameters (in R-Z -> 2 track parameters)

  // Set up partial derivative matrices diffs and diffsT -> line fit in s-Z to define these matrices
  for (auto i=iStart; i<=iEnd; i++) {
    if (hits.at(i)->isSensitive()) {
      // Partial derivatives for x = p[0] * y + p[1]
      diffs(i - offset, 0) = hits.at(i)->GetRPos();
      diffs(i - offset, 1) = 1;
    }
    else offset++;
  }
  
  // Transpose
  diffsT.Transpose(diffs);

  // Print
  //std::cout << "Diff matrix in R-Z: " << std::endl;
  //printMatrix(diffsT);
  
  // Get covariance matrix using global chi2 fit: C = cov(i,j) = (D^T * V^-1 * D)^-1
  covMatrixRZ = diffsT * V.Invert() * diffs;
  covMatrixRZ.Invert();
  
  return covMatrixRZ;
}
