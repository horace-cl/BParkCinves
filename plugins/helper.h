#ifndef PhysicsTools_BParkingNano_helpers
#define PhysicsTools_BParkingNano_helpers

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryVector/interface/PV3DBase.h"
#include "Math/LorentzVector.h"

#include <vector>
#include <algorithm>
#include <limits>
#include <memory>

#include "TVector3.h"
#include "TMatrixD.h"

typedef std::vector<reco::TransientTrack> TransientTrackCollection;

constexpr float K_MASS = 0.493677;
constexpr float PI_MASS = 0.139571;
constexpr float LEP_SIGMA = 0.0000001;
constexpr float K_SIGMA = 0.000016;
constexpr float PI_SIGMA = 0.000016;
constexpr float MUON_MASS = 0.10565837;
constexpr float ELECTRON_MASS = 0.000511;

inline std::pair<float, float> min_max_dr(const std::vector< edm::Ptr<reco::Candidate> > & cands) {
  float min_dr = std::numeric_limits<float>::max();
  float max_dr = 0.;
  for(size_t i = 0; i < cands.size(); ++i) {
    for(size_t j = i+1; j < cands.size(); ++j) {
      float dr = reco::deltaR(*cands.at(i), *cands.at(j));
      min_dr = std::min(min_dr, dr);
      max_dr = std::max(max_dr, dr);
    }
  }
  return std::make_pair(min_dr, max_dr);
}

template<typename FITTER, typename LORENTZ_VEC>
inline double cos_theta_2D(const FITTER& fitter, const reco::BeamSpot &bs, const LORENTZ_VEC& p4) {
  if(!fitter.success()) return -2;
  GlobalPoint point = fitter.fitted_vtx();
  auto bs_pos = bs.position(point.z());
  math::XYZVector delta(point.x() - bs_pos.x(), point.y() - bs_pos.y(), 0.);
  math::XYZVector pt(p4.px(), p4.py(), 0.);
  double den = (delta.R() * pt.R());
  return (den != 0.) ? delta.Dot(pt)/den : -2;
}

template<typename FITTER>
inline Measurement1D l_xy(const FITTER& fitter, const reco::BeamSpot &bs) {
  if(!fitter.success()) return {-2, -2};
  GlobalPoint point = fitter.fitted_vtx();
  GlobalError err = fitter.fitted_vtx_uncertainty();
  auto bs_pos = bs.position(point.z());
  GlobalPoint delta(point.x() - bs_pos.x(), point.y() - bs_pos.y(), 0.);  
  return {delta.perp(), sqrt(err.rerr(delta))};
}


inline GlobalPoint FlightDistVector (const reco::BeamSpot & bm, GlobalPoint Bvtx)
{
   GlobalPoint Dispbeamspot(-1*( (bm.x0()-Bvtx.x()) + (Bvtx.z()-bm.z0()) * bm.dxdz()),
			   -1*( (bm.y0()-Bvtx.y()) + (Bvtx.z()-bm.z0()) * bm.dydz()), 
                            0);                    
   return std::move(Dispbeamspot);
}


inline float CosA(GlobalPoint & dist, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> & Bp4)
{
    math::XYZVector vperp(dist.x(),dist.y(),0);
    math::XYZVector pperp(Bp4.Px(),Bp4.Py(),0); 
    return std::move(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
}


inline std::pair<double,double> computeDCA(const reco::TransientTrack& trackTT,
					   const reco::BeamSpot& beamSpot)
{
  double DCABS    = -1.;
  double DCABSErr = -1.;

  TrajectoryStateClosestToPoint theDCAXBS = 
    trackTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
  if (theDCAXBS.isValid()) {
    DCABS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
    DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  }

  return std::make_pair(DCABS,DCABSErr);
}


inline bool track_to_lepton_match(edm::Ptr<reco::Candidate> l_ptr, auto iso_tracks_id, unsigned int iTrk)
{
  for (unsigned int i = 0; i < l_ptr->numberOfSourceCandidatePtrs(); ++i) {
    if (! ((l_ptr->sourceCandidatePtr(i)).isNonnull() && 
           (l_ptr->sourceCandidatePtr(i)).isAvailable())
           )   continue;
    const edm::Ptr<reco::Candidate> & source = l_ptr->sourceCandidatePtr(i);
    if (source.id() == iso_tracks_id && source.key() == iTrk){
      return true;
    }        
  }
  return false;
}


inline float V0_Lifetime(TVector3 pv, TVector3 sv, TMatrixD EPV, TMatrixD ESV, float M, TVector3 pT, double &ct, double &ect)
{
  // NOTA1: Esta funcion calucla el tiempo de vida y su error usando la longitud propia de decaimiento transversal
  // Recuerde que estamos asumiendo que el error en el pt es despreciable por eso las matrices asociadas a este las definimos cero

  TVector3 svT(sv.X(),sv.Y(),0.0);
  TVector3 pvT(pv.X(),pv.Y(),0.0);
  TVector3 d = svT - pvT;

  TMatrixD VSV(2,2);
  VSV(0,0) = ESV(0,0);
  VSV(1,1) = ESV(1,1);
  VSV(1,0) = ESV(1,0);
  VSV(0,1) = VSV(1,0);

  TMatrixD VPV(2,2);
  VPV(0,0) = EPV(0,0);
  VPV(1,1) = EPV(1,1);
  VPV(1,0) = EPV(1,0);
  VPV(0,1) = VPV(1,0);

  TMatrixD VL(2,2); VL = VSV; VL+=VPV;

  TVector3 p = pT;

  TMatrixD VP(2,2);
  VP(0,0) = 0.0;
  VP(1,1) = 0.0;
  VP(0,1) = 0.0;
  VP(1,0) = 0.0;

  double Lxy = d.Dot(p)/p.Mag();
  double lf = Lxy*M/p.Mag();
  //cout<<" ---> "<<lf<<endl;
  ct = lf;
  
  //Ahora calaculamos el error en el tiempo de vida
  
  //computing Mass error
  //double sM2 = 0; //We assume 0 for now
  
  //computing Lxy error
  
  //Defining Matrix:
  TMatrixD A(2,2);
  TMatrixD B(2,2);
  TMatrixD C(2,2);
  TMatrixD EP(2,2);
  TMatrixD EL(2,2);
  
  //Aij = PiPj/p2
  //Bij = LiLj/Lxy2 (Li = SVi - PVi)
  //EPij = Vij(P)/p2
  //ELij = Vij(L)/Lxy^2;
  //Cij = LiPj/(pLxy)
  
  A(0,0) = p.X()*p.X()/p.Mag2();
  A(1,1) = p.Y()*p.Y()/p.Mag2();
  A(0,1) = p.X()*p.Y()/p.Mag2();
  A(1,0) = A(0,1);

  B(0,0) = d.X()*d.X()/(Lxy*Lxy);
  B(1,1) = d.Y()*d.Y()/(Lxy*Lxy);
  B(0,1) = d.X()*d.Y()/(Lxy*Lxy);
  B(1,0) = B(0,1);
  
  C(0,0) = d.X()*p.X()/(Lxy*p.Mag());
  C(1,1) = d.Y()*p.Y()/(Lxy*p.Mag());
  C(0,1) = d.X()*p.Y()/(Lxy*p.Mag());
  C(1,0) = d.Y()*p.X()/(Lxy*p.Mag());

  EP = VP;
  EP*= ((double)1.0/p.Mag2());
  EL = VL;
  EL*= ((double)1.0/(Lxy*Lxy));

  //Test
  //EL(0,1) = 0.0;
  //EL(1,0) = 0.0;

  //Calculando Sigma Lxy
  // Sigma Lxy^2 = Tr{A*EL + (B + 4*A - 4*C)*EP}
  // NOTA2: en nuestro caso basicamente es Sigma Lxy^2 = Tr{A*EL), dado que no consideramos el momentum P
  
  TMatrixD A1 = A;
  A1*=(double)4.0;
  A1+=B;
  TMatrixD C1 = C;
  C1*=(double)4.0;
  A1-=C1;
  
  TMatrixD A_EL(A,TMatrixD::kMult,EL);
  TMatrixD A1_EP(A1,TMatrixD::kMult,EP);
  TMatrixD SL = A_EL;SL+=A1_EP;
  double sLxy2 = SL(0,0) + SL(1,1); 
  
  return ect = (double) fabs(lf)*sqrt(sLxy2);
}
 
#endif
