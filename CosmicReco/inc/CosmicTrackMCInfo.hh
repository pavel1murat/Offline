#ifndef RecoDataProducts_CosmicTrackMCInfo_hh
#define RecoDataProducts_CosmicTrackMCInfo_hh

#include "TMath.h"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/DataProducts/inc/GenVector.hh"

namespace mu2e {
struct CosmicTrackMCInfo{

     double TrueTheta = 0;
     double TruePhi = 0;
     double TruePhiSIM = 0;
     double TrueThetaSIM = 0;
     double TrueMomentum = 0;

     std::vector<double> TrueTimeResiduals;
     std::vector<double> TrueDOCA;
     std::vector<double> Ambig;
     TrackEquation TrueFitEquation;
     TrackAxes TrueTrackCoordSystem;
     TrackParams RawTrueParams;
     CosmicTrackMCInfo();
     CosmicTrackMCInfo(double theta, double phi, TrackEquation eqn, TrackAxes axis, TrackParams para) : TrueTheta(theta), TruePhi(phi), TrueFitEquation(eqn), TrueTrackCoordSystem(axis), RawTrueParams(para) {};

};

} // end namespace mu2e

#endif
