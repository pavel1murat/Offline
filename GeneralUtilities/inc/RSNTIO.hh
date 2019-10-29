// Data structures used to dump info and then re-sample particles from
// a ROOT tree in some multi-stage job configurations.
//
// Andrei Gaponenko, 2015

#ifndef GeneralUtilities_inc_RSNTIO_hh
#define GeneralUtilities_inc_RSNTIO_hh

namespace mu2e {
  namespace IO {

    //================================================================
    struct StoppedParticleF {
      float x;
      float y;
      float z;
      float t;

      StoppedParticleF() : x(), y(), z(), t() {}

      static const char *branchDescription() {
        return "x/F:y/F:z/F:time/F";
      }

      static unsigned numBranchLeaves() { return 4; }
    };

    //================================================================
    struct StoppedParticleTauNormF {
      float x;
      float y;
      float z;
      float t;
      float tauNormalized;

      StoppedParticleTauNormF() : x(), y(), z(), t(), tauNormalized() {}

      static const char *branchDescription() {
        return "x/F:y/F:z/F:time/F:tauNormalized/F";
      }

      static unsigned numBranchLeaves() { return 5; }
    };

    //================================================================
    struct InFlightParticleD {
      double x;
      double y;
      double z;
      double time;
      double px;
      double py;
      double pz;
      int    pdgId;

      InFlightParticleD() : x(), y(), z(), time(), px(), py(), pz(), pdgId() {}

      static const char *branchDescription() {
        return "x/D:y/D:z/D:time/D:px/D:py/D:pz/D:pdgId/I";
      }

      static unsigned numBranchLeaves() { return 8; }
    };

    //================================================================
    //For resampling photon conversions
    struct ConversionPointF {
      float x;
      float y;
      float z;
      float time;
      float px;
      float py;
      float pz;
      float weight;
      float genEnergy;
      char  mat[40];

      ConversionPointF() : x(), y(), z(), time(), px(), py(), pz(), weight(), genEnergy(), mat() {}

      static const char *branchDescription() {
        return "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:weight/F:genEnergy/F:mat/C";
      }

      static unsigned numBranchLeaves() { return 10; }
    };

    //================================================================

  } // IO
} // mu2e

#endif/*GeneralUtilities_inc_RSNTIO_hh*/
