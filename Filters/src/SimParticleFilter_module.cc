// ======================================================================
//
// SimParticleMomFilter_module: allows filtering on the energy
//   of the SimParticle
//
// ======================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "MCDataProducts/inc/GenId.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"



#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class SimParticleFilter : public art::EDFilter {
    public:
      explicit SimParticleFilter(fhicl::ParameterSet const& pset);

    private:
      bool beginRun(art::Run& run) override;
      bool endRun(art::Run& run) override;
      bool endSubRun(art::SubRun& sr) override;
      bool filter(art::Event& event) override;
      
      art::InputTag _simpCollTag;
      double        _minP;
      PDGCode::type _pdgID;
      unsigned      _nevt, _npass;
  };

  SimParticleFilter::SimParticleFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _simpCollTag (pset.get<std::string>      ("simpCollTag"   ,"compressDigiMCs")),
    _minP        (pset.get<double>           ("minP          ",-1               )),
    _pdgID       (PDGCode::type(pset.get<int>("pdgID"         , 0)              )),
    _nevt(0), _npass(0){}

  bool SimParticleFilter::beginRun(art::Run& run) {
    return true;
  }

  bool SimParticleFilter::endRun(art::Run& run) {
    return true;
  }

  bool SimParticleFilter::endSubRun(art::SubRun& sr) {
    mf::LogInfo("INFO") << "GenEvents seen: "
			<< _nevt << ", GenEvents passed: " 
			<< _npass << " for " << sr.id()
			<< "\n";
    return true;

  }

  bool SimParticleFilter::filter(art::Event& event) {

    auto simpColl = event.getValidHandle<SimParticleCollection>(_simpCollTag);

    bool rc(false);

    // find highest momentum gen particle that passes cuts
    for ( const auto& i: *simpColl ) {
      const SimParticle* simp = &i.second;
      if (abs(simp->pdgId()) == _pdgID) {
        if (simp->startMomentum().vect().mag() > _minP) {
	  // there is at least one particle with requested PDG_ID and P > minP
	  rc = true; 
	  break;
	}
      }
    }

    ++_nevt;
    if (rc) ++_npass;

    return rc;
  }
}

using mu2e::SimParticleFilter;
DEFINE_ART_MODULE(SimParticleFilter);
