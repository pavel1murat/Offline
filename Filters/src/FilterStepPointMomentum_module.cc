// Pass events with at least one hit satisfying a min momentum cut.
//
// Andrei Gaponenko, 2013

#include <string>
#include <map>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "TH1.h"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {

  //================================================================
  class FilterStepPointMomentum : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags             inputTags_      ;
    double                cutMomentumMin_ ;
    double                cutMomentumMax_ ;

    SimParticleTimeOffset timeOffsets_    ; // time offsets
    double                tMin_           ; // if non-negative: hit timing cut-off
    double                tMax_           ; // 
    double                mbtime_         ; // microbunch length, ns
    uint                  volumeId_       ; // default : -1
    int                   pdgId_          ; // default :  0
    bool                  fillHistograms_ ; // 
    // statistics counters
    unsigned              numInputEvents_ ;
    unsigned              numPassedEvents_;

    struct Hist_t {
      TH1F*   time;                   // step point MC time
      TH1F*   tmax;                   // max time in the event
    } hist_ ;

    Hist_t fHist;
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> inputs { Name("inputs"),
          Comment("input StepPointMCCollections")
          };

      fhicl::Atom<double> cutMomentumMin { Name("cutMomentumMin"),
          Comment("The filter passes events if any of the step points satisties pmag>cutMomentumMin\n"
                  "By default cutMomentumMin=-inf\n"),
          -std::numeric_limits<double>::max()
          };

      fhicl::Atom<double> cutMomentumMax { Name("cutMomentumMax"),
          Comment("The filter passes events if any of the step points satisties pmag<cutMomentumMax\n"
                  "By default cutMomentumMax=inf\n"),
          std::numeric_limits<double>::max()
          };

      fhicl::Sequence<art::InputTag> timeOffsets { Name("TimeOffsets"), 
	  Comment("Sim Particle Time Offset Maps, default: NONE")
	  };

      fhicl::Atom<double> tMin { Name("tMin"),
          Comment("The filter passes events if any of the step points T>=tTMin, default:-1.e10\n") 
	  };

      fhicl::Atom<double> tMax { Name("tMax"),
          Comment("The filter passes events if any of the step points T<tTMin, default: 1.e10\n") 
	  };

      fhicl::Atom<int> volumeId { Name("volumeId"),
          Comment("if positive, defines the volume ID to look at, default: -1\n") 
	  };

      fhicl::Atom<int> pdgId { Name("pdgId"),
          Comment("if non-zero, defines the PDG ID to use, default: 0\n") 
	  };

      fhicl::Atom<bool> fillHistograms { Name("fillHistograms"),
          Comment("if true, fill histograms, default: false\n") 
	  };

    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit FilterStepPointMomentum(const Parameters& conf);

    void         fillHistograms(art::Event& event);

    virtual bool beginRun(art::Run&   run );
    virtual void beginJob() override;
    virtual bool filter  (art::Event& event) override;
    virtual void endJob  () override;
  };

  //================================================================
  FilterStepPointMomentum::FilterStepPointMomentum(const Parameters& conf)
    : art::EDFilter{conf}
    , cutMomentumMin_(conf().cutMomentumMin())
    , cutMomentumMax_(conf().cutMomentumMax())
    , timeOffsets_   (conf().timeOffsets()   )
    , tMin_          (conf().tMin()          )
    , tMax_          (conf().tMax()          )
    , volumeId_      (conf().volumeId()      )
    , pdgId_         (conf().pdgId()         )
    , fillHistograms_(conf().fillHistograms())
    , numInputEvents_ (0)
    , numPassedEvents_(0)
  {
    for(const auto& i : conf().inputs()) {
      inputTags_.emplace_back(i);
    }

  }

//-----------------------------------------------------------------------------
  void FilterStepPointMomentum::beginJob(){
    if (fillHistograms_) {
      art::ServiceHandle<art::TFileService> tfs;

      hist_.time  = tfs->make<TH1F>("time","step time"    ,200,-200,1800);
      hist_.tmax  = tfs->make<TH1F>("tmax","Max step time",200,-200,1800);
    }
  }

//-----------------------------------------------------------------------------
  bool FilterStepPointMomentum::beginRun( art::Run& run ){
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    mbtime_ = accPar->deBuncherPeriod; 
    return true;
  }

//-----------------------------------------------------------------------------
  void FilterStepPointMomentum::fillHistograms(art::Event& event) {
    float tmax (-1e6);

    for(const auto& cn : inputTags_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(cn);
      for (const auto& hit : *ih) {
	float t = hit.time()+timeOffsets_.totalTimeOffset(hit.simParticle());
	t  = fmod(t,mbtime_);
	if (t > tmax) tmax = t;

	hist_.time->Fill(t);
      }
    }
    hist_.tmax->Fill(tmax);
  }

  //================================================================
  bool FilterStepPointMomentum::filter(art::Event& event) {
    bool passed = false;

    if (tMin_ > 0) timeOffsets_.updateMap(event);

    for(const auto& cn : inputTags_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(cn);
      for(const auto& hit : *ih) {
	if ((pdgId_    != 0) and (hit.simParticle()->pdgId() != pdgId_   )) continue;
	if ((volumeId_ != 0) and (hit.volumeId()             != volumeId_)) continue;
	double p = hit.momentum().mag();
	if ((p <= cutMomentumMin_) or (p >= cutMomentumMax_))               continue;
	double t = hit.time(); 
	if (tMin_ > 0) {
	  // also check the hit time, wrap around the mucrobunch
	  
	  t  = t+timeOffsets_.totalTimeOffset(hit.simParticle());
	  t  = fmod(t,mbtime_);
	}
	  
	if ((t >= tMin_) && (t < tMax_)) {
	  passed = true;
	  break;
	}
      }
    }

    if (fillHistograms_) fillHistograms(event);

    ++numInputEvents_;
    if (passed) { ++numPassedEvents_; }
    return passed;
  }

  //================================================================
  void FilterStepPointMomentum::endJob() {
    mf::LogInfo("Summary")
      <<"FilterStepPointMomentum_module: passed "
      <<numPassedEvents_<<" / "<<numInputEvents_<<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointMomentum);
