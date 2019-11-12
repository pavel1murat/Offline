// Michael MacKenzie, 2019
// Samples photon conversions stored and outputs e+e- pair conversions
// Assumes G4 is part of the trigger path, will cause an exception otherwise

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"

// G4 includes.
#include "G4Material.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4GammaConversion.hh"
#include "G4VEmModel.hh"
#include "G4BetheHeitlerModel.hh"
#include "G4DynamicParticle.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"

namespace mu2e {

  //================================================================
  class GammaConversionGun : public art::EDProducer {
  public:

    typedef RootTreeSampler<IO::ConversionPointF> RTS;
    
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("verbosity level (>=0)"), 0};
      fhicl::Table<RTS::Config> stops{Name("gammaStops"), Comment("Gamma stops parameter set")};
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), Comment("Whether or not to make generation histograms (true/false)"), false};
      fhicl::Atom<double> xMin{Name("xMin"), Comment("Stop minimum x value (mm)"), -1.e9};
      fhicl::Atom<double> xMax{Name("xMax"), Comment("Stop maximum x value (mm)"),  1.e9};
      fhicl::Atom<double> yMin{Name("yMin"), Comment("Stop minimum y value (mm)"), -1.e9};
      fhicl::Atom<double> yMax{Name("yMax"), Comment("Stop maximum y value (mm)"),  1.e9};
      fhicl::Atom<double> zMin{Name("zMin"), Comment("Stop minimum z value (mm)"), -1.e9};
      fhicl::Atom<double> zMax{Name("zMax"), Comment("Stop maximum z value (mm)"),  1.e9};
      fhicl::Atom<double> rMin{Name("rMin"), Comment("Stop minimum radius in the DS (mm)"), -1.};
      fhicl::Atom<double> rMax{Name("rMax"), Comment("Stop maximum radius in the DS (mm)"),  1.e9};
      fhicl::Atom<double> czMin{Name("czMin"), Comment("Stop minimum cos(theta) value"), -1.};
      fhicl::Atom<double> czMax{Name("czMax"), Comment("Stop maximum cos(theta) value"),  1.};
      fhicl::Atom<double> pMin{Name("pMin"), Comment("Minimum photon daughter momentum (MeV/c)"), -1.};
      fhicl::Atom<double> pMax{Name("pMax"), Comment("Maximum photon generated energy (MeV/c) ( < 0 to ignore)"),  -1.};
      fhicl::Atom<std::string> defaultMat{Name("defaultMaterial"), Comment("Override ntuple material with a given material"),  ""};
      fhicl::Atom<double> testE{Name("testE"), Comment("Test photon energy to override ntuple energy with (MeV/c) ( < 0 to ignore)"), -1.};
      fhicl::Atom<double> xOffset{Name("solenoidXOffset"), Comment("X coordinate offset for radius calculations (mm)"), -3904.};
    };
    typedef art::EDProducer::Table<Config> Parameters;


    explicit GammaConversionGun(const Parameters& pset);
    virtual void produce(art::Event& event);
    int                 verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;

    CLHEP::RandFlat     randomFlat_;
    RTS stops_;
    // PairProduction pairProd_;

    bool doHistograms_;
    //for restricted space generation
    double xMin_;
    double xMax_;
    double yMin_;
    double yMax_;
    double zMin_;
    double zMax_;
    double rMin_; //defined from the ds axis
    double rMax_; 
    //minimum momentum of a daughter
    double pMin_;
    //maximum momentum of a photon
    double pMax_;
    //photon cos(theta) restriction
    double czMin_;
    double czMax_;

    std::string defaultMat_; //for testing the pair production spectrum for a given material
    double testE_; //for testing the pair production spectrum for a given Energy
    
    double xOffset_; //ds axis x offset

    G4ParticleDefinition* photon_; //for passing into the pair production spectrum
    G4BetheHeitlerModel* g4bhm_; //pair production spectrum
    
    TH1D* _hgencuts; //records number of events attempted
    TH1F* _hmomentum;
    TH1F* _hcos;
    TH1F* _hcosWt;
    TH1F* _hr; //x-y r from (x,y) = (-3904, 0)
    TH1F* _hz; 
    TH1F* _hElecMom  {nullptr};
    TH1F* _hPosiMom  {nullptr};
    TH1F* _hTotMom   {nullptr};
    TH1F* _hTotMomWt {nullptr};
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;                   // M(ee)/E(gamma)
    TH1F* _hy;				// splitting function
  };

  //================================================================
  GammaConversionGun::GammaConversionGun(const Parameters& pset)
    : EDProducer(pset)
    , verbosityLevel_  (pset().verbosityLevel())
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randomFlat_      (eng_)
    , stops_           (eng_, pset().stops())
    , doHistograms_    (pset().doHistograms())
    , xMin_            (pset().xMin())
    , xMax_            (pset().xMax())
    , yMin_            (pset().yMin())
    , yMax_            (pset().yMax())
    , zMin_            (pset().zMin())
    , zMax_            (pset().zMax())
    , rMin_            (pset().rMin())
    , rMax_            (pset().rMax())
    , pMin_            (pset().pMin())
    , pMax_            (pset().pMax())
    , czMin_           (pset().czMin())
    , czMax_           (pset().czMax())
    , defaultMat_      (pset().defaultMat())
    , testE_           (pset().testE())
    , xOffset_         (pset().xOffset())
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::EventWeight>();
    produces<GenParticle>("photon"); //store photon generation energy for RMC weights

    photon_ = 0;
    // objsize_ = sizeof(*photon_) + sizeof(*g4bhm_);

    if(verbosityLevel_ > 0) {
      std::cout<<"GammaConversionGun: using = "
               <<stops_.numRecords()
               <<" stopped gammas"
               <<std::endl;

      std::cout<<"GammaConversionGun: producing photon " << std::endl;
    }
    if(xMin_ >= xMax_ || yMin_ >= yMax_ || zMin_ >= zMax_ || rMin_ >= rMax_)
      throw cet::exception("BADCONFIG")
	<< "GammaConversionGun: Stop (x,y,z,r) restriction error! "
	<< xMin_ << " < x < " << xMax_ << ", "
	<< yMin_ << " < y < " << yMax_ << ", "
	<< zMin_ << " < z < " << zMax_ << ", "
	<< rMin_ << " < r < " << rMax_ << "\n";

    if(defaultMat_.size() > 0)
      std::cout << "GammaConversionGun: Overriding Ntuple defined material and instead using "
		<< defaultMat_.c_str() << std::endl;
    if(testE_ > 0) {
      std::cout << "GammaConversionGun: Overriding Ntuple defined energy and instead using "
		<< testE_ << std::endl;
      if(pMin_ > testE_)
	throw cet::exception("BADCONFIG")
	  << "GammaConversionGun: Test energy set to value below minimum electron/positron momentum!\n";
      if(pMax_ > 0. && pMax_ < testE_)
	throw cet::exception("BADCONFIG")
	  << "GammaConversionGun: Test energy set to value above maximum photon momentum!\n";

    }

    //essentially a weight histogram, so always store
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "GammaConversionGun" );
    _hgencuts  = tfdir.make<TH1D>( "hgencuts", "Attempts to generate vs cut number", 10,  0.5,  10.5  );
    if ( doHistograms_ ) {

      _hmomentum = tfdir.make<TH1F>("hmomentum", "Given photon momentum", 1000,  0.,  200.  );
      _hcos      = tfdir.make<TH1F>("hcos"     , "Given photon cos(#theta)", 1000,  -1.,  1.  );
      _hcosWt    = tfdir.make<TH1F>("hcosWt"   , "Given photon cos(#theta) weighted", 1000,  -1.,  1.  );
      _hr        = tfdir.make<TH1F>("hr"       , "Given photon radius from the DS axis", 1000,  0.,  2000.  );
      _hz        = tfdir.make<TH1F>("hz"       , "Given photon Z", 1000,  0., 15000.  );
      _hElecMom  = tfdir.make<TH1F>("hElecMom" , "Produced electron momentum", 2000,  0. , 200.);
      _hPosiMom  = tfdir.make<TH1F>("hPosiMom" , "Produced positron momentum", 2000,  0. , 200.);
      _hTotMom   = tfdir.make<TH1F>("hTotMom"   , "Produced total momentum", 2000.,  0. , 200.);
      _hTotMomWt = tfdir.make<TH1F>("hTotMomWt" , "Produced total momentum weighted", 2000.,  0. , 200.);
      _hMee      = tfdir.make<TH1F>("hMee"     , "M(e+e-) "           , 200,0.,200.);
      _hMeeVsE   = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"       , 200,0.,200.,200,0,200);
      _hMeeOverE = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "         , 200, 0.,1);
      _hy        = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);

    }

  }

  //================================================================
  void GammaConversionGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    
    //fields in the stop ntuple 
    double x, y, z, t, px, py, pz;
    char* mat = new char[40];
    double weight;
    double gen_energy;

    bool passed = false;
    CLHEP::Hep3Vector pos(0.,0.,0.);
    CLHEP::HepLorentzVector mome, momp;

    //Get the process information for photon conversions to sample from
    if(!photon_) {
      G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
      if(ptable) {
	photon_ = ptable->FindParticle(22);
	if(photon_) {
	  const G4ProcessManager* g4pm = photon_->GetProcessManager();
	  if(g4pm) {
	    G4VProcess* g4vp = g4pm->GetProcess("conv");
	    if(g4vp) {
	      G4GammaConversion* g4pprm = (G4GammaConversion*) g4vp;
	      if(g4pprm) {
		G4VEmModel* emmodel = g4pprm->EmModel(0); //BetheHeitler model, E_gamma < 80 GeV
		if(emmodel) {
		  g4bhm_ = (G4BetheHeitlerModel*) emmodel;
		  if(g4bhm_) {
		    printf("GammaConversionGun::produce : Successfully initialized Bethe Heitler Model\n");
		  } else
		    throw cet::exception("ERROR")
		      << "GammaConversionGun::produce : Failed to initialize Bethe Heitler Model\n";
		} else
		  throw cet::exception("ERROR")
		    << "GammaConversionGun::produce : Failed to initialize EM Model\n";
	      } else
		throw cet::exception("ERROR")
		  << "GammaConversionGun::produce : Failed to initialize GammaConversion process\n";
	    } else
	      throw cet::exception("ERROR")
		<< "GammaConversionGun::produce : Failed to initialize conversion process\n";
	  } else
	    throw cet::exception("ERROR")
	      << "GammaConversionGun::produce : Failed to initialize a process manager \n";
	} else
	  throw cet::exception("ERROR")
	    << "GammaConversionGun::produce : Failed to initialize a G4 photon \n";
      } else
	throw cet::exception("ERROR")
	  << "GammaConversionGun::produce : Failed to initialize a particle table\n";
    }
    
    //loop through generations until an event passes the generation cuts
    do {
      const auto& stop = stops_.fire();
      x  = stop.x;  y  = stop.y;  z  = stop.z;  t = stop.time;
      px = stop.px; py = stop.py; pz = stop.pz; 
      if(defaultMat_.size() == 0) sprintf(mat,"%s",stop.mat); 
      else                        sprintf(mat,"%s",defaultMat_.c_str()); 
      weight = stop.weight; gen_energy = stop.genEnergy;
      
      if(verbosityLevel_ > 2) {
	printf("Next Stop Attempt: (x,y,z,t) = (%.2f,%.2f,%.2f,%.2f), ",x,y,z,t);
	printf("(px,py,pz,gen_energy) = (%.2f,%.2f,%.2f,%.2f), ",px,py,pz,gen_energy);
	printf("material = %s\n", mat);
      }
      _hgencuts->Fill(1); //all generations
      passed = !(x < xMin_ || x > xMax_
		 || y < yMin_ || y > yMax_
		 || z < zMin_ || z > zMax_
		 || sqrt((x-xOffset_)*(x-xOffset_) + y*y) < rMin_ 
		 || sqrt((x-xOffset_)*(x-xOffset_) + y*y) > rMax_ );

      if(!passed) continue;
      _hgencuts->Fill(2); //passed spacial cut
      if(verbosityLevel_ > 2) 
	printf("Passed spacial cut\n");

      CLHEP::Hep3Vector mom(px, py, pz);
      if(testE_ > 0.) mom.setMag(testE_);
      
      passed = passed && (pMax_ < 0. || pMax_ > gen_energy);
      if(!passed) continue;
      _hgencuts->Fill(3); //passed maximum gen photon momentum cut
      if(verbosityLevel_ > 2) 
	printf("Passed gen energy cut\n");

      double cz = mom.cosTheta();
      passed = passed && cz >= czMin_ && cz <= czMax_;
      if(!passed) continue;
      _hgencuts->Fill(4); //Cos(theta) cut

      //can't make a daughter of pMin if energy below pmin + electron mass already
      //use slightly less than electron mass here to be safe
      double photonE = mom.mag();
      if(verbosityLevel_ > 2) 
	printf("Passed cos theta cut\nPhoton energy = %.2f\n", photonE);
      passed = passed && (photonE - 0.500) > pMin_;
      if(!passed) continue;
      _hgencuts->Fill(5); //passed initial minimum momentum cut
      if(verbosityLevel_ > 2) 
	printf("Passed min photon energy cut\n");

      //sample the spectrum
      std::vector<G4DynamicParticle*> *g4dpv = new std::vector<G4DynamicParticle*>();
      G4Material* material = findMaterialOrThrow(mat);
      G4MaterialCutsCouple* matcut = new G4MaterialCutsCouple(material);
      mom.setMag(1.); //just need momentum direction with the energy
      G4DynamicParticle* g4dp = new G4DynamicParticle(photon_, mom, photonE);
      g4bhm_->SampleSecondaries(g4dpv, matcut, g4dp,0.,0.);

      //get e+- pair
      mome = (*g4dpv)[0]->Get4Momentum();
      momp = (*g4dpv)[1]->Get4Momentum();

      //release memory
      delete g4dp;
      delete matcut;
      for(auto dp : *g4dpv)
	delete dp;
      delete g4dpv;

      passed = passed && (mome.vect().mag() > pMin_ || momp.vect().mag() > pMin_);
      if(!passed) continue;
      _hgencuts->Fill(6); //Final momentum cut
      if(verbosityLevel_ > 2) 
	printf("Passed min daughter energy cut\n");

      pos.setX(x); pos.setY(y); pos.setZ(z);

    } while(!passed);
    _hgencuts->Fill(10); //passing all cuts
    
                                                                                   //GenId = 44
    output->emplace_back(PDGCode::e_minus, GenId::gammaPairProduction, pos, mome, t);
    output->emplace_back(PDGCode::e_plus , GenId::gammaPairProduction, pos, momp, t);
    event.put(move(output));

    //add event weight and gen photon energy to the output
    std::unique_ptr<EventWeight> evtwt ( new EventWeight(weight) );
    event.put(move(evtwt));
    std::unique_ptr<mu2e::GenParticle> genenergy ( new mu2e::GenParticle(PDGCode::gamma, GenId::ExternalRMC, 
									 CLHEP::Hep3Vector(0.,0.,0.),
									 CLHEP::HepLorentzVector(0.,0.,gen_energy,gen_energy),0. ));
    event.put(move(genenergy),"photon"); 

    //release memory
    delete [] mat;
    if ( !doHistograms_ ) return;

    _hcos->Fill((mome+momp).vect().cosTheta());
    _hcosWt->Fill((mome+momp).vect().cosTheta(),weight);
    _hr->Fill(sqrt((pos.x()+3904.)*(pos.x()+3904.)+pos.y()*pos.y()));
    _hz->Fill(pos.z());
    _hElecMom ->Fill(mome.vect().mag());
    _hPosiMom ->Fill(momp.vect().mag());
    CLHEP::Hep3Vector p = mome.vect()+momp.vect();
    double momentum = p.mag();
    _hTotMom ->Fill(momentum);
    _hTotMomWt ->Fill(momentum, weight);

    double mee = (mome+momp).m();
    _hMee->Fill(mee);
    double energy = (mome+momp).e();
    _hMeeVsE->Fill(energy,mee);
    _hMeeOverE->Fill(mee/energy);

    double lepy = (mome.e()-momp.e())/energy;

    _hy->Fill(lepy);

    
    _hmomentum->Fill(energy);

  }


  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GammaConversionGun);
