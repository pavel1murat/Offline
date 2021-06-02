// Author: Shihua Huang, 2021
// Modified from Andrei Gaponenko, 2013
// This work used the V-A model of muon decay spectrum having radiative corrections to the one-loop order. 
//See reference:
//T Kinoshita and A Sirlin, "Radiative corrections to Fermi interactions", 1959

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"  
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/ConversionSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

#include "TH1.h"
#include "TH2.h"

using CLHEP::Hep3Vector;

namespace mu2e {

  //================================================================
  class StoppedParticleReactionGunMuPlus : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    PDGCode::type       pdgId_;
    double              mass_;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM};
    SpectrumVar       spectrumVariable_;
   
    CLHEP::Hep3Vector spin_{0,0,-1}; //default spin points to negative z
   
    BinnedSpectrum    spectrum_;
    GenId             genId_;
    int               verbosityLevel_;
    long int          seed_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral randSpectrum_;
    CLHEP::RandFlat    flat_;
    RandomUnitSphere   randomUnitSphere_;

    // CLHEP::RandExponential randTime_;
    double                                fixedTime_;
    double                                mmu_{105.658};//muon mass
    double                                me_{0.511};//electron mass
    double                                Eemax_{52.8302};//maximum electron energy
    double                                alpha_{0.00729735};//fine structure constant
     RootTreeSampler<IO::StoppedParticleF> stops_; //should this be used?
    // RootTreeSampler<IO::StoppedParticleTauNormF> stops_;
//-----------------------------------------------------------------------------
// histogramming
//-----------------------------------------------------------------------------
    bool    doHistograms_;
    TH1F*   _hEnergy[2];
    TH1F*   _hMom[2];//to be filled
    TH1F*   _hPdgId;
    TH1F*   _hGenId;
    TH1F*   _hTime;
    TH1F*   _hZ[2];
    TH2F*   _hRVsZ;
  
  private:
    static SpectrumVar    parseSpectrumVar(const std::string& name);
    double                generateEnergy();
    double                generateAngle(double const& Energy);
    double                Michel_non_polar(double const& x0);
    double                Michel_polar(double const& x0);
  public:
    explicit StoppedParticleReactionGunMuPlus(const fhicl::ParameterSet& pset);

    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedParticleReactionGunMuPlus::StoppedParticleReactionGunMuPlus(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    , psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , pdgId_(PDGCode::type(psphys_.get<int>("pdgId")))
    , mass_(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId_).ref().mass().value())
    , spectrumVariable_(parseSpectrumVar(psphys_.get<std::string>("spectrumVariable")))
    , spectrum_(BinnedSpectrum(psphys_))
    , genId_(GenId::findByName(psphys_.get<std::string>("genId")))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , seed_(pset.get<long int>("seed",art::ServiceHandle<SeedService>()->getSeed()))
    , eng_(createEngine(seed_))
    , randSpectrum_(eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , flat_(eng_)
    , randomUnitSphere_(eng_)
    , fixedTime_          (pset.get<double>("fixedTime"   ))
    , stops_(eng_, pset.get<fhicl::ParameterSet>("muonStops"))//needs to add momentum
    , doHistograms_       (pset.get<bool>  ("doHistograms",true))
  {
    produces<mu2e::GenParticleCollection>();

    if(genId_ == GenId::enum_type::unknown) {
      throw cet::exception("BADCONFIG")<<"StoppedParticleReactionGunMuPlus: unknown genId "
                                       <<psphys_.get<std::string>("genId", "StoppedParticleReactionGunMuPlus")
                                       <<"\n";
    }

    if (verbosityLevel_ > 0) {
      std::cout<<"StoppedParticleReactionGunMuPlus: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"StoppedParticleReactionGunMuPlus: using seed  = " << seed_  << std::endl;
      std::cout<<"StoppedParticleReactionGunMuPlus: using GenId = " << genId_ << std::endl;

      std::cout<<"StoppedParticleReactionGunMuPlus: producing particle "<< pdgId_ << ", mass = "<< mass_ << std::endl;

      std::cout <<"StoppedParticleReactionGunMuPlus: spectrum shape = "
		<<psphys_.get<std::string>("spectrumShape") << std::endl;
      if (psphys_.get<std::string>("spectrumShape")  == "tabulated")
	std::cout << " Spectrum file = "
		  << psphys_.get<std::string>("spectrumFileName")
		  << std::endl;
    }
    if (verbosityLevel_ > 1){
      std::cout <<"StoppedParticleReactionGunMuPlus: spectrum: " << std::endl;
      spectrum_.print();
    }

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      //      art::TFileDirectory tfdir = tfs->mkdir( "StoppedParticleReactionGun");
      _hEnergy[0] = tfs->make<TH1F>("hEnergy" , "Energy[0]"   , 2400,   0.0,   120);
      _hEnergy[1] = tfs->make<TH1F>("hEnergy2", "Energy[1]"   , 2400,   0.0,  1200);
      _hMom[0]       = tfs->make<TH1F>("hMom"    , "Momentum"    ,  500,   0.0,  1000);
      _hMom[1]       = tfs->make<TH1F>("hMomZ"    , "Momentum"    ,  500,   0.0,  1000); //want to know z momentum as well
      _hGenId     = tfs->make<TH1F>("hGenId"  , "Generator ID",  100,   0.0,   100);
      _hPdgId     = tfs->make<TH1F>("hPdgId"  , "PDG ID"      ,  500,  -250,   250);
      _hTime      = tfs->make<TH1F>("hTime"   , "Time"        ,  400,   0.0,  2000);
      _hZ[0]      = tfs->make<TH1F>("hZ"      , "Z"           ,  500,  5400,  6400);
      _hZ[1]      = tfs->make<TH1F>("hZ2"     , "Z[1]"        , 2500,     0, 25000);
      _hRVsZ      = tfs->make<TH2F>("hRvsZ"   , "R vs Z]"     , 2500,     0, 25000, 100,0,1000);

    }
  }


  //================================================================
  StoppedParticleReactionGunMuPlus::SpectrumVar StoppedParticleReactionGunMuPlus::parseSpectrumVar(const std::string& name) {
    if (name == "totalEnergy"  )  return TOTAL_ENERGY;
    if (name == "kineticEnergy")  return KINETIC_ENERY;
    if (name == "momentum"     )  return MOMENTUM;
    throw cet::exception("BADCONFIG")<<"StoppedParticleReactionGunMuPlus: unknown spectrum variable "<<name<<"\n";
  }


  //================================================================
  void StoppedParticleReactionGunMuPlus::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    const double energy = generateEnergy();//positron energy
    const double p      = sqrt(energy*energy-mass_*mass_);//positron momentum
    const double cosTheta = generateAngle(energy);//positron emission angle w.r.t the muon spin polarization

    double time = stop.t;

    if (fixedTime_ >= 0) time = fixedTime_;
    
    double phi = CLHEP::twopi*flat_.fire();
    // CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(p);
    CLHEP::Hep3Vector p3 = polar3Vector (p, cosTheta, phi);//momentum direction with respect to spin polarization
    double px;double py;double pz;
    if(std::sqrt(1-spin_.z()*spin_.z())>0){
    px=p3.x()*spin_.x()*spin_.z()/std::sqrt(1-spin_.z()*spin_.z())-p3.y()*spin_.y()/std::sqrt(1-spin_.z()*spin_.z())+p3.z()*spin_.x()/(1-spin_.z()*spin_.z());
    py=p3.x()*spin_.x()*spin_.z()/std::sqrt(1-spin_.z()*spin_.z())+p3.y()*spin_.x()/std::sqrt(1-spin_.z()*spin_.z())+p3.z()*spin_.y()/(1-spin_.z()*spin_.z());
    pz=-p3.x()/std::sqrt(1-spin_.z()*spin_.z())+p3.z()*spin_.z();
    }
    else if(spin_.z()==-1){
      px=-p3.x();py=p3.y();pz=-p3.z();
       }
    else{
      px=p3.x();py=p3.y();pz=p3.z();
    }
    CLHEP::Hep3Vector p3_lab={px,py,pz};
    CLHEP::HepLorentzVector fourmom(p3_lab, energy);
    output->emplace_back(pdgId_,
                         genId_,
                         pos,
                         fourmom,
                         time);

    event.put(std::move(output));
//-----------------------------------------------------------------------------
// if requested, fill histograms. Currently, the only one
//-----------------------------------------------------------------------------
    if (doHistograms_) {
      _hGenId->Fill(genId_.id());
      _hPdgId->Fill(pdgId_);
      _hEnergy[0]->Fill(energy);
      _hEnergy[1]->Fill(energy);
      _hMom[0]->Fill(p);
      _hMom[1]->Fill(p3_lab.z());
      _hTime->Fill(stop.t);
      _hZ[0]->Fill(pos.z());
      _hZ[1]->Fill(pos.z());

      float dx = stop.x+3904;
      float r  = sqrt(dx*dx+stop.y*stop.y);

      _hRVsZ->Fill(pos.z(),r);
    }
  }

//-----------------------------------------------------------------------------
// generate (pseudo-)random particle energy 
// the spectrum itself doesn't know whether is stored momentum, kinetic or full 
// energy
//-----------------------------------------------------------------------------
  double StoppedParticleReactionGunMuPlus::generateEnergy() {
    double res = spectrum_.sample(randSpectrum_.fire());

    if (res < 0.0) {
      throw cet::exception("BADE")<<"StoppedParticleReactionGunMuPlus: negative energy "<< res <<"\n";
    }

    switch(spectrumVariable_) {
    case TOTAL_ENERGY  : break;
    case KINETIC_ENERY : res += mass_; break;
    case MOMENTUM      : res = sqrt(res*res+mass_*mass_); break;
    }
    return res;
  }
   
 //-----------------------------------------------------------------------------
// generate particle emission angle w.r.t spin polarization from Michel decay
//
//-----------------------------------------------------------------------------
  double StoppedParticleReactionGunMuPlus::generateAngle(double const& Energy) {
    double x0=Energy/Eemax_;//e+ energy divided by maximum e+ energy

    double cosTheta=0;//costheta w.r.t to the polarization

    double cumulativeP = flat_.fire();//cumulative probability(P) for a certain costheta A: Integral dP/dA in [-1,A]=cumulativeP;

    double C1=Michel_non_polar(x0);//the non-polarization related differential cross-section knowing the positron four-vector

    double C2=Michel_polar(x0);//the non-polarization related differential cross-section knowing the positron four-vector

    cosTheta =(-1+std::sqrt(1-2*C2/C1*(1-C2/(2*C1)-2*cumulativeP)))/(C2/C1);
    return cosTheta;
  }

//-----------------------------------------------------------------------------
// generate the differential cross-section(un-normalized) related to the polarization
// energy
//-----------------------------------------------------------------------------
  double StoppedParticleReactionGunMuPlus::Michel_polar(double const& x0) {
    
    double omega=std::log(mmu_/me_);//logarithm factor
    double R=0;
  
   if (x0 < me_/Eemax_) {
      throw cet::exception("BADE")<<"StoppedParticleReactionGunMuPlus: negative energy "<< x0 <<"\n";
    }

    if(x0>0.94){
      R=2*(-252.4354+796.87534*x0-837.0617*x0*x0+294.2622*x0*x0*x0)-(CLHEP::pi*CLHEP::pi)/3-2+omega*(3/2+2*std::log((1-x0)/x0))-(2*std::log(x0)-1)*std::log(x0)+(3*std::log(x0)-1-1/x0)*std::log(1-x0);
  }
    else {
       R=2*(-0.01309+1.21057*x0-0.49195*x0*x0+0.86247*x0*x0*x0)-(CLHEP::pi*CLHEP::pi)/3-2+omega*(3/2+2*std::log((1-x0)/x0))-(2*std::log(x0)-1)*std::log(x0)+(3*std::log(x0)-1-1/x0)*std::log(1-x0);
    }

    double g=0;
    g=(2-4*x0)*R+(2-6*x0)*std::log(x0)-(1-x0)/(3*x0*x0)*((1+x0+34*x0*x0)*(omega+std::log(x0))+3-7*x0-32*x0*x0+4*(1-x0)*(1-x0)/x0*std::log(1-x0));
    
    double Beta=std::sqrt(x0*x0-(me_/Eemax_)*(me_/Eemax_))/x0;//beta factor of emitted positron
    double C2=-Beta*(1-2*x0+me_*me_/(mmu_*Eemax_)+alpha_/CLHEP::twopi*g);//positron having the extra minus sign

    return C2;
  //================================================================
}
//-----------------------------------------------------------------------------
// generate the differential cross-section(un-normalized) not related to the polarization
// energy
//-----------------------------------------------------------------------------
  double StoppedParticleReactionGunMuPlus::Michel_non_polar(double const& x0) {
    
    double omega=std::log(mmu_/me_);//logarithm factor
    double R=0;
  
   if (x0 < me_/Eemax_) {
      throw cet::exception("BADE")<<"StoppedParticleReactionGunMuPlus: negative energy "<< x0 <<"\n";
    }

    if(x0>0.94){
      R=2*(-252.4354+796.87534*x0-837.0617*x0*x0+294.2622*x0*x0*x0)-(CLHEP::pi*CLHEP::pi)/3-2+omega*(3/2+2*std::log((1-x0)/x0))-(2*std::log(x0)-1)*std::log(x0)+(3*std::log(x0)-1-1/x0)*std::log(1-x0);
  }
    else {
       R=2*(-0.01309+1.21057*x0-0.49195*x0*x0+0.86247*x0*x0*x0)-(CLHEP::pi*CLHEP::pi)/3-2+omega*(3/2+2*std::log((1-x0)/x0))-(2*std::log(x0)-1)*std::log(x0)+(3*std::log(x0)-1-1/x0)*std::log(1-x0);
    }

    double f=0;
    f=(6-4*x0)*R+(6-6*x0)*std::log(x0)+(1-x0)/(3*x0*x0)*((5+17*x0-34*x0*x0)*(omega+std::log(x0))-22*x0+34*x0*x0);
    
    double C1=3-2*x0-me_*me_/(mmu_*Eemax_)+alpha_/CLHEP::twopi*f;

    return C1;
  //================================================================
}
 // namespace mu2e
}
DEFINE_ART_MODULE(mu2e::StoppedParticleReactionGunMuPlus);
