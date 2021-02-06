///////////////////////////////////////////////////////////////////////////////
#include "Print/inc/CrvPhotonsPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"

void mu2e::CrvPhotonsPrinter::Print(art::Event const& event, std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CrvPhotonsCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CrvPhotonsCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::CrvPhotonsPrinter::Print(const art::Handle<CrvPhotonsCollection>& handle, std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void mu2e::CrvPhotonsPrinter::Print(const art::ValidHandle<CrvPhotonsCollection>& handle, std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void mu2e::CrvPhotonsPrinter::Print(const CrvPhotonsCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CrvPhotonsCollection has " << coll.size() << " photons\n";
  if(verbose()!=1) return;
  
  os << std::setiosflags(std::ios::fixed | std::ios::right);

  PrintListHeader();

  for(const auto& obj: coll) {
    int bar_index = obj.first.asInt();

    os << " bar_index:" << std::setw(5) << bar_index << "\n";
    
    // loop over photons 

    for (int sipm=0; sipm<4; sipm++) {
      const std::vector<CrvPhotons::SinglePhoton>* list_of_photons =  &obj.second.GetPhotons(sipm);
      int nph = list_of_photons->size();
      os << " sipm, nphotons" << std::setw(5) << sipm << "  " << nph << "\n";
      for (int i=0; i<nph; i++) {
	const CrvPhotons::SinglePhoton* ph = &list_of_photons->at(i);
	const mu2e::SimParticle* simp = ph->_step->simParticle().get();

	os << "i, time, eDep, simID, pdgID : " << i << " " << ph->_time << " " << ph->_step->eDep() << " " << simp->pdgId() << "\n";
      }
    }
  }

}

void mu2e::CrvPhotonsPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::CrvPhotonsPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "   ind      id    time\n";

}

