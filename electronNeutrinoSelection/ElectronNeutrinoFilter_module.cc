////////////////////////////////////////////////////////////////////////
// Class:       ElectronNeutrinoFilter
// Plugin Type: filter (art v2_05_01)
// File:        ElectronNeutrinoFilter_module.cc
//
// Generated at Tue Jul 31 16:09:58 2018 by Stefano Soleti using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ElectronEventSelectionAlg.h"

#include <memory>

namespace lee {
  class ElectronNeutrinoFilter;
}


class lee::ElectronNeutrinoFilter : public art::EDFilter {
public:
  explicit ElectronNeutrinoFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ElectronNeutrinoFilter(ElectronNeutrinoFilter const &) = delete;
  ElectronNeutrinoFilter(ElectronNeutrinoFilter &&) = delete;
  ElectronNeutrinoFilter & operator = (ElectronNeutrinoFilter const &) = delete;
  ElectronNeutrinoFilter & operator = (ElectronNeutrinoFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  bool endSubRun(art::SubRun &sr) override;

private:
  lee::ElectronEventSelectionAlg fElectronEventSelectionAlg;
  std::ofstream _run_subrun_list_file;

};

lee::ElectronNeutrinoFilter::ElectronNeutrinoFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.

  _run_subrun_list_file.open("run_subrun_list_filter.txt", std::ofstream::out | std::ofstream::trunc);

  this->reconfigure(p);
}

bool lee::ElectronNeutrinoFilter::endSubRun(art::SubRun &sr)
{
  _run_subrun_list_file << sr.run() << " " << sr.subRun() << std::endl;
  return true;
}

bool lee::ElectronNeutrinoFilter::filter(art::Event & e)
{
  std::cout << "[ElectronNeutrinoFilter] "
            << "RUN " << e.run() << " SUBRUN " << e.subRun() << " EVENT " << e.id().event()
            << std::endl;

  bool passed = fElectronEventSelectionAlg.eventSelected(e);
  std::cout << "[ElectronNeutrinoFilter] Passing filter? " << passed << std::endl;
  return passed;
}

void lee::ElectronNeutrinoFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fElectronEventSelectionAlg.reconfigure(p.get<fhicl::ParameterSet>("ElectronSelectionAlg"));
}

DEFINE_ART_MODULE(lee::ElectronNeutrinoFilter)
