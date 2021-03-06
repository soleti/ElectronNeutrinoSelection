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
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "TTree.h"
#include "TH1F.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "uboone/EventWeight/EventWeightTreeUtility.h"
#include "uboone/EventWeight/MCEventWeight.h"

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
    void respondToOpenInputFile(art::FileBlock const &fb) override;
    bool endSubRun(art::SubRun &sr) override;
    void clear();

  private:
    TTree *myPOTTTree;
    TTree *myTTree;

    lee::ElectronEventSelectionAlg fElectronEventSelectionAlg;
    std::ofstream _run_subrun_list_file;
    GeometryHelper geoHelper;
    PandoraInterfaceHelper pandoraHelper;

    float _lee_bins[12] = {200, 300, 375, 475, 550, 675, 800, 950, 1100, 1300, 1500, 3000};
    float _lee_scaling[13] = {0, 3.88549, 3.05421, 1.59615, 0.383725, 0, 0, 0, 0, 0, 0, 0, 0};

    TH1F *_h_lee_scaling = new TH1F("h_lee_scaling", "", 11, _lee_bins);
    int _selection_result;

    std::vector< std::vector<double> > _genie_weights;
    std::vector< std::string > _genie_names;

    std::vector< std::vector<double> > _flux_weights;
    std::vector< std::string > _flux_names;

    double _lee_weight;
    double _bnbweight;
    double _sum_pot;
    bool _passed;
    bool m_isOverlaidSample;
    bool m_isData;
    unsigned int _run_sr;
    unsigned int _subrun_sr;
    unsigned int _n_true_nu;
    int _ccnc;
    double _qsqr;
    double _pt;
    double _theta;
    double _w;
    int _nu_pdg;
    int _interaction_type;
    unsigned int _n_primaries;
    bool _true_nu_is_fiducial;
    double _nu_energy;
    double _pot;
    double _true_vx;
    double _true_vy;
    double _true_vz;
    double _true_vx_sce;
    double _true_vy_sce;
    double _true_vz_sce;

    std::vector<double> _nu_daughters_E;
    std::vector<int> _nu_daughters_pdg;
    std::vector<double> _flash_time;
    std::vector<double> _flash_pe;

    std::vector < std::vector<double> > _nu_daughters_p;
    std::vector < std::vector<double> > _nu_daughters_start_v;
    std::vector < std::vector<double> > _nu_daughters_end_v;

    std::string _mctruthLabel = "generator";
    std::string _mcparticleLabel = "largeant";
};
