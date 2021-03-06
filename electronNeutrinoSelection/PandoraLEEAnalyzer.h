/**
 * \class PandoraLEEAnalyzer
 *
 * \ingroup lee
 *
 * \brief Main analyzer class filling a ROOT TTree with information from the selected
 * electron neutrino candidate
 *
 *
 * \author Stefano Roberto Soleti <stefano.soleti@physics.ox.ac.uk>
 *
 * \version analyzer (art v2_05_00)
 *
 * \date 20/07/2018
 *
 * Contact: stefano.soleti@physics.ox.ac.uk
 *
 * Created on: Fri Jul  20 11:20:39 2018
 *
 */

/** \addtogroup lee
@{*/

#ifndef PANDORA_LEE_H
#define PANDORA_LEE_H

#include <math.h>
#include <fstream>
#include <algorithm>
#include <functional>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// uncomment the lines below as you use these objects

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "TEfficiency.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "ElectronEventSelectionAlg.h"
#include "PandoraInterfaceHelper.h"

#include "EnergyHelper.h"
#include "GeometryHelper.h"

#include "larcore/Geometry/WireGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "uboone/EventWeight/EventWeightTreeUtility.h"

namespace lee {

class PandoraLEEAnalyzer : public art::EDAnalyzer {
public:
  explicit PandoraLEEAnalyzer(fhicl::ParameterSet const &pset);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.
  virtual ~PandoraLEEAnalyzer();

  // Plugins should not be copied or assigned.
  PandoraLEEAnalyzer(PandoraLEEAnalyzer const &) = delete;
  PandoraLEEAnalyzer(PandoraLEEAnalyzer &&) = delete;
  PandoraLEEAnalyzer &operator=(PandoraLEEAnalyzer const &) = delete;
  PandoraLEEAnalyzer &operator=(PandoraLEEAnalyzer &&) = delete;

  /**
  * @brief Main analyzer method, runs for every event in the file
  *
  * @param[in] evt  current art::Event
  */
  void analyze(art::Event const &e) override;

  /**
  * @brief Method called at the end of each subrun, it stores the number of POT
  *
  * @param[in] sr  current art::SubRun
  */
  void endSubRun(const art::SubRun &sr);

  /**
  * @brief Assigns the values of the FHICL file
  *
  * @param[in] pset  set of FHICL parameters
  */
  void reconfigure(fhicl::ParameterSet const &pset) override;

  /**
  * @brief Clears filled variables
  */
  void clear();

  /**
  * @brief Return the longest reconstructed track
  *
  * @param[in] candidates  vector of neutrino candidates indexes
  * @param[in] evt         current art::Event
  *
  * @return            index of the chosen neutrino candidate
  */
  size_t choose_candidate(
      std::vector<size_t> &candidates,
      const art::Event &evt);

  /**
  * @brief Return the longest reconstructed track
  *
  * @param[in] tracks  vector of reconstructed tracks
  *
  * @return            longest reconstructed track object
  */
  art::Ptr<recob::Track> get_longest_track(
      std::vector<art::Ptr<recob::Track>> &tracks);

  /*
  * @brief Return the most energetic reconstructed shower
  *
  * @param[in] showers  vector of reconstructed showers
  *
  * @return             most energetic reconstructed shower object
  */
  art::Ptr<recob::Shower> get_most_energetic_shower(
      std::vector<art::Ptr<recob::Shower>> &showers);

  /**
  * @brief Determines if a PFParticle is matched with a MCParticle coming from
  * a neutrino interaction or a cosmic ray
  *
  * @param[in]  evt               current art Event
  * @param[out] neutrino_pdg      vector of PDG codes for neutrino-matched PFParticles
  * @param[out] neutrino_process  vector of interaction processes for neutrino-matched PFParticles
  * @param[out] neutrino_energy   vector of true energy of neutrino-matched MCParticles
  * @param[out] neutrino_pf       vector of neutrino-matched PFParticles
  * @param[out] cosmic_pdg        vector of PDG codes for cosmic-matched PFParticles
  * @param[out] cosmic_process    vector of interaction processes for cosmic-matched PFParticles
  * @param[out] cosmic_energy     vector of true energy of cosmic-matched MCParticles
  * @param[out] cosmic_pf         vector of cosmic-matched PFParticles
  */
  void categorizePFParticles(
      art::Event const &evt,

      std::vector<int> &neutrino_pdg,
      std::vector<std::string> &neutrino_process,
      std::vector<double> &neutrino_energy,
      std::vector<art::Ptr<recob::PFParticle>> &neutrino_pf,

      std::vector<int> &cosmic_pdg,
      std::vector<std::string> &cosmic_process,
      std::vector<double> &cosmic_energy,
      std::vector<art::Ptr<recob::PFParticle>> &cosmic_pf);

  void fillTrackFields(size_t pf_id,
                       recob::PFParticle const *pfparticle,
                       art::FindManyP<recob::Cluster> *clusters_per_pfpart,
                       art::FindManyP<recob::SpacePoint> *spcpnts_per_pfpart,
                       art::FindOneP<recob::Track> *track_per_pfpart,
                       art::FindManyP<anab::Calorimetry> *calos_per_track,
                       art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> *mcps_per_hit,
                       art::FindManyP<recob::Hit> *hits_per_spcpnts,
                       art::FindManyP<anab::ParticleID> *pid_per_track,
                       art::FindManyP<recob::Hit> *hits_per_cluster);

private:
  std::string m_hitfinderLabel;
  std::string _geantModuleLabel = "largeant";
  std::string m_pfp_producer;
  std::string m_pid_producer;
  std::string m_calorimetry_producer;

  std::string m_spacepointLabel;
  std::string _mctruthLabel = "generator";
  std::string m_hitmatching_producer;

  lee::ElectronEventSelectionAlg fElectronEventSelectionAlg;
  EnergyHelper energyHelper;
  GeometryHelper geoHelper;
  PandoraInterfaceHelper pandoraHelper;
  uboone::EWTreeUtil ewutil;
  trkf::TrackMomentumCalculator _trkmom;

  float _lee_bins[12] = {200, 300, 375, 475, 550, 675, 800, 950, 1100, 1300, 1500, 3000};
  float _lee_scaling[13] = {0, 3.88549, 3.05421, 1.59615, 0.383725, 0, 0, 0, 0, 0, 0, 0, 0};

  TH1F *_h_lee_scaling = new TH1F("h_lee_scaling", "", 11, _lee_bins);

  TFile *myTFile;
  TTree *myTTree;
  TTree *mySBNTTree;
  TTree *myPOTTTree;

  double _reco_nu_energy;
  bool m_showersAsTracks;
  bool m_isLEE;
  int _interaction_type;

  double m_fidvolXstart;
  double m_fidvolXend;

  double m_fidvolYstart;
  double m_fidvolYend;

  double m_fidvolZstart;
  double m_fidvolZend;
  bool m_useParticleID;
  bool m_isData;
  bool m_isCosmicInTime;
  bool m_printDebug;
  bool m_isOverlaidSample;
  bool m_save_flux_info;
  const int k_cosmic = 1;
  const int k_nu_e = 2;
  const int k_nu_mu = 3;
  const int k_nc = 4;
  const int k_dirt = 5;
  const int k_data = 6;
  const int k_other = 0;
  const int k_mixed = 7;
  const int k_nu_e_other = 8;
  std::vector<double> _energy;
  int _true_nu_is_fiducial;
  double _nu_energy;

  int _n_tracks;
  int _n_showers;
  int _n_showers_as_tracks;
  double _vx;
  double _vy;
  double _vz;
  double _interaction_length;
  double _true_vx;
  double _true_vy;
  double _true_vz;

  double _true_vx_sce;
  double _true_vy_sce;
  double _true_vz_sce;
  std::vector< std::vector<double> > _genie_weights;
  std::vector< std::string > _genie_names;

  std::vector< std::vector<double> > _flux_weights;
  std::vector< std::string > _flux_names;

  std::vector<double> _true_shower_x_sce;
  std::vector<double> _true_shower_y_sce;
  std::vector<double> _true_shower_z_sce;
  std::vector<int> _true_shower_pdg;
  std::vector<double> _true_shower_depE;

  int _nu_matched_tracks;
  int _nu_matched_showers;

  int _nu_pdg;
  int _ccnc;
  int _category;
  int _run;
  int _subrun;
  int _event;
  int _n_candidates;
  int _n_true_nu;
  int _run_sr;
  int _subrun_sr;
  int _n_matched;
  double _pot;
  int _event_passed;
  int _numu_passed;
  int _numu_cuts;
  double _distance;
  double _cosmic_fraction;

  std::vector<int> _flash_passed;
  std::vector<int> _track_passed;
  std::vector<int> _shower_passed;
  std::vector<int> _primary_indexes;
  std::vector<int> _number_tracks;
  std::vector<int> _number_showers;

  std::vector<int> _matched_showers;
  std::vector<int> _matched_tracks;

  std::vector<std::string> _matched_tracks_process;
  std::vector<double> _matched_tracks_energy;

  std::vector<std::string> _matched_showers_process;
  std::vector<double> _matched_showers_energy;

  std::map<std::string, std::vector<double>> _weights;

  int _n_primaries;
  int _chosen_candidate;

  float _leeweight;
  double _bnbweight;

  std::vector<std::vector<double>> _shower_dQdx_hits;
  std::vector<std::vector<double>> _shower_dEdx_hits;

  std::vector<std::vector<double>> _shower_dQdx;
  std::vector<std::vector<double>> _shower_dEdx;
  std::vector<std::vector<double>> _shower_dQdx_cali;
  std::vector<std::vector<double>> _shower_dEdx_cali;
  std::vector<std::vector<double>> _shower_pitches;
  std::vector<std::vector<int>> _shower_dQdx_hits_in_the_box;

  std::vector<std::vector<double>> _track_dQdx_hits;
  std::vector<std::vector<double>> _track_dEdx_hits;

  std::vector<std::vector<double>> _track_dQdx;
  std::vector<std::vector<double>> _track_dQdx_cali;

  std::vector<std::vector<double>> _track_dEdx;
  std::vector<double> _track_energy_length;

  std::vector<size_t> _nu_track_ids;
  std::vector<size_t> _nu_shower_as_track_ids;

  std::vector<size_t> _nu_shower_ids;

  std::vector< std::vector<size_t> > _nu_track_daughters;
  std::vector< std::vector<size_t> > _nu_shower_daughters;

  std::vector<double> _shower_open_angle;
  std::vector<double> _shower_length;
  std::vector<double> _shower_dir_x;
  std::vector<double> _shower_dir_y;
  std::vector<double> _shower_dir_z;

  std::vector<double> _shower_start_x;
  std::vector<double> _shower_start_y;
  std::vector<double> _shower_start_z;

  std::vector<double> _shower_theta;
  std::vector<double> _shower_phi;

  std::vector<std::vector<double>> _shower_energy;
  std::vector<std::vector<double>> _shower_energy_cali;
  std::vector<std::vector<double>> _track_energy_cali;

  std::vector<double> _track_dir_x;
  std::vector<double> _track_dir_y;
  std::vector<double> _track_dir_z;
  std::vector<int> _track_is_fiducial;
  std::vector<int> _shower_is_fiducial;

  std::vector<double> _track_res_mean;
  std::vector<double> _track_res_std;

  std::vector<double> _shower_res_mean;
  std::vector<double> _shower_res_std;

  std::vector<double> _track_start_x;
  std::vector<double> _track_start_y;
  std::vector<double> _track_start_z;

  std::vector<double> _track_end_x;
  std::vector<double> _track_end_y;
  std::vector<double> _track_end_z;

  std::vector<double> _track_theta;
  std::vector<double> _track_phi;

  std::vector<double> _track_length;
  std::vector<double> _track_id;

  std::vector<double> _track_bragg_p;
  std::vector<double> _track_bragg_mu;
  std::vector<double> _track_bragg_mip;

  std::vector<double> _track_pidchi;
  std::vector<double> _track_pidchipr;
  std::vector<double> _track_pidchika;
  std::vector<double> _track_pidchipi;
  std::vector<double> _track_pidchimu;
  std::vector<double> _track_pida;

  std::vector<double> _track_energy_dedx;
  std::vector<std::vector<double>> _track_energy_hits;

  std::vector<int> _nu_daughters_pdg;
  std::vector<double> _nu_daughters_E;

  std::vector<double> _nu_daughters_px;
  std::vector<double> _nu_daughters_py;
  std::vector<double> _nu_daughters_pz;

  std::vector<double> _nu_daughters_vx;
  std::vector<double> _nu_daughters_vy;
  std::vector<double> _nu_daughters_vz;

  std::vector<double> _nu_daughters_endx;
  std::vector<double> _nu_daughters_endy;
  std::vector<double> _nu_daughters_endz;

  std::vector<double> _flash_PE;
  std::vector<double> _flash_time;

  double _TPC_x;
  double _flash_x;

  std::vector<std::vector<double>> _shower_pca;
  std::vector<std::vector<double>> _track_pca;

  std::vector<std::vector<int>> _shower_nhits;
  std::vector<std::vector<int>> _track_nhits;

};

}

#endif // PANDORA_LEE_H
