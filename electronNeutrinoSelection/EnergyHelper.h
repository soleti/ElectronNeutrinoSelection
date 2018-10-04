/**
 * \class EnergyHelper
 *
 * \ingroup lee
 *
 * \brief Helper class with methods accessing calorimetric information
 * (energy and dE/dx)
 *
 *
 * \author Stefano Roberto Soleti <stefano.soleti@physics.ox.ac.uk>
 *
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

#ifndef ENERGYHELPER_H
#define ENERGYHELPER_H

#include "HelperBase.h"

#include "GeometryHelper.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "TPrincipal.h"
// #include "uboone/UBXSec/Algorithms/TrackQuality.h"
// #include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// #include "uboone/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

namespace lee {

class EnergyHelper : public HelperBase {
public:
  explicit EnergyHelper();
  ~EnergyHelper() = default;

  /**
    * @brief Configure all of the parameters of this class
    *
    * @param p fcl parameter set
    */
  void reconfigure(fhicl::ParameterSet const &pset);

  /**
   * @brief      Measure the dQdx of a shower
   *
   * @param[in]  shower_obj          Pointer to the shower object
   * @param[in]  clusters            Pointer to the cluster associated tot the shower
   * @param[in]  hits_per_cluster    Pointer to the hits-cluster association
   * @param[out] dqdx                Vector of the dQ/dx median values per plane
   * @param[out] dqdx_hits           Vector of the dQ/dx hits values per plane
   * @param[out] pitches             Vector of the pitches per plane
   */
  void dQdx(const recob::Shower *shower_obj,
            std::vector<art::Ptr<recob::Cluster>> *clusters,
            art::FindManyP<recob::Hit> *hits_per_cluster,
            art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> *mcps_per_hit,
            std::vector<double> &dqdx,
            std::vector<std::vector<double>> &dqdx_hits,
            std::vector<double> &pitches);

  /**
   * @brief      Return the value of a specific ParticleID algorithm for a single track
   *
   * @param[in]  pids            Pointer to the vector of ParticleID objects
   * @param[in]  trackID         TrackID value
   * @param[in]  AlgName         Name of the algorithm (e.g. Bragg peak, PIDA, etc.)
   * @param[in]  VariableType    Type of variable (e.g. likelihood)
   * @param[in]  TrackDirection  Assumed direction of the track (anab::kBackward or anab::kForward)
   * @param[in]  pdgCode         Assumed PDG code of the track
   * 
   * @return  Value of VariableType for the ParticleID AlgName algorithm,
   *          given its hypothetical PDG code
   */
  double PID(art::Ptr<anab::ParticleID> selected_pid,
                           std::string AlgName,
                           anab::kVariableType VariableType,
                           anab::kTrackDir TrackDirection,
                           int pdgCode);

  /**
   * @brief      Convert dQ/dx vector into dE/dx vector (in MeV)
   *
   * @param[out] dedx          Address of the dE/dx vector
   * @param[in]  dqdx          dQ/dx vector
   */
  void dEdx_from_dQdx(std::vector<double> &dedx,
                      std::vector<double> dqdx);

  /**
   * @brief      Principal Component Analysis of reconstructed clusters
   *
   * @param[in]  clusters          Pointer to the vector of reconstructed clusters
   * @param[in]  hits_pre_cluster  Pointer to the hits-cluster association
   * @param[out] pca_planes        Address of the two-dimensional vector storing the PCA eigenvalues per plane
   */
  void PCA(std::vector<art::Ptr<recob::Cluster>> *clusters,
           art::FindManyP<recob::Hit> *hits_per_cluster,
           std::vector<std::vector<double>> &pca_planes);

  /**
   * @brief      Calibration value for the energy of a reconstructed object
   *
   * @param[in]  spcpnts           Pointer to the vector of reconstructed space-points
   * @param[in]  hits_per_spcpnts  Pointer to the hits-spacepoints association
   * @param[out] cali_corr         Address of the two-dimensional vector storing calibration values
   */
  void get_cali(std::vector<art::Ptr<recob::SpacePoint>> *spcpnts,
                art::FindManyP<recob::Hit> *hits_per_spcpnts,
                std::vector<double> &cali_corr);

  /**
   * @brief      Measure calorimetric energy for a reconstructed object
   *
   * @param[in]  clusters          Pointer to the vector of reconstructed clusters
   * @param[in]  hits_per_cluster  Pointer to the hits-cluster association
   * @param[out] nHits             Address of the vector of the number of hits per plane
   * @param[out] pfenergy          Address of the vector of reconstructed energy per plane
   */
  void energy_from_hits(std::vector<art::Ptr<recob::Cluster>> *clusters,
                        art::FindManyP<recob::Hit> *hits_per_cluster,
                        std::vector<int> &nHits,
                        std::vector<double> &pfenergy);

  /**
   * @brief      Measure the spatial residuals of the hits in a reconstructed cluster along its direction
   *
   * @param[in]  clusters          Pointer to the vector of reconstructed clusters
   * @param[in]  hits_per_cluster  Pointer to the hits-cluster association
   * @param[out] mean_v            Address of the mean value
   * @param[out] std_v             Address of the standard deviation value
   */
  void cluster_residuals(std::vector<art::Ptr<recob::Cluster>> *clusters,
                         art::FindManyP<recob::Hit> *hits_per_cluster,
                         double &mean_v,
                         double &std_v);

  /**
   * @brief      Measure the dQ/dx and dE/dx of a track using the anab::Calorimetry information
   *
   * @param[in]  calos             Pointer to the vector of the calorimetry objects
   * @param[out] dqdx              Address of the vector of the dQ/dx values per plane
   * @param[out] dedx              Address of the vector of the dE/dx values per plane
   */
  void track_dQdx(std::vector<art::Ptr<anab::Calorimetry>> *calos,
                  std::vector<double> &dqdx,
                  std::vector<double> &dedx);

  /**
   * @brief      Calibration value for the dQ/dx of a reconstructed shower
   *
   * @param[in]  shower_obj        Pointer to the reconstructed shower
   * @param[out] dqdx_cali         Address of the vector of with the dQ/dx calibration values
   */
  void dQdx_cali(const recob::Shower *shower_obj,
                 std::vector<double> &dqdx_cali);

  bool is_hit_data(art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> *mcps_per_hit,
                                 size_t hit_key);

private:
  std::vector<double> _data_gain = {236.41, 228.83, 242.72}; // DocDB 14754
  std::vector<double> _mc_gain = {193.05, 196.85, 196.85};  // Plane 0, plane 1, plane 2
  std::vector<double> _gain;
  const lariov::TPCEnergyCalibProvider &_energy_calib_provider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
  const detinfo::DetectorProperties *_detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double _drift = _detprop->DriftVelocity() * 1e-3;
  double _readout_window = 4.8;
  double _from_tick_to_ns = _readout_window / _detprop->ReadOutWindowSize() * 1e6;
  double _wire_spacing = 0.3;
  double _work_function = 23 / 1e6;
  double _betap;
  double _alpha;
  double _recombination_factor;
  double _dQdx_rectangle_length;
  double _dQdx_rectangle_width;
  bool m_isOverlaidSample;
  GeometryHelper geo_helper;
};
} // namespace lee

#endif
