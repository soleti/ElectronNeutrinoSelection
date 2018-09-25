#ifndef ELECTRON_EVENT_SELECTION_ALG_CXX
#define ELECTRON_EVENT_SELECTION_ALG_CXX

#include "ElectronEventSelectionAlg.h"

namespace lee
{

void ElectronEventSelectionAlg::clear()
{
  _n_neutrino_candidates = 0.0;
  _TPC_x = std::numeric_limits<double>::lowest();
  _flash_x = std::numeric_limits<double>::lowest();

  _primary_indexes.clear();
  _neutrino_candidate_passed.clear();
  _op_flash_indexes.clear();
  _neutrino_vertex.clear();
  _n_showers.clear();
  _n_tracks.clear();
  _n_showers_as_tracks.clear();
  _pfp_id_showers_from_primary.clear();
  _pfp_id_tracks_from_primary.clear();
  _pfp_id_showers_as_tracks_from_primary.clear();
  _flash_PE.clear();
  _flash_time.clear();
}

TVector3 ElectronEventSelectionAlg::spaceChargeTrueToReco(const TVector3 &xyz)
{

  TVector3 correctedPoint(xyz);
  try
  {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    if (SCE->GetPosOffsets(xyz.X(), xyz.Y(), xyz.Z()).size() == 3)
    {
      correctedPoint.SetX(xyz.X() - SCE->GetPosOffsets(xyz.X(), xyz.Y(), xyz.Z())[0] + 0.7);
      correctedPoint.SetX(xyz.Y() + SCE->GetPosOffsets(xyz.X(), xyz.Y(), xyz.Z())[1]);
      correctedPoint.SetX(xyz.Z() + SCE->GetPosOffsets(xyz.X(), xyz.Y(), xyz.Z())[2]);
    }
    else
    {
      std::cout << "[PandoraLEE] "
                << "Space Charge service offset size not 3" << std::endl;
    }
  }
  catch (...)
  {
    std::cout << "[PandoraLEE] "
              << "Space Charge service error" << std::endl;
  }

  return correctedPoint;
}

void ElectronEventSelectionAlg::reconfigure(fhicl::ParameterSet const &p)
{
  // Implementation of optional member function here.
  m_printDebug = p.get<bool>("PrintDebug", false);

  m_nTracks = p.get<int>("nTracks", 1);
  m_fidvolXstart = p.get<double>("fidvolXstart", 10);
  m_fidvolXend = p.get<double>("fidvolXend", 10);

  m_fidvolYstart = p.get<double>("fidvolYstart", 20);
  m_fidvolYend = p.get<double>("fidvolYend", 20);

  m_fidvolZstart = p.get<double>("fidvolZstart", 10);
  m_fidvolZend = p.get<double>("fidvolZend", 50);
  m_pfp_producer = p.get<std::string>("PFParticleLabel", "pandoraNu::McRecoStage2");

  geoHelper.setFiducialVolumeCuts(m_fidvolXstart, m_fidvolXend, m_fidvolYstart,
                                  m_fidvolYend, m_fidvolZstart, m_fidvolZend);

  m_fractionsigmaflashwidth = p.get<double>("fractionsigmaflashwidth", 2.0);
  m_absoluteflashdist = p.get<double>("absoluteflashdist", 50.0);

  m_startbeamtime = p.get<double>("startbeamtime", 3.2);
  m_endbeamtime = p.get<double>("endbeamtime", 4.8);
  m_PE_threshold = p.get<double>("PE_threshold", 50.);
  m_cut_zwidth = p.get<double>("cut_zwidth", 105);
  m_cut_sigzwidth = p.get<double>("cut_sigzwidth", 1.0);
  m_cut_ywidth = p.get<double>("cut_ywidth", 95.);
  m_cut_sigywidth = p.get<double>("cut_sigywidth", 2.2);
  m_charge_light_ratio = p.get<double>("charge_light_ratio", 3.0);

  m_flashmatching = p.get<bool>("Flashmatching", true);
  m_FM_all = p.get<bool>("Flashmatching_first", true);
  _do_opdet_swap = p.get<bool>("DoOpDetSwap", false);
  _opdet_swap_map = p.get<std::vector<int>>("OpDetSwapMap");
  m_isCosmicInTime = p.get<bool>("isCosmicInTime", false);
  m_showersAsTracks = p.get<bool>("ShowersAsTracks", false);
  m_mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));
  fOpticalFlashFinderLabel = p.get<std::string>("OpticalFlashFinderLabel", "simpleFlashBeam");

}

const std::map<size_t, int> ElectronEventSelectionAlg::flashBasedSelection(const art::Event &evt,
                                                                           const std::vector<size_t> &pfplist,
                                                                           const art::ValidHandle<std::vector<recob::PFParticle>> pfparticle_handle)
{
  // All initializations
  std::vector<double> ChargeCenter;
  std::vector<flashana::QCluster_t> qcvec;
  std::vector<unsigned int> PFPIDvector; //links the pfp indices to the qvec indices.
  std::vector<double> chargexvector;
  std::vector<flashana::FlashMatch_t> matchvec;
  std::vector<double> scorevector;
  std::vector<double> TPC_x_vector;
  std::vector<unsigned int> TPCIDvector; //links the qvec indices to the matched ones, matched ones are score ordered already

  size_t chosen_index = -1;
  std::map<size_t, int> result;

  //Select the flash with the biggest PE inside the window

  art::InputTag optical_tag{fOpticalFlashFinderLabel};
  auto const &optical_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);

  int maxIndex = -1;
  double maxPE = 0;
  for (unsigned int ifl = 0; ifl < optical_handle->size(); ++ifl)
  {
    recob::OpFlash const &flash = optical_handle->at(ifl);
    _flash_PE.push_back(flash.TotalPE());
    _flash_time.push_back(flash.Time());

    
    if ((flash.Time() < m_endbeamtime && flash.Time() > m_startbeamtime))
    {
      double thisPE = flash.TotalPE();
      if (thisPE > maxPE)
      {
        maxPE = thisPE;
        maxIndex = ifl;
      }
    }
  }
  if (maxIndex == -1 || maxPE < m_PE_threshold)
  {
    std::cout << "[ElectronEventSelectionAlg] "
              << "No flash in event within window and over " << m_PE_threshold << "PE!" << std::endl;
  }

  // Else means we have a good flash
  else
  {

    // Store what I want to know about the flash
    recob::OpFlash const &flash = optical_handle->at(maxIndex);
    ::flashana::Flash_t f;
    f.x = f.x_err = 0;
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.pe_v.resize(m_geo->NOpDets());
    f.pe_err_v.resize(m_geo->NOpDets());
    f.time = flash.Time();
    for (unsigned int i = 0; i < f.pe_v.size(); i++)
    {
      unsigned int opdet = m_geo->OpDetFromOpChannel(i);
      if (_do_opdet_swap && evt.isRealData())
      {
//         std::cout << "[ElectronEventSelectionAlg] Switching the PMT mapping before flashmatching!" << std::endl;
        opdet = _opdet_swap_map.at(opdet);
      }
      f.pe_v[opdet] = flash.PE(i);
      f.pe_err_v[opdet] = sqrt(flash.PE(i));
    }

    // Loop over the neutrino candidates to do prematching cuts
    art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

    for (size_t pfpindex : pfplist)
    {
      ChargeCenter = pandoraHelper.calculateChargeCenter(pfpindex, pfparticle_handle, evt, m_pfp_producer);

      // candidates that fail the prematching cuts do not need to be passed to the manager
      bool prematching_cuts;

      prematching_cuts = (m_cut_zwidth > std::abs(ChargeCenter[2] - f.z)) && // Cut in Z direction
                         (m_cut_sigzwidth > std::abs(ChargeCenter[2] - f.z) / f.z_err) &&
                         (m_cut_ywidth > std::abs(ChargeCenter[1] - f.y)) && // Cut in Y direction
                         (m_cut_sigywidth > std::abs(ChargeCenter[1] - f.y) / f.y_err) &&
                         (m_charge_light_ratio < std::abs(ChargeCenter[3] / maxPE));

      if (prematching_cuts)
      {
        std::vector<size_t> daughters;
        pandoraHelper.traversePFParticleTree(pfparticle_handle, pfpindex, daughters, m_pfp_producer);
        qcvec.emplace_back(collect3DHits(evt, daughters));
        PFPIDvector.emplace_back(pfpindex);

        chargexvector.emplace_back(ChargeCenter[0]);
        if (m_printDebug) {
          std::cout << "[ElectronEventSelectionAlg] "
                    << "Neutrino candidate " << pfpindex << " passed (prematching cuts)." << std::endl;
        }
      }
      else
      {
        if (m_printDebug) {
          std::cout << "[ElectronEventSelectionAlg] "
                    << "Neutrino candidate " << pfpindex << " rejected (prematching cuts)." << std::endl;
        }
      }
    } // Loop over the list of neutrino candidates

    // If there are no postmatching cuts and one remaining candidate, don't run the matching!
    if (qcvec.size() == 0)
    {
      std::cout << "[ElectronEventSelectionAlg] "
                << "All neutrino candidate rejected (prematching cuts)." << std::endl;
    }

    //  Commented out for now to force flashmatching to make sure x-values are calculated.
    else if (qcvec.size() ==1 )
    {
      chosen_index = PFPIDvector[0];
    }
    else
    { // Else run flashmatching in the remaining cases
      ::flashana::Flash_t flashRecoCopy = f;
      m_mgr.Reset();
      m_mgr.Emplace(std::move(flashRecoCopy));
      for (auto cluster : qcvec)
      {
        flashana::QCluster_t clusterCopy = cluster;
        m_mgr.Emplace(std::move(clusterCopy));
      }

      matchvec = m_mgr.Match();

      if (matchvec.size() == 0)
      {
        std::cout << "[ElectronEventSelectionAlg] "
                  << "Flashmatchmanager unable to match!" << std::endl;
      }
      else
      { // Else means we have a match
        for (auto match : matchvec)
        {
          //Postmatching cuts.
          if (2. * m_geo->DetHalfWidth() < chargexvector[match.tpc_id])
          {
            std::cout << "[ElectronEventSelectionAlg] "
                      << "The x-position of the TPC object is outside the detector!" << std::endl;
          }
          else
          {
            // std::cout << "Match score " << match.score << std::endl;
            // std::cout << "flash_x " << match.tpc_point.x << std::endl;
            // std::cout << "TPC_x " << chargexvector[match.tpc_id] << std::endl;

            scorevector.emplace_back(match.score);
            TPC_x_vector.emplace_back(match.tpc_point.x);
            TPCIDvector.emplace_back(match.tpc_id);
          }
        }

        chosen_index = PFPIDvector[TPCIDvector[0]];

        _TPC_x = chargexvector[TPCIDvector[0]];
        _flash_x = TPC_x_vector[0];
        std::cout << "[ElectronEventSelectionAlg] "
                  << "Candidate " << PFPIDvector[TPCIDvector[0]]
                  << " passed optical selection! TPC_X: " << _TPC_x << " flash_x " << _flash_x << std::endl;
      } // Else means we have a match

    } // Else run flashmatching in the remaining cases

  } // Else means we have a good flash

  for (unsigned int i = 0; i < pfplist.size(); i++)
  {
    if (pfplist[i] == chosen_index)
    {
      result[pfplist[i]] = maxIndex; //MaxIndex is the index corresponding to the flash
    }
    else
    {
      result[pfplist[i]] = -1;
    }
  }
  return result;
} // End of flashbased selection function

const flashana::QCluster_t ElectronEventSelectionAlg::collect3DHits(
    const art::Event &evt,
    const std::vector<size_t> &pfplist)
{

  flashana::QCluster_t cluster;

  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);

  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, m_pfp_producer);

  for (auto &pfpindex : pfplist)
  {
    //unsigned short pdgcode = pfparticle_handle->at(pfpindex).PdgCode();
    // double lycoef = m_ly_map[pdgcode];
    double lycoef = 1.0;

    std::vector<flashana::Hit3D_t> hitlist;
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pfpindex);

    // Loop over the spacepoints and get the associated hits:
    for (auto &_sps : spcpnts)
    {
      std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
      // Add the hits to the weighted average, if they are collection hits:
      for (auto &hit : hits)
      {
        if (hit->View() == geo::kZ)
        {
          // Collection hits only
          auto xyz = _sps->XYZ();
          flashana::Hit3D_t hit3D; //Collection plane hits
          hit3D.x = xyz[0];
          hit3D.y = xyz[1];
          hit3D.z = xyz[2];
          hit3D.plane = 2;
          double q = hit->Integral();

          hit3D.q = q;
          hitlist.emplace_back(hit3D);
        }
      }
    }
    cluster += ((flashana::LightCharge *)(m_mgr.GetCustomAlgo("LightCharge")))->FlashHypothesisCharge(hitlist, lycoef);
  }
  return cluster;
}

const std::map<size_t, int> ElectronEventSelectionAlg::opticalfilter(const art::Event &evt,
                                                                     const std::vector<size_t> &pfplist,
                                                                     const art::ValidHandle<std::vector<recob::PFParticle>> pfparticle_handle)
{
  // All initializations
  std::vector<double> ChargeCenter;
  std::map<size_t, int> result;

  art::InputTag optical_tag{fOpticalFlashFinderLabel};
  auto const &optical_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);

  // Loop over pfp neutrino candidates and flashes.
  for (size_t pfp_i : pfplist)
  {
    result[pfp_i] = -1;
    for (unsigned int ifl = 0; ifl < optical_handle->size(); ++ifl)
    {
      recob::OpFlash const &flash = optical_handle->at(ifl);
      _flash_PE.push_back(flash.TotalPE());
      _flash_time.push_back(flash.Time());

      // Request flash in time window
      if ((flash.Time() < m_endbeamtime && flash.Time() > m_startbeamtime))
      {

        ChargeCenter = pandoraHelper.calculateChargeCenter(pfp_i, pfparticle_handle, evt, m_pfp_producer);

        // Cut on the z position
        double absolute = std::abs(flash.ZCenter() - ChargeCenter[2]);
        double sigma = absolute / flash.ZWidth();
        // std::cout << "z_diff: " << absolute << " sigma: " << sigma << std::endl;

        if (absolute < m_absoluteflashdist || sigma < 1. / m_fractionsigmaflashwidth)
        {
          result[pfp_i] = ifl;
          // std::cout << "candidate " << pfp_i << " passed opt cut with flash " << ifl << std::endl;
        }
      }
    }
  }
  return result;
}

bool ElectronEventSelectionAlg::eventSelected(const art::Event &evt)
{

  clear();

  std::cout << "[ElectronEventSelectionAlg] Selection start" << std::endl;

  if (m_printDebug) {
      std::cout << "[ElectronEventSelectionAlg] START BEAM TIME " << m_startbeamtime << std::endl; 
  }

  // Get the list of pfparticles:
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);

  // Are there any pfparticles?
  if (pfparticle_handle->size() == 0)
  {
    std::cout << "[ElectronEventSelectionAlg] "
              << "NO RECO DATA PRODUCTS" << std::endl;
    return false;
  }

  // Get the list of primary pfparticles that are also neutrinos
  for (size_t _i_pfp = 0; _i_pfp < pfparticle_handle->size(); _i_pfp++)
  {
    if ((abs(pfparticle_handle->at(_i_pfp).PdgCode()) == 12 ||
         abs(pfparticle_handle->at(_i_pfp).PdgCode()) == 14 ||
         abs(pfparticle_handle->at(_i_pfp).PdgCode()) == 16) &&
        pfparticle_handle->at(_i_pfp).IsPrimary())
    {
      _primary_indexes.push_back(_i_pfp);
    }
  }

  // If there are no particles flagged as primary, return false
  if (_primary_indexes.size() == 0)
  {
    return false;
  }

  _n_neutrino_candidates = _primary_indexes.size();
  std::cout << "[ElectronEventSelectionAlg] "
            << "Primary PFParticles " << _n_neutrino_candidates << std::endl;
  // For each of the primary particles, determine if it and it's daughters pass
  // the cuts:

  // Need associations from pfparticle to vertex
  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt,
                                                 m_pfp_producer);

  for (auto &_i_primary : _primary_indexes)
  {
    unsigned int numDaughters = pfparticle_handle->at(_i_primary).NumDaughters();
    
    if (m_printDebug) {
      std::cout << "[ElectronEventSelectionAlg] "
                << "Primary PDG " << pfparticle_handle->at(_i_primary).PdgCode()
                << std::endl;
      std::cout << "[ElectronEventSelectionAlg] "
                << "N. of Daughters "
                << numDaughters << std::endl;
    }

    _neutrino_candidate_passed[_i_primary] = true;
    _neutrino_vertex[_i_primary] = TVector3(0, 0, 0);
    _n_showers[_i_primary] = 0;
    _pfp_id_showers_from_primary[_i_primary] = std::vector<size_t>();
    _n_tracks[_i_primary] = 0;
    _pfp_id_tracks_from_primary[_i_primary] = std::vector<size_t>();
    _n_showers_as_tracks[_i_primary] = 0;
    _pfp_id_showers_as_tracks_from_primary[_i_primary] = std::vector<size_t>();

    // Get the neutrino vertex and check if it's fiducial:
    std::vector<double> neutrino_vertex;
    neutrino_vertex.resize(3);
    try
    {
      auto const &neutrino_vertex_obj = vertex_per_pfpart.at(_i_primary);
      neutrino_vertex_obj->XYZ(
          &neutrino_vertex[0]); // PFParticle neutrino vertex coordinates

      // Save it as a TVector3:
      _neutrino_vertex.at(_i_primary).SetX(neutrino_vertex[0]);
      _neutrino_vertex.at(_i_primary).SetY(neutrino_vertex[1]);
      _neutrino_vertex.at(_i_primary).SetZ(neutrino_vertex[2]);

      if (!geoHelper.isFiducial(_neutrino_vertex.at(_i_primary)))
      {
        _neutrino_candidate_passed[_i_primary] = false;
        std::cout << "[ElectronEventSelectionAlg] "
                  << "Neutrino vertex not within fiducial volume" << std::endl;
      }
    }
    catch (...)
    {
      std::cout << "[ElectronEventSelectionAlg] "
                << "NO VERTEX AVAILABLE " << std::endl;
      _neutrino_candidate_passed[_i_primary] = false;
    }

    int showers = 0;                // number of showers in the hierarchy of the pfp neutrino candidate.
    int tracks = 0;
    int shower_daughters = 0;       // number of showers that are direct daughters of the pfp neutrino candidate.
    int track_daughters = 0;
    int showers_as_tracks = 0;

    std::vector<size_t> daughters_id;
    pandoraHelper.traversePFParticleTree(pfparticle_handle, _i_primary, daughters_id, m_pfp_producer);

    for (auto const &pfdaughter : daughters_id)
    {
      
      if (m_printDebug) {
        std::cout << "[ElectronEventSelectionAlg] "
                  << "Daughter ID: " << pfdaughter << " PDG "
                  << pfparticle_handle->at(pfdaughter).PdgCode() << std::endl;
      }

      if (pfparticle_handle->at(pfdaughter).PdgCode() == 11)
      {
          art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt,
                                                         m_pfp_producer);
          auto const &shower_obj = shower_per_pfpart.at(pfdaughter);

          if (shower_obj.isNull())
            continue;

          if(pfparticle_handle->at(pfdaughter).Parent()==_i_primary)
          {
            shower_daughters++;
          }
          _pfp_id_showers_from_primary[_i_primary].push_back(pfdaughter);
          showers++;
      }

      if (pfparticle_handle->at(pfdaughter).PdgCode() == 11 && m_showersAsTracks)
      {
        art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

        auto const &shower_as_track_obj = track_per_pfpart.at(pfdaughter);

        if (shower_as_track_obj.isNull())
          continue;

        _pfp_id_showers_as_tracks_from_primary[_i_primary].push_back(pfdaughter);
        showers_as_tracks++;
      }

      if (pfparticle_handle->at(pfdaughter).PdgCode() == 13)
      {
          art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
          auto const &track_obj = track_per_pfpart.at(pfdaughter);
          if (track_obj.isNull())
            continue;

          if(pfparticle_handle->at(pfdaughter).Parent()==_i_primary)
          {
            track_daughters++;
          }
          _pfp_id_tracks_from_primary[_i_primary].push_back(pfdaughter);
          tracks++;
      }

    }
    _n_showers_as_tracks[_i_primary] = showers_as_tracks;

    _n_tracks[_i_primary] = tracks;
    _n_showers[_i_primary] = showers;

    if (m_printDebug) {
      std::cout << "[ElectronEventSelectionAlg] "
                << "Showers: " << showers << ", tracks: " << tracks << std::endl;
    }

    // Cut on the topology to select 1e Np like signal
    // Np: at least N direct track daughters
    // 1e: at least one direct shower or 1 direct track daughter with a shower connected to it

    if ((track_daughters + tracks + showers + shower_daughters) < 1) {
      _neutrino_candidate_passed[_i_primary] = false;
    }


    if (track_daughters < m_nTracks){
      // There are less direct daughter tracks than we want protons, FAIL   
      std::cout << "[ElectronEventSelectionAlg] There are less direct daughter tracks than we want protons, FAIL" << std::endl;             
      _neutrino_candidate_passed[_i_primary] = false;
    }
    if (showers==0){
      // There are no showers in the complete hierarchy, FAIL (this is probably where we lose LEE efficiency)
      std::cout << "[ElectronEventSelectionAlg] There are no showers in the complete hierarchy, FAIL" << std::endl;             

      _neutrino_candidate_passed[_i_primary] = false;
    }
    if (shower_daughters==0 && showers==1 && track_daughters < (m_nTracks+1) ){
        std::cout << "[ElectronEventSelectionAlg] There are no direct showers, but there is one shower in the hierarchy, this means we require N+1 direct tracks, otherwise, FAIL" << std::endl;

      // There are no direct showers, but there is one shower in the hierarchy, this means we require N+1 direct tracks, otherwise, FAIL
      _neutrino_candidate_passed[_i_primary] = false;
    }
    // Else, the pfp neutrino candidate passes the required topology!
  }

  std::map<size_t, int> optical_map;
  std::vector<size_t> pfplist;

  if (m_FM_all)
  {
    // Do flashmatching for all candidates
    pfplist = _primary_indexes;
  }
  else
  {
    // Fill PFPlist with particles that passed the track shower and fiducial cut requirements.

    for (auto val : _neutrino_candidate_passed)
    {
      if (val.second)
      {
        pfplist.push_back(val.first);
      }
      else
      {
        optical_map[val.first] = -1;
      }
    }
  }

  if (pfplist.size() > 0)
  {
    if (m_flashmatching)
    {
      _op_flash_indexes = flashBasedSelection(evt, pfplist, pfparticle_handle);
    }
    else
    {
      _op_flash_indexes = opticalfilter(evt, pfplist, pfparticle_handle);
    }
  }

  _op_flash_indexes.insert(optical_map.begin(), optical_map.end());

  for (auto &val : _neutrino_candidate_passed)
  {
    if (val.second)
    {
      if (_op_flash_indexes[val.first] < 0)
      {
        val.second = false;
      }
    }
  }

  // Last, determine if any primary particles passed:
  for (auto val : _neutrino_candidate_passed)
  {
    if (m_printDebug) {
        std::cout << "[ElectronEventSelectionAlg] candidate " << val.second << std::endl;
    }
    if (val.second)
    {
      std::cout << "[ElectronEventSelectionAlg] "
                << "EVENT SELECTED" << std::endl;
      return true;
    }
  }

  return false;
}

} // lee

#endif // ELECTRON_EVENT_SELECTION_ALG_CXX
