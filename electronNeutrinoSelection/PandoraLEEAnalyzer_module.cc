////////////////////////////////////////////////////////////////////////
// Class:       PandoraLEEAnalyzer
// Module Type: analyzer
// File:        PandoraLEEAnalyzer_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "PandoraLEEAnalyzer.h"

lee::PandoraLEEAnalyzer::PandoraLEEAnalyzer(fhicl::ParameterSet const &pset)
    : EDAnalyzer(pset) // ,
// More initializers here.
{
  for (int i = 0; i < _h_lee_scaling->GetNbinsX(); i++) {
    _h_lee_scaling->SetBinContent(i+1, _lee_scaling[i]);
  }
  // create output tree
  art::ServiceHandle<art::TFileService> tfs;
  // myTFile = new TFile("PandoraLEEAnalyzerOutput.root", "RECREATE");
  myTTree = tfs->make<TTree>("pandoratree", "PandoraAnalysis Tree");

  myPOTTTree = tfs->make<TTree>("pot", "POT Tree");

  myTTree->Branch("weights", "std::map< std::string, std::vector< double > >", &_weights);
  myTTree->Branch("flux_weights", "std::map< std::string, std::vector< double > >", &_flux_weights);

  myTTree->Branch("category", &_category, "category/I");

  myTTree->Branch("n_tracks", &_n_tracks, "n_tracks/I");
  myTTree->Branch("n_showers", &_n_showers, "n_showers/I");
  myTTree->Branch("ccnc", &_ccnc, "ccnc/I");
  myTTree->Branch("cosmic_fraction", &_cosmic_fraction, "cosmic_fraction/d");

  myTTree->Branch("vx", &_vx, "vx/d");
  myTTree->Branch("vy", &_vy, "vy/d");
  myTTree->Branch("vz", &_vz, "vz/d");

  myTTree->Branch("true_vx", &_true_vx, "true_vx/d");
  myTTree->Branch("true_vy", &_true_vy, "true_vy/d");
  myTTree->Branch("true_vz", &_true_vz, "true_vz/d");

  myTTree->Branch("true_shower_x_sce", "std::vector< double >", &_true_shower_x_sce);
  myTTree->Branch("true_shower_y_sce", "std::vector< double >", &_true_shower_y_sce);
  myTTree->Branch("true_shower_z_sce", "std::vector< double >", &_true_shower_z_sce);
  myTTree->Branch("true_shower_pdg", "std::vector< int >", &_true_shower_pdg);
  myTTree->Branch("true_shower_depE", "std::vector< double >", &_true_shower_depE);

  myTTree->Branch("true_vx_sce", &_true_vx_sce, "true_vx_sce/d");
  myTTree->Branch("true_vy_sce", &_true_vy_sce, "true_vy_sce/d");
  myTTree->Branch("true_vz_sce", &_true_vz_sce, "true_vz_sce/d");

  myTTree->Branch("nu_E", &_nu_energy, "nu_E/d");
  myTTree->Branch("passed", &_event_passed, "passed/I");
  myTTree->Branch("numu_passed", &_numu_passed, "numu_passed/I");
  myTTree->Branch("numu_cuts", &_numu_cuts, "numu_cuts/I");

  myTTree->Branch("n_total_candidates", &_n_total_candidates, "n_total_candidates/i");
  myTTree->Branch("candidate_vx", "std::vector< double >", &_candidate_vx);
  myTTree->Branch("candidate_vy", "std::vector< double >", &_candidate_vy);
  myTTree->Branch("candidate_vz", "std::vector< double >", &_candidate_vz);

  myTTree->Branch("n_candidates", &_n_candidates, "n_candidates/I");
  myTTree->Branch("n_true_nu", &_n_true_nu, "n_true_nu/I");
  myTTree->Branch("distance", &_distance, "distance/d");
  myTTree->Branch("true_nu_is_fiducial", &_true_nu_is_fiducial,
                  "true_nu_is_fiducial/I");

  myTTree->Branch("n_matched", &_n_matched, "n_matched/I");
  myTTree->Branch("nu_matched_tracks", &_nu_matched_tracks,
                  "nu_matched_tracks/I");
  myTTree->Branch("nu_matched_showers", &_nu_matched_showers,
                  "nu_matched_showers/I");

  myTTree->Branch("nu_daughters_pdg", "std::vector< int >", &_nu_daughters_pdg);
  myTTree->Branch("nu_daughters_E", "std::vector< double >", &_nu_daughters_E);

  myTTree->Branch("nu_daughters_vx", "std::vector< double >",
                  &_nu_daughters_vx);
  myTTree->Branch("nu_daughters_vy", "std::vector< double >",
                  &_nu_daughters_vy);
  myTTree->Branch("nu_daughters_vz", "std::vector< double >",
                  &_nu_daughters_vz);

  myTTree->Branch("nu_daughters_endx", "std::vector< double >",
                  &_nu_daughters_endx);
  myTTree->Branch("nu_daughters_endy", "std::vector< double >",
                  &_nu_daughters_endy);
  myTTree->Branch("nu_daughters_endz", "std::vector< double >",
                  &_nu_daughters_endz);

  myTTree->Branch("nu_daughters_px", "std::vector< double >",
                  &_nu_daughters_px);
  myTTree->Branch("nu_daughters_py", "std::vector< double >",
                  &_nu_daughters_py);
  myTTree->Branch("nu_daughters_pz", "std::vector< double >",
                  &_nu_daughters_pz);

  myTTree->Branch("nu_track_ids", "std::vector< size_t >",
                  &_nu_track_ids);
  myTTree->Branch("nu_shower_ids", "std::vector< size_t >",
                  &_nu_shower_ids);

  myTTree->Branch("nu_shower_daughters", "std::vector< std::vector< size_t > >",
                  &_nu_shower_daughters);
  myTTree->Branch("nu_track_daughters", "std::vector< std::vector< size_t > >",
                  &_nu_track_daughters);

  myTTree->Branch("event", &_event, "event/I");
  myTTree->Branch("run", &_run, "run/I");
  myTTree->Branch("subrun", &_subrun, "subrun/I");
  myTTree->Branch("leeweight", &_leeweight, "leeweight/f");

  myTTree->Branch("bnbweight", &_bnbweight, "bnbweight/d");

  myTTree->Branch("chosen_candidate", &_chosen_candidate, "chosen_candidate/I");
  myTTree->Branch("n_primaries", &_n_primaries, "n_primaries/I");

  myTTree->Branch("primary_indexes", "std::vector< int >", &_primary_indexes);
  myTTree->Branch("number_tracks", "std::vector< int >", &_number_tracks);
  myTTree->Branch("number_showers", "std::vector< int >", &_number_showers);

  myTTree->Branch("flash_time", "std::vector< double >", &_flash_time);
  myTTree->Branch("flash_PE", "std::vector< double >", &_flash_PE);
  myTTree->Branch("flash_x", &_flash_x, "flash_x/d");
  myTTree->Branch("TPC_x", &_TPC_x, "TPC_x/d");

  myTTree->Branch("flash_passed", "std::vector< int >", &_flash_passed);
  myTTree->Branch("track_passed", "std::vector< int >", &_track_passed);
  myTTree->Branch("shower_passed", "std::vector< int >", &_shower_passed);

  myTTree->Branch("shower_distance", "std::vector< double >", &_shower_distance);
  myTTree->Branch("shower_dir_x", "std::vector< double >", &_shower_dir_x);
  myTTree->Branch("shower_dir_y", "std::vector< double >", &_shower_dir_y);
  myTTree->Branch("shower_dir_z", "std::vector< double >", &_shower_dir_z);

  myTTree->Branch("shower_start_x", "std::vector< double >", &_shower_start_x);
  myTTree->Branch("shower_start_y", "std::vector< double >", &_shower_start_y);
  myTTree->Branch("shower_start_z", "std::vector< double >", &_shower_start_z);

  myTTree->Branch("shower_theta", "std::vector< double >", &_shower_theta);
  myTTree->Branch("shower_phi", "std::vector< double >", &_shower_phi);
  myTTree->Branch("shower_n_clusters", "std::vector< double >", &_shower_n_clusters);

  myTTree->Branch("shower_energy", "std::vector< std::vector< double > >", &_shower_energy);
  // myTTree->Branch("track_energy_dedx", "std::vector< double >", &_track_energy_dedx);
  myTTree->Branch("track_energy_hits", "std::vector< std::vector< double > >", &_track_energy_hits);

  myTTree->Branch("track_distance", "std::vector< double >", &_track_distance);
  myTTree->Branch("track_dir_x", "std::vector< double >", &_track_dir_x);
  myTTree->Branch("track_dir_y", "std::vector< double >", &_track_dir_y);
  myTTree->Branch("track_dir_z", "std::vector< double >", &_track_dir_z);

  myTTree->Branch("track_start_x", "std::vector< double >", &_track_start_x);
  myTTree->Branch("track_start_y", "std::vector< double >", &_track_start_y);
  myTTree->Branch("track_start_z", "std::vector< double >", &_track_start_z);

  myTTree->Branch("track_is_fiducial", "std::vector< int >",
                  &_track_is_fiducial);
  myTTree->Branch("shower_is_fiducial", "std::vector< int >",
                  &_shower_is_fiducial);

  myTTree->Branch("track_end_x", "std::vector< double >", &_track_end_x);
  myTTree->Branch("track_end_y", "std::vector< double >", &_track_end_y);
  myTTree->Branch("track_end_z", "std::vector< double >", &_track_end_z);

  myTTree->Branch("track_theta", "std::vector< double >", &_track_theta);
  myTTree->Branch("track_phi", "std::vector< double >", &_track_phi);

  myTTree->Branch("track_len", "std::vector< double >", &_track_length);
  myTTree->Branch("track_id", "std::vector< double >", &_track_id);

  myTTree->Branch("track_bragg_p", "std::vector< double >", &_track_bragg_p);
  myTTree->Branch("track_bragg_mu", "std::vector< double >", &_track_bragg_mu);
  myTTree->Branch("track_bragg_mip", "std::vector< double >", &_track_bragg_mip);
  myTTree->Branch("track_pida", "std::vector< double >", &_track_pida);
  myTTree->Branch("track_pid_chipr", "std::vector< double >", &_track_pidchipr);
  myTTree->Branch("track_pid_chipi", "std::vector< double >", &_track_pidchipi);
  myTTree->Branch("track_pid_chika", "std::vector< double >", &_track_pidchika);
  myTTree->Branch("track_pid_chimu", "std::vector< double >", &_track_pidchimu);

  myTTree->Branch("track_res_mean", "std::vector< double >", &_track_res_mean);
  myTTree->Branch("track_res_std", "std::vector< double >", &_track_res_std);
  myTTree->Branch("track_n_clusters", "std::vector< double >", &_track_n_clusters);

  myTTree->Branch("shower_res_mean", "std::vector< double >", &_shower_res_mean);
  myTTree->Branch("shower_res_std", "std::vector< double >", &_shower_res_std);
  myTTree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");

  myPOTTTree->Branch("run", &_run_sr, "run/I");
  myPOTTTree->Branch("subrun", &_subrun_sr, "subrun/I");
  myPOTTTree->Branch("pot", &_pot, "pot/d");

  myTTree->Branch("interaction_type", &_interaction_type, "interaction_type/I");

  myTTree->Branch("shower_dQdx", "std::vector< std::vector< double > >",
                  &_shower_dQdx);
  myTTree->Branch("shower_dEdx", "std::vector< std::vector< double > >",
                  &_shower_dEdx);

  myTTree->Branch("shower_dQdx_cali", "std::vector< std::vector< double > >", &_shower_dQdx_cali);

  myTTree->Branch("track_dQdx", "std::vector< std::vector< double > >",
                  &_track_dQdx);

  myTTree->Branch("track_dEdx", "std::vector< std::vector< double > >",
                  &_track_dEdx);

  myTTree->Branch("shower_open_angle", "std::vector< double >",
                  &_shower_open_angle);

  myTTree->Branch("shower_length", "std::vector< double >",
                  &_shower_length);

  myTTree->Branch("shower_dQdx_hits", "std::vector< std::vector< double > >",
                  &_shower_dQdx_hits);

  myTTree->Branch("shower_pitches", "std::vector< std::vector< double > >",
                  &_shower_pitches);

  myTTree->Branch("shower_dQdx_hits_in_the_box", "std::vector< std::vector< int > >",
                  &_shower_dQdx_hits_in_the_box);

  myTTree->Branch("shower_dEdx_hits", "std::vector< std::vector< double > >",
                  &_shower_dEdx_hits);

  myTTree->Branch("track_dQdx_hits", "std::vector< std::vector< double > >",
                  &_track_dQdx_hits);
  myTTree->Branch("track_dEdx_hits", "std::vector< std::vector< double > >",
                  &_track_dEdx_hits);

  myTTree->Branch("matched_tracks", "std::vector< int >",
                  &_matched_tracks);
  myTTree->Branch("matched_tracks_energy", "std::vector< double >",
                  &_matched_tracks_energy);
  myTTree->Branch("matched_tracks_process", "std::vector< std::string >",
                  &_matched_tracks_process);

  myTTree->Branch("matched_showers", "std::vector< int >",
                  &_matched_showers);
  myTTree->Branch("matched_showers_process", "std::vector< std::string >",
                  &_matched_showers_process);
  myTTree->Branch("matched_showers_energy", "std::vector< double >",
                  &_matched_showers_energy);

  myTTree->Branch("shower_pca", "std::vector< std::vector< double > >",
                  &_shower_pca);

  myTTree->Branch("shower_nhits", "std::vector< std::vector<int> >",
                  &_shower_nhits);

  myTTree->Branch("track_pca", "std::vector< std::vector< double > >",
                  &_track_pca);

  myTTree->Branch("track_nhits", "std::vector< std::vector<int> >",
                  &_track_nhits);

  myTTree->Branch("shower_energy_cali", "std::vector< std::vector< double > >", &_shower_energy_cali);
  myTTree->Branch("track_energy_cali", "std::vector< std::vector< double > >", &_track_energy_cali);

  this->reconfigure(pset);
}

lee::PandoraLEEAnalyzer::~PandoraLEEAnalyzer()
{
  std::cout << "[PandoraLEE] "
            << "End!" << std::endl;
}

art::Ptr<recob::Shower> lee::PandoraLEEAnalyzer::get_most_energetic_shower(
    std::vector<art::Ptr<recob::Shower>> &showers)
{
  art::Ptr<recob::Shower> most_energetic_shower;

  double max_energy = std::numeric_limits<double>::lowest();
  for (auto const &shower : showers)
  {
    if (shower->Energy()[shower->best_plane()] > max_energy)
    {
      most_energetic_shower = shower;
      max_energy = shower->Energy()[shower->best_plane()];
    }
  }
  return most_energetic_shower;
}

art::Ptr<recob::Track> lee::PandoraLEEAnalyzer::get_longest_track(
    std::vector<art::Ptr<recob::Track>> &tracks)
{
  art::Ptr<recob::Track> longest_track;

  double max_length = std::numeric_limits<double>::lowest();
  for (auto const &track : tracks)
  {
    try
    {
      if (track->Length() > max_length)
      {
        longest_track = track;
        max_length = track->Length();
      }
    }
    catch (...)
    {
      std::cout << "[PandoraLEE] "
                << "Error getting longest track " << track << std::endl;
    }
  }
  return longest_track;
}

size_t
lee::PandoraLEEAnalyzer::choose_candidate(std::vector<size_t> &candidates,
                                          const art::Event &evt)
{
  if(candidates.size()==1){
    return candidates[0];
  }

  size_t chosen_candidate = 0;
  double most_z = -1;
  double longest_track_dir;

  for (auto const &ic : candidates)
  {
    std::vector<art::Ptr<recob::Track>> nu_tracks;
    std::cout << "[PandoraLEE] Candidate " << ic << std::endl;
    std::vector<size_t> _nu_track_ids = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(ic);
    pandoraHelper.get_daughter_tracks(_nu_track_ids, evt, nu_tracks, m_pfp_producer);
    longest_track_dir = get_longest_track(nu_tracks)->StartDirection().Z();

    if (longest_track_dir > most_z)
    {
      chosen_candidate = ic;
      most_z = longest_track_dir;
    }
  }

  return chosen_candidate;
}

void lee::PandoraLEEAnalyzer::endSubRun(const art::SubRun &sr)
{

  _run_sr = sr.run();
  _subrun_sr = sr.subRun();

  art::Handle<sumdata::POTSummary> potListHandle;
  if (!m_isData)
  {
    if (sr.getByLabel(_mctruthLabel, potListHandle))
      _pot = potListHandle->totpot;
    else
      _pot = 0.;
  }
  else
  {
    if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle))
      _pot = potListHandle->totpot;
    else
      _pot = 0.;
  }

  myPOTTTree->Fill();
}
void lee::PandoraLEEAnalyzer::clear()
{
  _track_res_mean.clear();
  _track_res_std.clear();
  _shower_res_mean.clear();
  _shower_res_std.clear();
  _weights.clear();
  _flux_weights.clear();

  _shower_pitches.clear();
  _shower_dQdx_hits_in_the_box.clear();
  _ccnc = std::numeric_limits<int>::lowest();
  _interaction_type = std::numeric_limits<int>::lowest();
  _track_bragg_p.clear();
  _track_bragg_mu.clear();
  _track_bragg_mip.clear();
  _track_pida.clear();
  _track_pidchipr.clear();
  _track_pidchika.clear();
  _track_pidchipi.clear();
  _track_pidchimu.clear();
  _track_pidchi.clear();
  _shower_pca.clear();
  _track_pca.clear();
  _shower_nhits.clear();
  _track_nhits.clear();
  _matched_tracks.clear();
  _matched_tracks_process.clear();
  _matched_tracks_energy.clear();

  _matched_showers.clear();
  _matched_showers_process.clear();
  _matched_showers_energy.clear();

  _shower_dQdx_hits.clear();
  _shower_dEdx_hits.clear();
  _shower_dQdx.clear();
  _shower_dEdx.clear();

  _shower_dQdx_cali.clear();
  _shower_dEdx_cali.clear();

  _track_dQdx_hits.clear();
  _track_dEdx_hits.clear();
  _track_dQdx.clear();

  _track_dEdx.clear();


  _shower_open_angle.clear();
  _shower_length.clear();
  _shower_distance.clear();
  _shower_dir_x.clear();
  _shower_dir_y.clear();
  _shower_dir_z.clear();
  _shower_n_clusters.clear();

  _n_total_candidates = std::numeric_limits<int>::lowest();
  _candidate_vx.clear();
  _candidate_vy.clear();
  _candidate_vz.clear();

  _nu_track_ids.clear();
  _nu_shower_ids.clear();
  _nu_track_daughters.clear();
  _nu_shower_daughters.clear();

  _shower_start_x.clear();
  _shower_start_y.clear();
  _shower_start_z.clear();

  _track_is_fiducial.clear();
  _shower_is_fiducial.clear();

  _track_distance.clear();
  _track_start_x.clear();
  _track_start_y.clear();
  _track_start_z.clear();

  _track_n_clusters.clear();

  _true_shower_x_sce.clear();
  _true_shower_y_sce.clear();
  _true_shower_z_sce.clear();
  _true_shower_pdg.clear();
  _true_shower_depE.clear();

  _track_end_x.clear();
  _track_end_y.clear();
  _track_end_z.clear();

  _shower_theta.clear();
  _shower_phi.clear();

  _track_dir_x.clear();
  _track_dir_y.clear();
  _track_dir_z.clear();

  _track_theta.clear();
  _track_phi.clear();

  _shower_energy.clear();
  _track_energy_dedx.clear();
  _track_energy_hits.clear();

  _track_length.clear();
  _track_id.clear();

  _n_true_pions = 0;
  _n_true_protons = 0;
  _nu_pdg = 0;

  _flash_passed.clear();
  _track_passed.clear();
  _shower_passed.clear();
  _primary_indexes.clear();
  _number_tracks.clear();
  _number_showers.clear();

  _flash_PE.clear();
  _flash_time.clear();

  _chosen_candidate = std::numeric_limits<int>::lowest();
  _n_primaries = 0;

  _shower_energy_cali.clear();
  _track_energy_cali.clear();

  _true_nu_is_fiducial = std::numeric_limits<int>::lowest();
  _nu_energy = std::numeric_limits<double>::lowest();
  _nu_theta = std::numeric_limits<double>::lowest();
  _nu_phi = std::numeric_limits<double>::lowest();

  _n_tracks = std::numeric_limits<int>::lowest();
  _n_showers = std::numeric_limits<int>::lowest();

  _vx = std::numeric_limits<double>::lowest();
  _vy = std::numeric_limits<double>::lowest();
  _vz = std::numeric_limits<double>::lowest();

  _true_vx = std::numeric_limits<double>::lowest();
  _true_vy = std::numeric_limits<double>::lowest();
  _true_vz = std::numeric_limits<double>::lowest();

  _true_vx_sce = std::numeric_limits<double>::lowest();
  _true_vy_sce = std::numeric_limits<double>::lowest();
  _true_vz_sce = std::numeric_limits<double>::lowest();

  _nu_matched_tracks = std::numeric_limits<int>::lowest();
  _nu_matched_showers = std::numeric_limits<int>::lowest();

  _category = std::numeric_limits<int>::lowest();
  _run = std::numeric_limits<int>::lowest();
  _subrun = std::numeric_limits<int>::lowest();
  _event = std::numeric_limits<int>::lowest();
  _n_candidates = std::numeric_limits<int>::lowest();
  _n_true_nu = std::numeric_limits<int>::lowest();
  _run_sr = std::numeric_limits<int>::lowest();
  _subrun_sr = std::numeric_limits<int>::lowest();
  _n_matched = std::numeric_limits<int>::lowest();
  _pot = std::numeric_limits<double>::lowest();
  _event_passed = 0;
  _numu_passed = 0;
  _numu_cuts = 0;
  _cosmic_fraction = std::numeric_limits<double>::lowest();
  _distance = std::numeric_limits<double>::lowest();

  _flash_x = std::numeric_limits<double>::lowest();
  _TPC_x = std::numeric_limits<double>::lowest();

  _nu_daughters_E.clear();
  _nu_daughters_pdg.clear();

  _nu_daughters_px.clear();
  _nu_daughters_py.clear();
  _nu_daughters_pz.clear();

  _nu_daughters_vx.clear();
  _nu_daughters_vy.clear();
  _nu_daughters_vz.clear();

  _nu_daughters_endx.clear();
  _nu_daughters_endy.clear();
  _nu_daughters_endz.clear();
  _leeweight = std::numeric_limits<int>::lowest();
  _bnbweight = std::numeric_limits<int>::lowest();
}

void lee::PandoraLEEAnalyzer::categorizePFParticles(
    art::Event const &evt,
    std::vector<int> &neutrino_pdg,
    std::vector<std::string> &neutrino_process,
    std::vector<double> &neutrino_energy,
    std::vector<art::Ptr<recob::PFParticle>> &neutrino_pf,

    std::vector<int> &cosmic_pdg,
    std::vector<std::string> &cosmic_process,
    std::vector<double> &cosmic_energy,
    std::vector<art::Ptr<recob::PFParticle>> &cosmic_pf)
{

  lar_pandora::PFParticlesToMCParticles matchedParticles;
  pandoraHelper.Configure(evt, m_pfp_producer, m_spacepointLabel,
                          m_hitfinderLabel, _geantModuleLabel, m_hitmatching_producer);
  pandoraHelper.GetRecoToTrueMatches(matchedParticles);

  // art::ServiceHandle<cheat::BackTracker> bt;

  for (lar_pandora::PFParticlesToMCParticles::const_iterator
           iter = matchedParticles.begin(),
           iterEnd = matchedParticles.end();
       iter != iterEnd; ++iter)
  {

    art::Ptr<simb::MCParticle> mc_par = iter->second; // The MCParticle
    art::Ptr<recob::PFParticle> pf_par = iter->first; // The matched PFParticle

    const auto mc_truth =
        pandoraHelper.TrackIDToMCTruth(evt, _geantModuleLabel, mc_par->TrackId());

    if (mc_truth->Origin() == simb::kBeamNeutrino)
    {
      neutrino_pf.push_back(pf_par);
      neutrino_pdg.push_back(mc_par->PdgCode());
      if (mc_par->EndProcess() != NULL)
      {
        neutrino_process.push_back(mc_par->EndProcess());
      }
      neutrino_energy.push_back(mc_par->E());
    }

    if (mc_truth->Origin() == simb::kCosmicRay)
    {
      cosmic_pf.push_back(pf_par);
      cosmic_pdg.push_back(mc_par->PdgCode());
      if (mc_par->EndProcess() != NULL)
      {
        cosmic_process.push_back(mc_par->EndProcess());
      }
      cosmic_energy.push_back(mc_par->E());
    }
  }
}


void lee::PandoraLEEAnalyzer::analyze(art::Event const &evt)
{
  clear();
  auto const *sce_service = lar::providerFrom<spacecharge::SpaceChargeService>();

  _run = evt.run();
  _subrun = evt.subRun();
  _event = evt.id().event();

  std::cout << "[PandoraLEEAnalyzer] "
            << "RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event
            << std::endl;

  std::vector<size_t> nu_candidates;

  _event_passed = (int)fElectronEventSelectionAlg.eventSelected(evt);

  if (_event_passed && !evt.isRealData() && m_save_flux_info) {
    std::cout<< "[PandoraLEEAnalyzer] Saving Flux info..." << std::endl;
    art::Handle<std::vector<simb::GTruth>> gTruthHandle;
    evt.getByLabel(_mctruthLabel, gTruthHandle);
    if (!gTruthHandle.isValid())
      return;
    std::vector<art::Ptr<simb::GTruth>> gTruthVec;
    art::fill_ptr_vector(gTruthVec, gTruthHandle);
    if (gTruthVec.size() == 0)
    {
      std::cout << "\n[NUMUSEL] No GTruth Information" << std::endl;
      return;
    }

    art::Handle<std::vector<simb::MCFlux>> mcFluxHandle;
    evt.getByLabel(_mctruthLabel, mcFluxHandle);
    if (!mcFluxHandle.isValid())
      return;
    std::vector<art::Ptr<simb::MCFlux>> mcFluxVec;
    art::fill_ptr_vector(mcFluxVec, mcFluxHandle);
    if (mcFluxVec.size() == 0)
    {
      std::cout << "\n[NUMUSEL] No MCFlux Information" << std::endl;
      return;
    }

    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    evt.getByLabel(_mctruthLabel, mcTruthHandle);
    if (!mcTruthHandle.isValid())
      return;
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVec;
    art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
    if (mcTruthVec.size() == 0)
    {
      std::cout << "\n[NUMUSEL] No MCTruth Information" << std::endl;
      return;
    }

    const art::Ptr<simb::MCFlux> mcFlux = mcFluxVec.at(0);
    const art::Ptr<simb::GTruth> gTruth = gTruthVec.at(0);
    const art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);

    ewutil.WriteTree(evt, mcFlux, mcTruth, gTruth);
  }

  _flash_PE = fElectronEventSelectionAlg.get_flash_PE();
  _flash_time = fElectronEventSelectionAlg.get_flash_time();

  _flash_x = fElectronEventSelectionAlg.get_flash_x();
  _TPC_x = fElectronEventSelectionAlg.get_TPC_x();
  _category = 0;
  std::vector<double> true_neutrino_vertex(3);
  std::cout << "[PandoraLEEAnalyzer] Is real data? " << evt.isRealData() << std::endl;

  art::Handle<std::vector<ubana::SelectionResult>> selection_h;
  evt.getByLabel("UBXSec", selection_h);

  if (!selection_h.isValid() || selection_h->empty())
  {
    std::cout << "[PandoraLEEAnalyzer] SelectionResult handle is not valid or empty." << std::endl;
  }

  std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
  if (selection_h.isValid())
    art::fill_ptr_vector(selection_v, selection_h);

  if (selection_v.size() > 0) {
    _numu_passed = (int)selection_v.at(0)->GetSelectionStatus();
    if (selection_v.at(0)->GetSelectionStatus())
    {
        std::cout << "[PandoraLEEAnalyzer] Event is selected by UBXSec" << std::endl;
    }
    else
    {
        std::cout << "[PandoraLEEAnalyzer] Event is not selected by UBXSec" << std::endl;
        std::cout << "[PandoraLEEAnalyzer] Failure reason " << selection_v.at(0)->GetFailureReason() << std::endl;
    }
    std::map<std::string, bool> failure_map = selection_v.at(0)->GetCutFlowStatus();
    for (auto iter : failure_map)
    {
      if (m_printDebug) {
          std::cout << "[PandoraLEEAnalyzer] UBXSec Cut: " << iter.first << "  >>>  " << (iter.second ? "PASSED" : "NOT PASSED") << std::endl;
      }
      if (iter.second) {
          _numu_cuts += 1;
      }
    }
  }



  if ((!evt.isRealData() || m_isOverlaidSample) && !m_isCosmicInTime)
  {
    // nu_e flux must be corrected by event weight
    art::InputTag eventweight_tag("eventweight");
    auto const &eventweights_handle =
        evt.getValidHandle<std::vector<evwgh::MCEventWeight>>(eventweight_tag);
    if (!eventweights_handle.isValid()) {
      std::cout << "[PandoraLEEAnalyzer] No MCEventWeight data product" << std::endl;
      _bnbweight = 1;
    } else {
      auto const &eventweights(*eventweights_handle);
      if (eventweights.size() > 0) {
        for (auto last : eventweights.at(0).fWeight) {
          if (last.first.find("bnbcorrection") != std::string::npos && std::isfinite(last.second.at(0))) {
              _bnbweight = last.second.at(0);
          } else {
            _bnbweight = 1;
          }
        }
      } else {
        _bnbweight = 1;
      }
    }

    try {
      art::InputTag genie_eventweight_tag("genieeventweightmultisim");
      auto const &genie_eventweights_handle = evt.getValidHandle<std::vector<evwgh::MCEventWeight>>(genie_eventweight_tag);
      if (!genie_eventweights_handle.isValid())
      {
        std::cout << "[PandoraLEEAnalyzer] GENIE MCEventWeight handle not valid" << std::endl;
      }
      else
      {
        auto const &genie_eventweights(*genie_eventweights_handle);
        _weights = genie_eventweights.at(0).fWeight;
      }
    } catch (...) {
      std::cout << "[PandoraLEEAnalyzer] No GENIE MCEventWeight data product" << std::endl;
    }

    try
    {
      art::InputTag flux_eventweight_tag("fluxeventweightmultisim");
      auto const &flux_eventweights_handle = evt.getValidHandle<std::vector<evwgh::MCEventWeight>>(flux_eventweight_tag);
      if (!flux_eventweights_handle.isValid())
      {
        std::cout << "[PandoraLEEAnalyzer] Flux MCEventWeight handle not valid" << std::endl;
      }
      else
      {
        auto const &flux_eventweights(*flux_eventweights_handle);
        _flux_weights = flux_eventweights.at(0).fWeight;
      }
    }
    catch (...)
    {
      std::cout << "[PandoraLEEAnalyzer] No Flux MCEventWeight data product" << std::endl;
    }

    auto const &generator_handle = evt.getValidHandle<std::vector<simb::MCTruth>>(_mctruthLabel);
    auto const &generator(*generator_handle);
    _n_true_nu = generator.size();
    _true_nu_is_fiducial = 0;

    bool there_is_a_neutrino = false;

    for (auto &gen : generator)
    {
      if (gen.Origin() == simb::kBeamNeutrino)
      {
        there_is_a_neutrino = true;
        _nu_pdg = gen.GetNeutrino().Nu().PdgCode();
        _nu_energy = gen.GetNeutrino().Nu().E();
        _nu_theta = gen.GetNeutrino().Nu().Theta();
        _nu_phi = gen.GetNeutrino().Nu().Phi();

        if (abs(_nu_pdg) == 12) {
            int n_bin = _h_lee_scaling->FindBin(_nu_energy*1000);
            _leeweight = _lee_scaling[n_bin];
            std::cout << "LEE WEIGHT " << _nu_energy << " " << _leeweight << std::endl;
        }

        _ccnc = gen.GetNeutrino().CCNC();

        if (_ccnc == simb::kNC)
        {
          _category = k_nc;
        }

        true_neutrino_vertex[0] = gen.GetNeutrino().Nu().Vx();
        true_neutrino_vertex[1] = gen.GetNeutrino().Nu().Vy();
        true_neutrino_vertex[2] = gen.GetNeutrino().Nu().Vz();
        _true_vx = true_neutrino_vertex[0];
        _true_vy = true_neutrino_vertex[1];
        _true_vz = true_neutrino_vertex[2];
        _true_nu_is_fiducial = (int)geoHelper.isFiducial(true_neutrino_vertex);

        _interaction_type = gen.GetNeutrino().Mode();

        if (sce_service->GetPosOffsets(_true_vx, _true_vy, _true_vz).size() == 3)
        {
          _true_vx_sce =
              _true_vx - sce_service->GetPosOffsets(_true_vx, _true_vy, _true_vz)[0] + 0.7;
          _true_vy_sce =
              _true_vy + sce_service->GetPosOffsets(_true_vx, _true_vy, _true_vz)[1];
          _true_vz_sce =
              _true_vz + sce_service->GetPosOffsets(_true_vx, _true_vy, _true_vz)[2];
        }
        else
        {
          std::cout << "[PandoraLEEAnalyzer] "
                    << "Space Charge service offset size < 3" << std::endl;
          continue;
        }

        if (!geoHelper.isActive(true_neutrino_vertex))
        {
          _category = k_dirt;
        }
        break; // In case of events with more than one neutrino (2% of the total) we take for the moment only the first one
      }
    }

    if (!there_is_a_neutrino)
      _category = k_cosmic;

    auto const &mcparticles_handle = evt.getValidHandle< std::vector< simb::MCParticle > >("largeant");
    auto const &mcparticles(*mcparticles_handle);

    for (auto &mcparticle : mcparticles)
    {
      if (!(mcparticle.Process() == "primary" &&
            mcparticle.T() != 0 &&
            mcparticle.StatusCode() == 1))
        continue;

      const auto mc_truth =
            pandoraHelper.TrackIDToMCTruth(evt, "largeant", mcparticle.TrackId());
      if (mc_truth->Origin() == simb::kBeamNeutrino)
      {
        _nu_daughters_E.push_back(mcparticle.E());
        _nu_daughters_pdg.push_back(mcparticle.PdgCode());

        if (abs(mcparticle.PdgCode()) == 211)
        {
          _n_true_pions += 1;
        }
        else if (abs(mcparticle.PdgCode()) == 2212)
        {
          _n_true_protons += 1;
        }
        _nu_daughters_px.push_back(mcparticle.Px());
        _nu_daughters_py.push_back(mcparticle.Py());
        _nu_daughters_pz.push_back(mcparticle.Pz());

        _nu_daughters_vx.push_back(mcparticle.Vx());
        _nu_daughters_vy.push_back(mcparticle.Vy());
        _nu_daughters_vz.push_back(mcparticle.Vz());

        _nu_daughters_endx.push_back(mcparticle.EndX());
        _nu_daughters_endy.push_back(mcparticle.EndY());
        _nu_daughters_endz.push_back(mcparticle.EndZ());
      }
    }

    //Insert block to save the start point of the MCshower object for all showers that have a neutrino as mother and a kbeamneutrino as origin
    auto const &mcshower_handle = evt.getValidHandle<std::vector<sim::MCShower>>("mcreco");
    auto const &mcshowers(*mcshower_handle);

    for (auto const &mcshower: mcshowers)
    {
      int pdg_mother = mcshower.MotherPdgCode();
      int origin = mcshower.Origin();

      if ((pdg_mother == 22 || pdg_mother == 11) && origin == 1)
      {
        _true_shower_pdg.push_back(mcshower.AncestorPdgCode());
        _true_shower_depE.push_back(mcshower.DetProfile().E());

        double x_det = mcshower.Start().X();
        double y_det = mcshower.Start().Y();
        double z_det = mcshower.Start().Z();

        if (pdg_mother == 22)
        { //For photons take the end of the shower
          x_det = mcshower.End().X();
          y_det = mcshower.End().Y();
          z_det = mcshower.End().Z();
        }

        _true_shower_x_sce.push_back(x_det - sce_service->GetPosOffsets(x_det, y_det, z_det)[0] + 0.7);
        _true_shower_y_sce.push_back(y_det + sce_service->GetPosOffsets(x_det, y_det, z_det)[1]);
        _true_shower_z_sce.push_back(z_det + sce_service->GetPosOffsets(x_det, y_det, z_det)[2]);

        if (m_printDebug) {
          std::cout << "[PandoraLEEAnalyzer] "
                    << "MCShower End: (" << x_det - sce_service->GetPosOffsets(x_det, y_det, z_det)[0] + 0.7
                    << "," << y_det + sce_service->GetPosOffsets(x_det, y_det, z_det)[1]
                    << "," << z_det + sce_service->GetPosOffsets(x_det, y_det, z_det)[2] << ")" << std::endl;

          std::cout << "[PandoraLEEAnalyzer] "
                    << "TrueVTX: (" << _true_vx_sce << "," << _true_vy_sce << "," << _true_vz_sce << ")" << std::endl;
        }
      }
    }

    if (_category != k_cosmic && _category != k_dirt && _category != k_nc)
    {
      if (abs(_nu_pdg) == 12)
      {
        _category = k_nu_e;
      }
      if (abs(_nu_pdg) == 14)
      {
        _category = k_nu_mu;
      }
    }

    std::cout << "[PandoraLEEAnalyzer] "
            << "True neutrino PDG " << _nu_pdg << std::endl;
    std::cout << "[PandoraLEEAnalyzer] "
              << "Nu energy " << _nu_energy << std::endl;
    std::cout << "[PandoraLEEAnalyzer] "
              << "Nu theta " << _nu_theta << std::endl;
    std::cout << "[PandoraLEEAnalyzer] "
              << "Nu phi " << _nu_phi << std::endl;
  }
  else
  {
    _category = k_data;
  }


  std::vector<art::Ptr<recob::PFParticle>> neutrino_pf;
  std::vector<art::Ptr<recob::PFParticle>> cosmic_pf;

  std::vector<int> neutrino_pdg;
  std::vector<std::string> neutrino_process;
  std::vector<double> neutrino_energy;

  std::vector<int> cosmic_pdg;
  std::vector<std::string> cosmic_process;
  std::vector<double> cosmic_energy;

  if ((!evt.isRealData() || m_isOverlaidSample) && !m_isCosmicInTime)
  {

    categorizePFParticles(evt,
                          neutrino_pdg, neutrino_process, neutrino_energy, neutrino_pf,
                          cosmic_pdg, cosmic_process, cosmic_energy, cosmic_pf);

    _n_matched = neutrino_pf.size();
  }

  // save some information of all the neutrino n_candidates
  _n_total_candidates = fElectronEventSelectionAlg.get_n_neutrino_candidates();
  for (auto &inu : fElectronEventSelectionAlg.get_primary_indexes())
  {
    TVector3 _candidate_v = fElectronEventSelectionAlg.get_neutrino_vertex().at(inu);
    _candidate_vx.push_back(_candidate_v.X());
    _candidate_vy.push_back(_candidate_v.Y());
    _candidate_vz.push_back(_candidate_v.Z());
  }

  fElectronEventSelectionAlg.get_neutrino_vertex();

  // save some information of the passed candidates
  size_t ipf_candidate = std::numeric_limits<size_t>::lowest();
  if (_event_passed)
  {

    for (auto &inu : fElectronEventSelectionAlg.get_primary_indexes())
    {
      if (fElectronEventSelectionAlg.get_neutrino_candidate_passed().at(inu))
      {
        nu_candidates.push_back(inu);
      }
    }

    _n_candidates = nu_candidates.size();

    std::cout << "[PandoraLEE] "
              << "EVENT PASSED" << std::endl;

    auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
    auto const &track_handle = evt.getValidHandle<std::vector<recob::Track>>(m_pfp_producer);
    auto const &spcpnts_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);
    auto const &cluster_handle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pfp_producer);

    art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
    art::FindManyP<recob::Hit> hits_per_spcpnts(spcpnts_handle, evt, m_pfp_producer);
    art::FindManyP<anab::Calorimetry> calos_per_track(track_handle, evt, m_calorimetry_producer);
    art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
    art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, m_pfp_producer);
    art::FindManyP<anab::ParticleID> pid_per_track(track_handle, evt, m_pid_producer);
    art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

    ipf_candidate = choose_candidate(nu_candidates, evt);
    std::cout << "[PandoraLEE] "
              << "Neutrino candidate " << ipf_candidate << std::endl;
    _chosen_candidate = ipf_candidate;

    auto const &vertex_obj = vertex_per_pfpart.at(ipf_candidate);

    double reco_neutrino_vertex[3];
    vertex_obj->XYZ(reco_neutrino_vertex);
    _vx = reco_neutrino_vertex[0];
    _vy = reco_neutrino_vertex[1];
    _vz = reco_neutrino_vertex[2];
    vector<double> reco_vertex = {_vx, _vy, _vz};

    if (_category != k_data)
    {
      TVector3 v_reco_vertex(_vx, _vy, _vz);
      TVector3 sce_true_vertex(_true_vx_sce, _true_vy_sce, _true_vz_sce);
      _distance = geoHelper.distance(v_reco_vertex, sce_true_vertex);
    }

    _n_tracks = fElectronEventSelectionAlg.get_n_tracks().at(ipf_candidate);

    if (_n_tracks > 0)
    {
      art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
      _nu_track_ids = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(ipf_candidate);

      for (auto &pf_id : _nu_track_ids)
      {
        std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pf_id);
        std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pf_id);
        auto const &track_obj = track_per_pfpart.at(pf_id);

        recob::PFParticle const &pfparticle = pfparticle_handle->at(pf_id);
        _nu_track_daughters.push_back(pfparticle.Daughters());

        std::vector<double> dqdx(3, std::numeric_limits<double>::lowest());
        std::vector<double> dedx(3, std::numeric_limits<double>::lowest());

        if (calos_per_track.isValid())
        {
          std::vector<art::Ptr<anab::Calorimetry>> calos = calos_per_track.at(track_obj->ID());
          energyHelper.track_dQdx(&calos, dqdx, dedx);
        } else {
          std::cout << "[PandoraLEEAnalyzer] calos_per_track not valid" << std::endl;
        }

        _track_dQdx.push_back(dqdx);

        std::vector<double> track_cali;

        energyHelper.get_cali(&spcpnts, &hits_per_spcpnts, track_cali);

        _track_energy_cali.push_back(track_cali);
        _track_dEdx.push_back(dedx);

        double mean = std::numeric_limits<double>::lowest();
        double stdev = std::numeric_limits<double>::lowest();
        energyHelper.cluster_residuals(&clusters, &hits_per_cluster, mean, stdev);
        _track_res_mean.push_back(mean);
        _track_res_std.push_back(stdev);

        _track_n_clusters.push_back(clusters.size());

        if (m_useParticleID) {
          art::Ptr<anab::ParticleID> selected_pid = pid_per_track.at(track_obj->ID())[0];

          double bragg_p = std::max(energyHelper.PID(selected_pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212),
                                    energyHelper.PID(selected_pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212));

          double bragg_mu = std::max(energyHelper.PID(selected_pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13),
                                     energyHelper.PID(selected_pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13));

          double bragg_mip = energyHelper.PID(selected_pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0);

          double pidchipr = energyHelper.PID(selected_pid, "Chi2", anab::kGOF, anab::kForward, 2212);
          double pidchimu = energyHelper.PID(selected_pid, "Chi2", anab::kGOF, anab::kForward, 13);
          double pidchipi = energyHelper.PID(selected_pid, "Chi2", anab::kGOF, anab::kForward, 211);
          double pidchika = energyHelper.PID(selected_pid, "Chi2", anab::kGOF, anab::kForward, 321);

          double pida_mean = energyHelper.PID(selected_pid, "PIDA_mean", anab::kPIDA, anab::kForward, 0);

          _track_bragg_p.push_back(bragg_p);
          _track_bragg_mu.push_back(bragg_mu);
          _track_bragg_mip.push_back(bragg_mip);
          _track_pida.push_back(pida_mean);
          _track_pidchipr.push_back(pidchipr);
          _track_pidchimu.push_back(pidchimu);
          _track_pidchipi.push_back(pidchipi);
          _track_pidchika.push_back(pidchika);
        }
        _matched_tracks.push_back(std::numeric_limits<int>::lowest());
        _matched_tracks_process.push_back("");
        _matched_tracks_energy.push_back(std::numeric_limits<double>::lowest());

        std::vector<double> start_point;
        std::vector<double> end_point;

        start_point.push_back(track_obj->Start().X());
        start_point.push_back(track_obj->Start().Y());
        start_point.push_back(track_obj->Start().Z());

        end_point.push_back(track_obj->End().X());
        end_point.push_back(track_obj->End().Y());
        end_point.push_back(track_obj->End().Z());

        _track_is_fiducial.push_back((int)(geoHelper.isFiducial(start_point) &&
                                           geoHelper.isFiducial(end_point)));

        std::vector<double> this_energy;
        std::vector<int> this_nhits;

        energyHelper.energy_from_hits(&clusters, &hits_per_cluster, this_nhits, this_energy);

        _track_energy_hits.push_back(this_energy);

        _track_length.push_back(track_obj->Length());
        _track_id.push_back(track_obj->ID());

        _track_distance.push_back(geoHelper.distance(start_point, reco_vertex));
        _track_dir_x.push_back(track_obj->StartDirection().X());
        _track_dir_y.push_back(track_obj->StartDirection().Y());
        _track_dir_z.push_back(track_obj->StartDirection().Z());

        _track_start_x.push_back(track_obj->Start().X());
        _track_start_y.push_back(track_obj->Start().Y());
        _track_start_z.push_back(track_obj->Start().Z());

        _track_end_x.push_back(track_obj->End().X());
        _track_end_y.push_back(track_obj->End().Y());
        _track_end_z.push_back(track_obj->End().Z());

        _track_theta.push_back(track_obj->Theta());
        _track_phi.push_back(track_obj->Phi());

        std::vector<std::vector<double>> pca;
        pca.resize(2, std::vector<double>(3));

        energyHelper.PCA(&clusters, &hits_per_cluster, pca);

        _track_pca.push_back(pca[0]);
        _track_nhits.push_back(this_nhits);
      }
    }

    try
    {
      _nu_shower_ids = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary().at(
              ipf_candidate);
      _n_showers = fElectronEventSelectionAlg.get_n_showers().at(ipf_candidate);
    }
    catch (...)
    {
      std::cout << "[PandoraLEE] Error getting associated showers to neutrino "
                   "candidate"
                << std::endl;
    }

    art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

    for (auto &pf_id : _nu_shower_ids)
    {
      auto const &shower_obj = shower_per_pfpart.at(pf_id);
      if (shower_obj.isNull()) {
        std::cout << "[PandoraLEEAnalyzer] Shower pointer " << pf_id << " is null, exiting" << std::endl;
        continue;
      }

      std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pf_id);

      recob::PFParticle const &pfparticle = pfparticle_handle->at(pf_id);
      _nu_shower_daughters.push_back(pfparticle.Daughters());

      double mean = std::numeric_limits<double>::lowest();
      double stdev = std::numeric_limits<double>::lowest();
      energyHelper.cluster_residuals(&clusters, &hits_per_cluster, mean, stdev);
      _shower_res_mean.push_back(mean);
      _shower_res_std.push_back(stdev);

      std::vector<double> pitches = {0, 0, 0);
      std::vector<double> dQdx_hits_in_the_box = {-1, -1, -1);

      std::vector<double> dqdx(3, std::numeric_limits<double>::lowest());
      std::vector<double> dedx(3, std::numeric_limits<double>::lowest());
      std::vector<std::vector<double>> dqdx_hits_shower(3, std::vector<double>());

      std::vector<double> dqdx_cali(3, std::numeric_limits<double>::lowest());

      _matched_showers.push_back(std::numeric_limits<int>::lowest());
      _matched_showers_process.push_back("");

      _matched_showers_energy.push_back(std::numeric_limits<double>::lowest());
      energyHelper.dQdx(&(*shower_obj),
                        &clusters,
                        &hits_per_cluster,
                        dqdx,
                        dqdx_hits_shower,
                        pitches,
                        dQdx_hits_in_the_box);
      energyHelper.dQdx_cali(&(*shower_obj), dqdx_cali);

      _shower_dQdx_hits.push_back(dqdx_hits_shower[2]);
      _shower_pitches.push_back(pitches);
      _shower_n_clusters.push_back(clusters.size());

      std::vector<std::vector<double>> dedx_hits_shower(3, std::vector<double>());

      energyHelper.dEdx_from_dQdx(dedx_hits_shower[0], dqdx_hits_shower[0]);
      energyHelper.dEdx_from_dQdx(dedx_hits_shower[1], dqdx_hits_shower[1]);
      energyHelper.dEdx_from_dQdx(dedx_hits_shower[2], dqdx_hits_shower[2]);

      _shower_dEdx_hits.push_back(dedx_hits_shower[2]);

      for (size_t i_pl=0; i_pl < 3; i_pl++) {
        if (dedx_hits_shower[i_pl].size()) {
            std::nth_element(dedx_hits_shower[i_pl].begin(), dedx_hits_shower[i_pl].begin() + dedx_hits_shower[i_pl].size() / 2, dedx_hits_shower[i_pl].end());
            dedx[i_pl] = dedx_hits_shower[i_pl][dedx_hits_shower[i_pl].size() / 2];
        }
      }

      _shower_dQdx.push_back(dqdx);
      _shower_dEdx.push_back(dedx);

      _shower_dQdx_hits_in_the_box.push_back(dQdx_hits_in_the_box);
      _shower_dQdx_cali.push_back(dqdx_cali);

      _shower_dir_x.push_back(shower_obj->Direction().X());
      _shower_dir_y.push_back(shower_obj->Direction().Y());
      _shower_dir_z.push_back(shower_obj->Direction().Z());

      _shower_open_angle.push_back(shower_obj->OpenAngle());
      double shower_length = shower_obj->Length();
      _shower_length.push_back(shower_length);

      std::vector<double> start_point;
      std::vector<double> end_point;

      start_point.resize(3);
      end_point.resize(3);
      for (int ix = 0; ix < 3; ix++)
      {
        start_point[ix] = shower_obj->ShowerStart()[ix];
        end_point[ix] =
            shower_obj->ShowerStart()[ix] + shower_obj->Direction()(ix) * shower_length;
      }

      _shower_is_fiducial.push_back((int)(geoHelper.isFiducial(start_point) &&
                                          geoHelper.isFiducial(end_point)));

      _shower_distance.push_back(geoHelper.distance(start_point, reco_vertex));
      _shower_start_x.push_back(shower_obj->ShowerStart().X());
      _shower_start_y.push_back(shower_obj->ShowerStart().Y());
      _shower_start_z.push_back(shower_obj->ShowerStart().Z());

      _shower_phi.push_back(shower_obj->Direction().Phi());
      _shower_theta.push_back(shower_obj->Direction().Theta());

      std::vector<double> this_energy;
      std::vector<int> this_nhits;

      energyHelper.energy_from_hits(&clusters, &hits_per_cluster, this_nhits, this_energy);
      std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pf_id);

      std::vector<double> shower_cali;
      energyHelper.get_cali(&spcpnts, &hits_per_spcpnts, shower_cali);

      _shower_energy_cali.push_back(shower_cali);
      _shower_energy.push_back(this_energy);

      std::vector<std::vector<double>> pca;
      pca.resize(2, std::vector<double>(3));

      energyHelper.PCA(&clusters, &hits_per_cluster, pca);

      _shower_pca.push_back(pca[0]);
      _shower_nhits.push_back(this_nhits);

    }

    /* For each neutrino PFParticle checks how many daughter tracks and showers
    with true neutrino origin we have */
    _nu_matched_showers = 0;
    bool shower_cr_found = false;

    _nu_matched_tracks = 0;
    bool track_cr_found = false;
    double cr_hits = 0;
    double nu_hits = 0;

    for (auto &inu : fElectronEventSelectionAlg.get_primary_indexes())
    {
      _primary_indexes.push_back(inu);

      _flash_passed.push_back(fElectronEventSelectionAlg.get_op_flash_indexes().at(inu));
      _number_tracks.push_back(fElectronEventSelectionAlg.get_n_tracks().at(inu));
      _number_showers.push_back(fElectronEventSelectionAlg.get_n_showers().at(inu));

      std::vector<size_t> pfp_showers_id = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary().at(inu);

      int pass_shower = 0;

      for (size_t ish = 0; ish < pfp_showers_id.size(); ish++)
      {
        for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
        {

          if (pfp_showers_id[ish] == neutrino_pf[ipf].key())
          {
            pass_shower += 1;

            if (inu == ipf_candidate)
            {
              _nu_matched_showers++;
              nu_hits += _shower_nhits[ish][2] + _shower_nhits[ish][1] + _shower_nhits[ish][0];
              _matched_showers[ish] = neutrino_pdg[ipf];
              _matched_showers_process[ish] = neutrino_process[ipf];
              _matched_showers_energy[ish] = neutrino_energy[ipf];
            }
          }
        }

        if (inu == ipf_candidate)
        {
          for (size_t ipf = 0; ipf < cosmic_pf.size(); ipf++)
          {
            if (pfp_showers_id[ish] == cosmic_pf[ipf].key())
            {
              shower_cr_found = true;
              cr_hits += _shower_nhits[ish][2] + _shower_nhits[ish][1] + _shower_nhits[ish][0];
              _matched_showers[ish] = cosmic_pdg[ipf];
              _matched_showers_process[ish] = cosmic_process[ipf];
              _matched_showers_energy[ish] = cosmic_energy[ipf];
            }
          }
        }
      }

      _shower_passed.push_back(pass_shower);

      if (_number_tracks.back() > 0)
      {
        std::vector<size_t> pfp_tracks_id = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(inu);

        int pass_track = 0;

        for (size_t itr = 0; itr < pfp_tracks_id.size(); itr++)
        {

          for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
          {
            if (pfp_tracks_id[itr] == neutrino_pf[ipf].key())
            {
              pass_track += 1;
              if (inu == ipf_candidate)
              {
                _nu_matched_tracks++;
                nu_hits += _track_nhits[itr][2] + _track_nhits[itr][1] + _track_nhits[itr][0];
                _matched_tracks[itr] = neutrino_pdg[ipf];
                _matched_tracks_process[itr] = neutrino_process[ipf];
                _matched_tracks_energy[itr] = neutrino_energy[ipf];
              }
            }
          }

          if (inu == ipf_candidate)
          {
            for (size_t ipf = 0; ipf < cosmic_pf.size(); ipf++)
            {
              if (pfp_tracks_id[itr] == cosmic_pf[ipf].key())
              {
                track_cr_found = true;
                cr_hits += _track_nhits[itr][2] + _track_nhits[itr][1] + _track_nhits[itr][0];
                _matched_tracks[itr] = cosmic_pdg[ipf];
                _matched_tracks_process[itr] = cosmic_process[ipf];
                _matched_tracks_energy[itr] = cosmic_energy[ipf];
              }
            }
          }
        }

        _track_passed.push_back(pass_track);
      }
    }

    _n_primaries = _primary_indexes.size();
    _cosmic_fraction = cr_hits / (cr_hits + nu_hits);

    if (
        !track_cr_found && _nu_matched_tracks == 0 && !shower_cr_found && _nu_matched_showers == 0 && _category != k_dirt && _category != k_data)
    {
      _category = k_other;
      std::cout << "[PandoraLEE] "
                << "***NOT NEUTRINO NOR COSMIC***" << std::endl;
    }

    if ((track_cr_found || shower_cr_found) &&
        (_nu_matched_tracks > 0 || _nu_matched_showers > 0))
    {
      _category = k_mixed;
    }

    if ((track_cr_found || shower_cr_found) && _category != k_mixed)
    {
      _category = k_cosmic;
    }

    _re_category = _category;
    if (abs(_nu_pdg) == 12)
    {
      if (_n_true_protons>0 && _n_true_pions==0)
      {
        if (_category==2) _re_category = 11;
        else if (_category==7) _re_category = 14;
      }
      else
      {
        if (_category==2) _re_category = 12;
        else if (_category==7) _re_category = 15;
      }
    }
    else if (abs(_nu_pdg) == 14)
    {
      if (_category==7) _re_category = 16;
    }

    std::cout << "[PandoraLEE] "
              << "Category " << _category << std::endl;
  }
  else
  {
    std::cout << "[PandoraLEE] "
              << "EVENT NOT PASSED" << std::endl;

    _n_candidates = 0;
    for (auto &inu : fElectronEventSelectionAlg.get_primary_indexes())
    {
      _primary_indexes.push_back(inu);

      _flash_passed.push_back(fElectronEventSelectionAlg.get_op_flash_indexes().at(inu));
      _number_tracks.push_back(fElectronEventSelectionAlg.get_n_tracks().at(inu));
      _number_showers.push_back(fElectronEventSelectionAlg.get_n_showers().at(inu));

      std::vector<size_t> pfp_showers_id = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary().at(inu);

      int pass_shower = 0;

      for (size_t ish = 0; ish < pfp_showers_id.size(); ish++)
      {
        for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
        {
          if (pfp_showers_id[ish] == neutrino_pf[ipf].key())
          {
            pass_shower += 1;
          }
        }
      }

      _shower_passed.push_back(pass_shower);

      if (_number_tracks.back() > 0)
      {
        std::vector<size_t> pfp_tracks_id = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(inu);
        int pass_track = 0;

        for (size_t itr = 0; itr < pfp_tracks_id.size(); itr++)
        {

          for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
          {
            if (pfp_tracks_id[itr] == neutrino_pf[ipf].key())
            {
              pass_track += 1;
            }
          }
        }

        _track_passed.push_back(pass_track);
      }
    }

    _n_primaries = _primary_indexes.size();
  }

  myTTree->Fill();
  std::cout << "[PandoraLEE] "
            << "END ANALYZER" << std::endl;

} // end analyze function

//------------------------------------------------------------------------------------------------------------------------------------

void lee::PandoraLEEAnalyzer::reconfigure(fhicl::ParameterSet const &pset)
{

  // TODO: add an external fcl file to change configuration
  // add what you want to read, and default values of your labels etc. example:
  //  m_particleLabel = pset.get<std::string>("PFParticleModule","pandoraNu");
  fElectronEventSelectionAlg.reconfigure(pset.get<fhicl::ParameterSet>("ElectronSelectionAlg"));
  energyHelper.reconfigure(pset.get<fhicl::ParameterSet>("EnergyHelper"));

  m_fidvolXstart = pset.get<double>("fidvolXstart", 0);
  m_fidvolXend = pset.get<double>("fidvolXend", 0);

  m_fidvolYstart = pset.get<double>("fidvolYstart", 0);
  m_fidvolYend = pset.get<double>("fidvolYend", 0);

  m_fidvolZstart = pset.get<double>("fidvolZstart", 0);
  m_fidvolZend = pset.get<double>("fidvolZend", 0);

  geoHelper.setFiducialVolumeCuts(m_fidvolXstart, m_fidvolXend, m_fidvolYstart,
                                  m_fidvolYend, m_fidvolZstart, m_fidvolZend);

  m_hitfinderLabel = pset.get<std::string>("HitFinderLabel", "pandoraCosmicHitRemoval::McRecoStage2");
  m_pid_producer = pset.get<std::string>("ParticleIDModuleLabel", "pandoraNucalipid::McRecoCali");
  m_calorimetry_producer = pset.get<std::string>("CalorimetryLabel", "pandoraNucali::McRecoCali");
  m_hitmatching_producer = pset.get<std::string>("HitMatchingLabel", "crHitRemovalTruthMatch::McRecoStage2");
  m_pfp_producer = pset.get<std::string>("PFParticleLabel", "pandoraNu::McRecoStage2");
  m_spacepointLabel = pset.get<std::string>("SpacePointLabel", "pandoraNu::McRecoStage2");
  m_spacepointLabel = pset.get<std::string>("SpacePointLabel", "pandoraNu::McRecoStage2");
  m_printDebug = pset.get<bool>("PrintDebug", false);
  m_useParticleID = pset.get<bool>("UseParticleID", true);

  m_isData = pset.get<bool>("isData", false);
  m_isCosmicInTime = pset.get<bool>("isCosmicInTime", false);
  m_isOverlaidSample = pset.get<bool>("isOverlaidSample", false);
  m_save_flux_info = pset.get<bool>("saveFluxInfo", false);

}

//---------------------------------------------------------------------------------------------------------------------------
// add other functions here

DEFINE_ART_MODULE(lee::PandoraLEEAnalyzer)
