////////////////////////////////////////////////////////////////////////
// Class:       ElectronNeutrinoFilter
// Plugin Type: filter (art v2_05_01)
// File:        ElectronNeutrinoFilter_module.cc
//
// Generated at Tue Jul 31 16:09:58 2018 by Stefano Soleti using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "ElectronNeutrinoFilter.h"

lee::ElectronNeutrinoFilter::ElectronNeutrinoFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  art::ServiceHandle<art::TFileService> tfs;
  myPOTTTree = tfs->make<TTree>("pot", "POT Tree");
  myPOTTTree->Branch("pot", &_pot, "pot/D");
  myPOTTTree->Branch("run", &_run_sr, "run/i");
  myPOTTTree->Branch("subrun", &_subrun_sr, "subrun/i");

  myTTree = tfs->make<TTree>("filtertree", "Filter Tree");

  myTTree->Branch("bnbweight", &_bnbweight, "bnbweight/d");
  // myTTree->Branch("true_neutrino_vertex", "std::vector< double >", &_true_neutrino_vertex);
  // myTTree->Branch("true_neutrino_vertex_sce", "std::vector< double >", &_true_neutrino_vertex_sce);
  myTTree->Branch("passed", &_passed, "passed/O");
  myTTree->Branch("n_total_candidates", &_n_total_candidates, "n_total_candidates/I");

  myTTree->Branch("flash_time", "std::vector< double >", &_flash_time);
  myTTree->Branch("flash_PE", "std::vector< double >", &_flash_PE);
  myTTree->Branch("flash_x", &_flash_x, "flash_x/d");
  myTTree->Branch("TPC_x", &_TPC_x, "TPC_x/d");

  myTTree->Branch("n_true_nu", &_n_true_nu, "n_true_nu/I");
  myTTree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  myTTree->Branch("ccnc", &_ccnc, "ccnc/I");
  myTTree->Branch("interaction_type", &_interaction_type, "interaction_type/I");

  myTTree->Branch("nu_energy", &_nu_energy, "nu_energy/D");
  myTTree->Branch("nu_theta", &_nu_theta, "nu_theta/d");
  myTTree->Branch("nu_phi", &_nu_phi, "nu_phi/d");
  myTTree->Branch("nu_T", &_nu_T, "nu_T/d");
  myTTree->Branch("true_vx", &_true_vx, "true_vx/d");
  myTTree->Branch("true_vy", &_true_vy, "true_vy/d");
  myTTree->Branch("true_vz", &_true_vz, "true_vz/d");
  myTTree->Branch("true_nu_is_fiducial", &_true_nu_is_fiducial, "true_nu_is_fiducial/O");
  myTTree->Branch("true_nu_is_active", &_true_nu_is_active, "true_nu_is_active/O");

  myTTree->Branch("n_true_pions", &_n_true_pions, "n_true_pions/I");
  myTTree->Branch("n_true_protons", &_n_true_protons, "n_true_protons/I");
  myTTree->Branch("n_true_protons_above40", &_n_true_protons_above40, "n_true_protons_above40/I");
  myTTree->Branch("n_true_protons_above21", &_n_true_protons_above21, "n_true_protons_above21/I");
  myTTree->Branch("true_daughter_E", &_true_daughter_E, "true_daughter_E/D");
  myTTree->Branch("true_daughter_theta", &_true_daughter_theta, "true_daughter_theta/D");
  myTTree->Branch("true_daughter_phi", &_true_daughter_phi, "true_daughter_phi/D");
  myTTree->Branch("true_daughter_T", &_true_daughter_T, "true_daughter_T/D");

  _run_subrun_list_file.open("run_subrun_list_filter.txt", std::ofstream::out | std::ofstream::trunc);

  this->reconfigure(p);
}

bool lee::ElectronNeutrinoFilter::endSubRun(art::SubRun &sr)
{
  _run_subrun_list_file << sr.run() << " " << sr.subRun() << std::endl;

  _run_sr = sr.run();
  _subrun_sr = sr.subRun();

  art::Handle<sumdata::POTSummary> potListHandle;
  if (!m_isData || m_isOverlaidSample)
  {
    if (sr.getByLabel(_mctruthLabel, potListHandle)) {
      _pot = potListHandle->totpot;

    } else {
        _pot = 0.;
    }
  }
  else
  {
    if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle)) {
      _pot = potListHandle->totpot;
    } else {
      _pot = 0.;
    }
  }

  myPOTTTree->Fill();
  return true;
}

void lee::ElectronNeutrinoFilter::clear()
{
  // _true_neutrino_vertex.reserve(3);
  // _true_neutrino_vertex_sce.reserve(3);
  _true_neutrino_vertex.clear();
  _true_neutrino_vertex_sce.clear();

  _nu_daughters_p.clear();
  _nu_daughters_start_v.clear();
  _nu_daughters_end_v.clear();

  _nu_daughters_E.clear();
  _nu_daughters_pdg.clear();

  _n_true_nu = std::numeric_limits<unsigned int>::lowest();
  _ccnc = std::numeric_limits<int>::lowest();
  _nu_pdg = std::numeric_limits<int>::lowest();
  _interaction_type = std::numeric_limits<int>::lowest();
  _n_total_candidates = std::numeric_limits<int>::lowest();

  _true_nu_is_fiducial = false;
  _true_nu_is_active = false;
  _nu_energy = std::numeric_limits<double>::lowest();
  _nu_theta = std::numeric_limits<double>::lowest();
  _nu_phi = std::numeric_limits<double>::lowest();
  _nu_T = std::numeric_limits<double>::lowest();
  _true_vx = std::numeric_limits<double>::lowest();
  _true_vy = std::numeric_limits<double>::lowest();
  _true_vz = std::numeric_limits<double>::lowest();

  _lee_weight = 0;
  _passed = false;

  _bnbweight = 1;

  _n_true_pions = 0;
  _n_true_protons = 0;
  _n_true_protons_above40 = 0;
  _n_true_protons_above21 = 0;
  _true_daughter_E = std::numeric_limits<double>::lowest();
  _true_daughter_theta = std::numeric_limits<double>::lowest();
  _true_daughter_phi = std::numeric_limits<double>::lowest();
  _true_daughter_T = std::numeric_limits<double>::lowest();

  _pot = std::numeric_limits<double>::lowest();
  _run_sr = std::numeric_limits<unsigned int>::lowest();
  _subrun_sr = std::numeric_limits<unsigned int>::lowest();
}

bool lee::ElectronNeutrinoFilter::filter(art::Event &e)
{
  clear();
  std::cout << "[ElectronNeutrinoFilter] "
            << "RUN " << e.run() << " SUBRUN " << e.subRun() << " EVENT " << e.id().event()
            << std::endl;

  _passed = fElectronEventSelectionAlg.eventSelected(e);

  _flash_PE = fElectronEventSelectionAlg.get_flash_PE();
  _flash_time = fElectronEventSelectionAlg.get_flash_time();

  _flash_x = fElectronEventSelectionAlg.get_flash_x();
  _TPC_x = fElectronEventSelectionAlg.get_TPC_x();

  _n_total_candidates = fElectronEventSelectionAlg.get_n_neutrino_candidates();

  std::cout << "[ElectronNeutrinoFilter] Passing filter? " << _passed << std::endl;

  bool is_data = e.isRealData();
  if (m_isOverlaidSample) {
    is_data = false;
  }

  if (is_data) {
    myTTree->Fill();
    return _passed;
  }

  if (!e.isRealData() || m_isOverlaidSample)
  {
    // nu_e flux must be corrected by event weight
    art::InputTag eventweight_tag("eventweight");
    auto const &eventweights_handle =
        e.getValidHandle<std::vector<evwgh::MCEventWeight>>(eventweight_tag);
    if (!eventweights_handle.isValid())
    {
      std::cout << "[PandoraLEEAnalyzer] No MCEventWeight data product" << std::endl;
      _bnbweight = 1;
    }
    else
    {
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
  }

  auto const &generator_handle = e.getValidHandle<std::vector<simb::MCTruth>>(_mctruthLabel);
  auto const &generator(*generator_handle);
  _n_true_nu = generator.size();
  auto const *sce_service = lar::providerFrom<spacecharge::SpaceChargeService>();

  for (auto &gen : generator)
  {
    if (gen.Origin() == simb::kBeamNeutrino)
    {
      _nu_pdg = gen.GetNeutrino().Nu().PdgCode();
      _nu_energy = gen.GetNeutrino().Nu().E();
      _nu_theta = gen.GetNeutrino().Nu().Momentum().Theta();
      _nu_phi = gen.GetNeutrino().Nu().Momentum().Phi();
      _nu_T = gen.GetNeutrino().Nu().T();

      if (abs(_nu_pdg) == 12)
      {
        int n_bin = _h_lee_scaling->FindBin(_nu_energy * 1000);
        _lee_weight = _lee_scaling[n_bin];
      }

      _ccnc = gen.GetNeutrino().CCNC();

      _true_neutrino_vertex.push_back(gen.GetNeutrino().Nu().Vx());
      _true_neutrino_vertex.push_back(gen.GetNeutrino().Nu().Vy());
      _true_neutrino_vertex.push_back(gen.GetNeutrino().Nu().Vz());
      _true_vx = _true_neutrino_vertex[0];
      _true_vy = _true_neutrino_vertex[1];
      _true_vz = _true_neutrino_vertex[2];
      _true_nu_is_fiducial = geoHelper.isFiducial(_true_neutrino_vertex);
      _true_nu_is_active = geoHelper.isActive(_true_neutrino_vertex);
      _interaction_type = gen.GetNeutrino().Mode();

      std::vector<double> sce_offsets = sce_service->GetPosOffsets(_true_neutrino_vertex[0], _true_neutrino_vertex[1], _true_neutrino_vertex[2]);

      if (sce_offsets.size() == 3)
      {
        std::cout << _true_neutrino_vertex.size() << std::endl;

        _true_neutrino_vertex_sce.push_back(_true_neutrino_vertex[0] - sce_offsets[0] + 0.7);
        _true_neutrino_vertex_sce.push_back(_true_neutrino_vertex[1] + sce_offsets[1]);
        _true_neutrino_vertex_sce.push_back(_true_neutrino_vertex[2] + sce_offsets[2]);
      }
      else
      {
        std::cout << "[PandoraLEEAnalyzer] "
                  << "Space Charge service offset size < 3" << std::endl;
        continue;
      }
      break; // In case of events with more than one neutrino (2% of the total) we take for the moment only the first one
    }
  }

  auto const &mcparticles_handle = e.getValidHandle<std::vector<simb::MCParticle>>(_mcparticleLabel);
  auto const &mcparticles(*mcparticles_handle);

  for (auto &mcparticle : mcparticles)
  {
    if (!(mcparticle.Process() == "primary" &&
          mcparticle.T() != 0 &&
          mcparticle.StatusCode() == 1))
      continue;

    const auto mc_truth = pandoraHelper.TrackIDToMCTruth(e, _mcparticleLabel, mcparticle.TrackId());
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
        if (mcparticle.E() > (0.938272 + 0.040))
        {
          _n_true_protons_above40 += 1;
        }
        if (mcparticle.E() > (0.938272 + 0.02108))
        {
          _n_true_protons_above21 += 1;
        }
      }

      std::vector<double> p;
      p.push_back(mcparticle.Px());
      p.push_back(mcparticle.Py());
      p.push_back(mcparticle.Pz());

      _nu_daughters_p.push_back(p);

      std::vector<double> start_v;
      start_v.push_back(mcparticle.Vx());
      start_v.push_back(mcparticle.Vy());
      start_v.push_back(mcparticle.Vz());

      _nu_daughters_start_v.push_back(start_v);

      std::vector<double> end_v;
      end_v.push_back(mcparticle.EndX());
      end_v.push_back(mcparticle.EndY());
      end_v.push_back(mcparticle.EndZ());

      _nu_daughters_end_v.push_back(end_v);

      if (_ccnc == 0)
      {
        if (mcparticle.PdgCode() == _pdg_daughter[_nu_pdg])
        {
          _true_daughter_E = mcparticle.E();
          _true_daughter_theta = mcparticle.Momentum().Theta();
          _true_daughter_phi = mcparticle.Momentum().Phi();
          _true_daughter_T = mcparticle.T();
        }
      }

    }
  }

  myTTree->Fill();

  return _passed;
}

void lee::ElectronNeutrinoFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fElectronEventSelectionAlg.reconfigure(p.get<fhicl::ParameterSet>("ElectronSelectionAlg"));
  m_isOverlaidSample = p.get<bool>("isOverlaidSample", false);
  m_isData = p.get<bool>("isData", false);

  m_fidvolXstart = p.get<double>("fidvolXstart", 0);
  m_fidvolXend = p.get<double>("fidvolXend", 0);

  m_fidvolYstart = p.get<double>("fidvolYstart", 0);
  m_fidvolYend = p.get<double>("fidvolYend", 0);

  m_fidvolZstart = p.get<double>("fidvolZstart", 0);
  m_fidvolZend = p.get<double>("fidvolZend", 0);

  geoHelper.setFiducialVolumeCuts(m_fidvolXstart, m_fidvolXend, m_fidvolYstart,
                                  m_fidvolYend, m_fidvolZstart, m_fidvolZend);
}

DEFINE_ART_MODULE(lee::ElectronNeutrinoFilter)
