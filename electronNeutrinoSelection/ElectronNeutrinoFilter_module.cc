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

  myTTree->Branch("true_neutrino_vertex", "std::vector< double >", &_true_neutrino_vertex);
  myTTree->Branch("true_neutrino_vertex_sce", "std::vector< double >", &_true_neutrino_vertex_sce);
  myTTree->Branch("passed", &_passed, "passed/O");
  myTTree->Branch("nu_energy", &_nu_energy, "nu_energy/D");
  myTTree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");

  _run_subrun_list_file.open("run_subrun_list_filter.txt", std::ofstream::out | std::ofstream::trunc);

  this->reconfigure(p);
}

void lee::ElectronNeutrinoFilter::respondToOpenInputFile(art::FileBlock const &fb)
{
    _sum_pot = 0;
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
      _pot = potListHandle->totpot - _sum_pot;
      if (m_isOverlaidSample) {
          _sum_pot += _pot;
      }
      std::cout << "Subrun POT " << _pot << " " << potListHandle->totpot << " " << _sum_pot << std::endl;
    } else {
        _pot = 0.;
    }
  }
  else
  {
    if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle)) {
      _pot = potListHandle->totpot;
      std::cout << "Subrun POT " << _pot << std::endl;
    } else {
      _pot = 0.;
    }
  }

  myPOTTTree->Fill();
  return true;
}

void lee::ElectronNeutrinoFilter::clear()
{
  _true_neutrino_vertex.reserve(3);
  _true_neutrino_vertex_sce.reserve(3);
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
  _true_nu_is_fiducial = false;
  _nu_energy = std::numeric_limits<double>::lowest();
  _lee_weight = 0;
  _passed = false;
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
  std::cout << "[ElectronNeutrinoFilter] Passing filter? " << _passed << std::endl;

  bool is_data = e.isRealData();
  if (m_isOverlaidSample) {
    is_data = false;
  }

  if (is_data) {
    myTTree->Fill();
    return _passed;
  }

  auto const &generator_handle = e.getValidHandle<std::vector<simb::MCTruth>>(_mctruthLabel);
  auto const &generator(*generator_handle);
  _n_true_nu = generator.size();

  _true_nu_is_fiducial = 0;
  auto const *sce_service = lar::providerFrom<spacecharge::SpaceChargeService>();

  for (auto &gen : generator)
  {
    if (gen.Origin() == simb::kBeamNeutrino)
    {
      _nu_pdg = gen.GetNeutrino().Nu().PdgCode();
      _nu_energy = gen.GetNeutrino().Nu().E();

      if (abs(_nu_pdg) == 12)
      {
        int n_bin = _h_lee_scaling->FindBin(_nu_energy * 1000);
        _lee_weight = _lee_scaling[n_bin];
      }

      _ccnc = gen.GetNeutrino().CCNC();

      _true_neutrino_vertex[0] = gen.GetNeutrino().Nu().Vx();
      _true_neutrino_vertex[1] = gen.GetNeutrino().Nu().Vy();
      _true_neutrino_vertex[2] = gen.GetNeutrino().Nu().Vz();

      _interaction_type = gen.GetNeutrino().Mode();

      std::vector<double> sce_offsets = sce_service->GetPosOffsets(_true_neutrino_vertex[0], _true_neutrino_vertex[1], _true_neutrino_vertex[2]);


      if (sce_offsets.size() == 3)
      {
        _true_neutrino_vertex_sce[0] = _true_neutrino_vertex[0] - sce_offsets[0] + 0.7;
        _true_neutrino_vertex_sce[1] = _true_neutrino_vertex[1] + sce_offsets[1];
        _true_neutrino_vertex_sce[2] = _true_neutrino_vertex[2] + sce_offsets[2];
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
    // if (mc_truth.isNull())
    //   continue;

    if (mc_truth->Origin() == simb::kBeamNeutrino)
    {
      _nu_daughters_E.push_back(mcparticle.E());
      _nu_daughters_pdg.push_back(mcparticle.PdgCode());

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

    }

  }

  myTTree->Fill();

  return _passed;
}

void lee::ElectronNeutrinoFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fElectronEventSelectionAlg.reconfigure(p.get<fhicl::ParameterSet>("ElectronSelectionAlg"));
  pandoraHelper.reconfigure(p.get<fhicl::ParameterSet>("PandoraInterfaceHelper"));

  m_isOverlaidSample = p.get<bool>("isOverlaidSample", false);
  m_isData = p.get<bool>("isData", false);
}

DEFINE_ART_MODULE(lee::ElectronNeutrinoFilter)
