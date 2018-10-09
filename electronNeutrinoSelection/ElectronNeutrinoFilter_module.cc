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


  myTTree->Branch("true_vx", &_true_vx, "true_vx/d");
  myTTree->Branch("true_vy", &_true_vy, "true_vy/d");
  myTTree->Branch("true_vz", &_true_vz, "true_vz/d");
  myTTree->Branch("ccnc", &_ccnc, "ccnc/I");
  myTTree->Branch("theta", &_theta, "theta/d");
  myTTree->Branch("w", &_w, "w/d");
  myTTree->Branch("qsqr", &_qsqr, "qsqr/d");
  myTTree->Branch("pt", &_pt, "pt/d");
  myTTree->Branch("flash_time", "std::vector< double >", &_flash_time);
  myTTree->Branch("flash_pe", "std::vector< double >", &_flash_pe);
  myTTree->Branch("n_primaries", &_n_primaries, "n_primaries/i");

  myTTree->Branch("passed", &_passed, "passed/O");
  myTTree->Branch("nu_energy", &_nu_energy, "nu_energy/D");
  myTTree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  myTTree->Branch("nu_daughters_E", "std::vector< double >", &_nu_daughters_E);

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

  _true_vx = std::numeric_limits<double>::lowest();
  _true_vy = std::numeric_limits<double>::lowest();
  _true_vz = std::numeric_limits<double>::lowest();
  _true_vx_sce = std::numeric_limits<double>::lowest();
  _true_vy_sce = std::numeric_limits<double>::lowest();
  _true_vz_sce = std::numeric_limits<double>::lowest();

  _flash_time.clear();
  _flash_pe.clear();
  _nu_daughters_p.clear();
  _nu_daughters_start_v.clear();
  _nu_daughters_end_v.clear();

  _nu_daughters_E.clear();
  _nu_daughters_pdg.clear();

  _n_true_nu = std::numeric_limits<unsigned int>::lowest();
  _ccnc = std::numeric_limits<int>::lowest();
  _qsqr = std::numeric_limits<double>::lowest();
  _w = std::numeric_limits<double>::lowest();
  _pt = std::numeric_limits<double>::lowest();
  _theta = std::numeric_limits<double>::lowest();

  _nu_pdg = std::numeric_limits<int>::lowest();
  _interaction_type = std::numeric_limits<int>::lowest();
  _true_nu_is_fiducial = false;
  _nu_energy = std::numeric_limits<double>::lowest();
  _lee_weight = 0;
  _passed = false;
  _n_primaries = 0;
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

  _n_primaries = fElectronEventSelectionAlg.get_primary_indexes().size();
  std::cout << "[ElectronNeutrinoFilter] N primaries " << _n_primaries << std::endl;


  _flash_time = fElectronEventSelectionAlg.get_flash_time();
  std::cout << "[ElectronNeutrinoFilter] Flashes " << _flash_time.size() << std::endl;

  _flash_pe = fElectronEventSelectionAlg.get_flash_PE();

  for (size_t i_fl = 0; i_fl < _flash_time.size(); i_fl++)
    std::cout << "[ElectronNeutrinoFilter] Flash time " << _flash_time[i_fl] << std::endl;

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
  auto const *theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const *detClocks = lar::providerFrom<detinfo::DetectorClocksService>();

  for (auto &gen : generator)
  {
    if (gen.Origin() == simb::kBeamNeutrino)
    {
      _nu_pdg = gen.GetNeutrino().Nu().PdgCode();
      _nu_energy = gen.GetNeutrino().Nu().E();
      _pt = gen.GetNeutrino().Pt();
      _qsqr = gen.GetNeutrino().QSqr();
      _theta = gen.GetNeutrino().Theta();
      _w = gen.GetNeutrino().W();

      if (abs(_nu_pdg) == 12)
      {
        int n_bin = _h_lee_scaling->FindBin(_nu_energy * 1000);
        _lee_weight = _lee_scaling[n_bin];
      }

      _ccnc = gen.GetNeutrino().CCNC();

      _true_vx = gen.GetNeutrino().Nu().Vx();
      _true_vy = gen.GetNeutrino().Nu().Vy();
      _true_vz = gen.GetNeutrino().Nu().Vz();

      _interaction_type = gen.GetNeutrino().Mode();


      double g4Ticks = detClocks->TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + theDetector->GetXTicksOffset(0, 0, 0) - theDetector->TriggerOffset();
      _true_vx_sce =
          _true_vx - sce_service->GetPosOffsets(geo::Point_t(_true_vx, _true_vy, _true_vz)).X() + theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0);
      _true_vy_sce =
          _true_vy + sce_service->GetPosOffsets(geo::Point_t(_true_vx, _true_vy, _true_vz)).Y();
      _true_vz_sce =
          _true_vz + sce_service->GetPosOffsets(geo::Point_t(_true_vx, _true_vy, _true_vz)).Z();


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
