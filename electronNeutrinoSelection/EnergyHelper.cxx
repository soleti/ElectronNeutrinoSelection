#ifndef ENERGYHELPER_CXX
#define ENERGYHELPER_CXX

#include "EnergyHelper.h"

namespace lee
{

EnergyHelper::EnergyHelper() {

}

void EnergyHelper::reconfigure(fhicl::ParameterSet const &pset) {
  _recombination_factor = pset.get<double>("RecombinationFactor", 0.62);
  _dQdx_rectangle_length = pset.get<double>("dQdxRectangleLength", 4);
  _dQdx_rectangle_width = pset.get<double>("dQdxRectangleWidth", 1);
  _gain = pset.get<std::vector<double>>("Gains", _mc_gain);
  _betap = pset.get<double>("RecombinationBeta", 0.212);
  _alpha = pset.get<double>("RecombinationAlpha", 0.93);
}


void EnergyHelper::track_dQdx(std::vector<art::Ptr<anab::Calorimetry>> *calos,
                              std::vector<double> &dqdx,
                              std::vector<double> &dedx)
{

  int planenum = -1;
  for (auto c : *calos) {
    if (!c)
      continue;

    planenum = c->PlaneID().Plane;
    std::vector<double> dqdxs = c->dQdx();
    std::vector<double> dedxs = c->dEdx();

    if (dqdxs.size() <= 1) {
      continue;
    }

    size_t n = dqdxs.size() / 2;
    std::nth_element(dqdxs.begin(), dqdxs.begin() + n, dqdxs.end());
    dqdx[planenum] = dqdxs[n];
    if (dedxs.size() == dqdxs.size()) {
      dedx[planenum] = dedxs[n];
    }
  }

}

void EnergyHelper::cluster_residuals(std::vector<art::Ptr<recob::Cluster>> *clusters,
                                     art::FindManyP<recob::Hit> *hits_per_cluster,
                                     double &mean_v,
                                     double &std_v)
{

  std::vector<double> distances;

  for (auto _cl: *clusters) {

    if (_cl->Plane().Plane != 2) continue;

    TVector3 start_cluster(_cl->StartWire() * _wire_spacing, _from_tick_to_ns * _drift * _cl->StartTick(), 0);
    TVector3 end_cluster(_cl->EndWire() * _wire_spacing, _from_tick_to_ns * _drift * _cl->EndTick(), 0);
    TVector3 line(start_cluster-end_cluster);

    std::vector< art::Ptr<recob::Hit> > hits = hits_per_cluster->at(_cl.key());

    for (auto &hit : hits)
    {
      double w = hit->WireID().Wire * _wire_spacing;
      double t = _from_tick_to_ns * _drift * hit->PeakTime();
      TVector3 hit_v(w, t, 0);
      TVector3 num(line.Cross(start_cluster - hit_v));
      double side = (hit_v[0] - start_cluster[0])*(end_cluster[1] - start_cluster[1]) - (hit_v[1] - start_cluster[1])*(end_cluster[0] - start_cluster[0]);
      int sign_side = (side > 0) - (side < 0);
      distances.push_back(sign_side * num.Mag() / line.Mag());
    }

  }

  double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
  double mean = sum / distances.size();

  std::vector<double> diff(distances.size());
  std::transform(distances.begin(), distances.end(), diff.begin(), [mean](double x) { return x - mean; });

  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdv = std::sqrt(sq_sum / distances.size());

  if (!isnan(mean) && !isnan(stdv)) {
    mean_v = mean;
    std_v = stdv;
  }

}

void EnergyHelper::energy_from_hits(std::vector<art::Ptr<recob::Cluster>> *clusters,
                                    art::FindManyP<recob::Hit> *hits_per_cluster,
                                    std::vector<int>    &nHits,
                                    std::vector<double> &pfenergy)
{

  nHits.resize(3);
  pfenergy.resize(3);

  for (auto _cl: *clusters) {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster->at(_cl.key());

    for (auto &hit : hits)
    {
      auto plane_nr = hit->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;

      // https://arxiv.org/pdf/1704.02927.pdf
      pfenergy[plane_nr] += 1.01 * hit->Integral() * _gain[plane_nr] * _work_function / _recombination_factor / 1000; // convert MeV to GeV
      nHits[plane_nr]++;
    }
  }
}

void EnergyHelper::PCA(std::vector<art::Ptr<recob::Cluster>> *clusters,
                       art::FindManyP<recob::Hit> *hits_per_cluster,
                       std::vector<std::vector<double>> &pca_planes)
{

  for (auto _cl: *clusters)
  {
    TPrincipal fPrincipal(2, "D");

    std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster->at(_cl.key());
    for (auto &hit : hits)
    {

      double data[2];
      double w = hit->WireID().Wire * _wire_spacing;
      double t = _from_tick_to_ns * _drift * hit->PeakTime();
      data[0] = w;
      data[1] = t;
      fPrincipal.AddRow(data);
    }

    fPrincipal.MakePrincipals();
    pca_planes[0][_cl->Plane().Plane] = (*fPrincipal.GetEigenValues())[0];
    pca_planes[1][_cl->Plane().Plane] = (*fPrincipal.GetEigenValues())[1];
  }

}

void EnergyHelper::get_cali(
    std::vector<art::Ptr<recob::SpacePoint>> *spcpnts,
    art::FindManyP<recob::Hit> *hits_per_spcpnts,
    std::vector<double> &cali_corr)
{
  cali_corr.resize(3);
  std::vector<float> total_charge(3, 0);

  for (auto _sps : *spcpnts)
  {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts->at(_sps.key());
    const double *xyz = _sps->XYZ();

    for (auto &hit : hits)
    {
      auto plane_nr = hit->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;
      total_charge[plane_nr] += hit->Integral();
      float yzcorrection = _energy_calib_provider.YZdqdxCorrection(plane_nr, xyz[1], xyz[2]);
      float xcorrection = _energy_calib_provider.XdqdxCorrection(plane_nr, xyz[0]);
      if (!yzcorrection)
        yzcorrection = 1.0;
      if (!xcorrection)
        xcorrection = 1.0;
      cali_corr[plane_nr] += yzcorrection * xcorrection * hit->Integral();
    }
  }

  for (unsigned short i = 0; i < 3; ++i)
  {
    if (total_charge[i] > 0)
    {
      cali_corr[i] /= total_charge[i];
    }
    else
    {
      cali_corr[i] = 1;
    }
  }
}

double EnergyHelper::PID(art::Ptr<anab::ParticleID> selected_pid,
                         std::string AlgName,
                         anab::kVariableType VariableType,
                         anab::kTrackDir TrackDirection,
                         int pdgCode)
{

    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = selected_pid->ParticleIDAlgScores();
    for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
    {
      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
      int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);
      if (planeid < 0 || planeid > 2)
      {
        std::cout << "[EnergyHelper::PID] No information for planeid " << planeid << std::endl;
        continue;
      }

      if (AlgScore.fAlgName == AlgName)
      {
        if (anab::kVariableType(AlgScore.fVariableType) == VariableType
            && anab::kTrackDir(AlgScore.fTrackDir) == TrackDirection)
        {
          if (AlgScore.fAssumedPdg == pdgCode) {
              double alg_value = AlgScore.fValue;
              return alg_value;
          }
        }
      }
    }
    return std::numeric_limits<double>::lowest();
}

void EnergyHelper::dQdx_cali(const recob::Shower *shower_obj,
                             std::vector<double> &dqdx_cali)
{
  TVector3 pfp_dir;

  // Field needed for calibration factor
  double x_start, y_start, z_start;
  double x_middle, y_middle, z_middle;
  double x_end, y_end, z_end;
  float start_corr, middle_corr, end_corr;

  pfp_dir.SetX(shower_obj->Direction().X());
  pfp_dir.SetY(shower_obj->Direction().Y());
  pfp_dir.SetZ(shower_obj->Direction().Z());

  x_start = shower_obj->ShowerStart().X();
  y_start = shower_obj->ShowerStart().Y();
  z_start = shower_obj->ShowerStart().Z();

  pfp_dir.SetMag(_dQdx_rectangle_length / 2.); //Go 2cm along the direction of the object.
  x_middle = x_start + pfp_dir.X();
  y_middle = y_start + pfp_dir.Y();
  z_middle = z_start + pfp_dir.Z();
  x_end = x_middle + pfp_dir.X();
  y_end = y_middle + pfp_dir.Y();
  z_end = z_middle + pfp_dir.Z();
  pfp_dir.SetMag(1.); //Normalise again for safety (not needed).

  for (int plane_nr = 0; plane_nr < 3; ++plane_nr)
  {
    float yzcorrection_start = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_start, z_start);
    float xcorrection_start = _energy_calib_provider.XdqdxCorrection(plane_nr, x_start);
    if (!yzcorrection_start)
      yzcorrection_start = 1.0;
    if (!xcorrection_start)
      xcorrection_start = 1.0;
    start_corr = yzcorrection_start * xcorrection_start;

    float yzcorrection_middle = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_middle, z_middle);
    float xcorrection_middle = _energy_calib_provider.XdqdxCorrection(plane_nr, x_middle);
    if (!yzcorrection_middle)
      yzcorrection_middle = 1.0;
    if (!xcorrection_middle)
      xcorrection_middle = 1.0;
    middle_corr = yzcorrection_middle * xcorrection_middle;

    float yzcorrection_end = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_end, z_end);
    float xcorrection_end = _energy_calib_provider.XdqdxCorrection(plane_nr, x_end);
    if (!yzcorrection_end)
      yzcorrection_end = 1.0;
    if (!xcorrection_end)
      xcorrection_end = 1.0;
    end_corr = yzcorrection_end * xcorrection_end;
    //std::cout << "[EnergyHelper] dqdx_cali " << start_corr << middle_corr << end_corr << std::endl;
    dqdx_cali[plane_nr] = (start_corr + middle_corr + end_corr) / 3;
  }
}

void EnergyHelper::dQdx(const recob::Shower *shower_obj,
                        std::vector<art::Ptr<recob::Cluster>> *clusters,
                        art::FindManyP<recob::Hit> *hits_per_cluster,
                        std::vector<double> &dqdx,
                        std::vector<std::vector<double>> &dqdx_hits,
                        std::vector<double> &pitches,
                        std::vector<int> &dqdx_hits_in_the_box)
{

  double tolerance = 0.001;

  TVector3 pfp_dir;

  pfp_dir.SetX(shower_obj->Direction().X());
  pfp_dir.SetY(shower_obj->Direction().Y());
  pfp_dir.SetZ(shower_obj->Direction().Z());

  for (int i=0; i<3; i++)
  {
    pitches[i] = geo_helper.getPitch(pfp_dir, i);
  }

  for (auto _cl: *clusters)
  {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster->at(_cl.key());

    std::vector<double> cluster_axis;
    std::vector<double> cluster_start;
    std::vector<double> cluster_end;
    double start_x = _detprop->ConvertTicksToX(_cl->StartTick(), _cl->Plane());
    double end_x = _detprop->ConvertTicksToX(_cl->EndTick(), _cl->Plane());
    // std::cout << "start " << start_x << ", end: " << end_x << std::endl;
    double pitch = pitches[_cl->Plane().Plane];

    if (pitch >= 0)
    {
      std::reverse(hits.begin(), hits.end());
      cluster_axis = {cos(_cl->StartAngle()),
                      sin(_cl->StartAngle())};

      cluster_start = {_cl->StartWire() * _wire_spacing - tolerance * cos(_cl->StartAngle()),
                       start_x - tolerance * sin(_cl->StartAngle())};
      cluster_end = {_cl->EndWire() * _wire_spacing, end_x};
    }
    else
    {
      cluster_axis = {-1. * cos(_cl->StartAngle()),
                      -1. * sin(_cl->StartAngle())};
      cluster_start = {_cl->EndWire() * _wire_spacing + tolerance * cos(_cl->StartAngle()),
                       end_x + tolerance * sin(_cl->StartAngle())};
      cluster_end = {_cl->StartWire() * _wire_spacing, start_x};
    }

    double cluster_length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) +
                                 pow(cluster_end[1] - cluster_start[1], 2));
    if (cluster_length <= 0)
      continue;


    // Build rectangle 4 x 1 cm around the cluster axis
    std::vector<std::vector<double>> points;
    geo_helper.buildRectangle(_dQdx_rectangle_length, _dQdx_rectangle_width,
                             cluster_start, cluster_axis, points);

    std::vector<double> dqdxs;

    for (auto &hit : hits)
    {
      std::vector<double> hit_pos = {hit->WireID().Wire * _wire_spacing,
                                     _detprop->ConvertTicksToX(hit->PeakTime(),
                                     _cl->Plane())};

      bool is_within = geo_helper.isInside(hit_pos, points);

      if (is_within)
      {
        double q = hit->Integral() * _gain[_cl->Plane().Plane];
        dqdxs.push_back(q / fabs(pitch));
        dqdx_hits[_cl->Plane().Plane].push_back(q / fabs(pitch));
      }
    }

    dqdx_hits_in_the_box[_cl->Plane().Plane] = dqdxs.size();

    // Get the median
    if (dqdxs.size() > 0)
    {
      std::nth_element(dqdxs.begin(), dqdxs.begin() + dqdxs.size() / 2, dqdxs.end());
      dqdx[_cl->Plane().Plane] = dqdxs[dqdxs.size() / 2];
      pitches[_cl->Plane().Plane] = pitch;
    }
  }
}

void EnergyHelper::dEdx_from_dQdx(std::vector<double> &dedx,
                                  std::vector<double> dqdx)
{
  for (size_t i = 0; i < dqdx.size(); i++)
  {
    dedx.push_back(charge2energy(dqdx[i]));
  }
}

double EnergyHelper::charge2energy(double charge)
{
  double Rho = 1.4;
  double Efield = 0.273;

  return (exp(charge*(_betap/(Rho*Efield))*_work_function)-_alpha)/(_betap/(Rho*Efield));
}

void EnergyHelper::energy_from_hits_new_method(std::vector<art::Ptr<recob::Cluster>> *clusters,
                                    art::FindManyP<recob::Hit> *hits_per_cluster,
                                    std::vector<int>    &nHits,
                                    std::vector<double> &pfenergy)
{
  nHits.resize(3);
  pfenergy.resize(3);

  for (auto _cl: *clusters) {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster->at(_cl.key());

    for (auto &hit : hits)
    {
      auto plane_nr = hit->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;

      pfenergy[plane_nr] += charge2energy(hit->Integral());
      nHits[plane_nr]++;
    }
  }
}

} // namespace lee

#endif
