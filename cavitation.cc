//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "cavitation.h"
#include "bubblehist.h"
#include "communicator.h"
#include "mt.h"
//----------------------------------------------------------------------
Cavitation cavitation;
//----------------------------------------------------------------------
int
GetSeed(int rank) {
  std::ifstream ifs("seed.dat");
  int seed = 0;
  ifs >> seed;
  ifs.close();
  Communicator::Barrier();
  if (0 == rank) {
    std::ofstream ofs("seed.dat");
    ofs << seed + 1 << std::endl;
  }
  return seed;
}
//----------------------------------------------------------------------
void
Cavitation::Run(MDUnit *mdu, Parameter &param) {
  const double aimed_density = param.GetDoubleDef("AimedDensity", 0.658);
  const double density = param.GetDoubleDef("Density", 0.70);
  const double sigma = pow(aimed_density / density, 1.0 / 3.0);
  const double temperature = param.GetDoubleDef("AimedTemperature", 0.9);
  double *L = mdu->GetSystemSize();
  char filename[256];
  char output_filename[256];
  int seed = 1;
  if (param.GetBooleanDef("ReadSeed", true)) {
    seed = GetSeed(mdu->GetRank());
  }
  MT::SetSeed(seed);
  sprintf(filename, "S%03d.hist", seed);
  sprintf(output_filename, "S%03d.dat", seed);
  mdu->SetPeriodic(true);
  MakeConfigurationExtract(mdu, param);
  mdu->MakePairList();
  mdu->SetInitialVelocity(param.GetDoubleDef("InitialVelocity", 1.0));
  mdu->SetControlTemperature(true);
  mdu->ShowSystemInfo();
  mout << "#time,max_size,max_size2,ratio,temperature,pressure" << std::endl;

  Thermalize(mdu, param);
  if (param.GetBooleanDef("ThermalizeOnly", false)) {
    return;
  }

  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  mdu->SetControlTemperature(false);
  mdu->ChangeScale(sigma);
  mdu->ShowSystemInfo();
  mdu->StartMeasurement();
  BubbleHist bhist(mdu->GetRank(), mdu->GetProcessNumber(), mdu->GetRect(), mdu->GetSystemSize());

  const int RELAXATION_LOOP = param.GetIntegerDef("RelaxationLoop", 100);

  double start_time = mdu->GetSimulationTime();

  const double save_grid_size = param.GetDoubleDef("SaveGridSize", 2.0);

  for (int i = 0; i < RELAXATION_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mdu->MakePairList();
//      mdu->SaveDensityBinary(save_grid_size);
      bhist.Analyse(mdu->GetVariables());
      double max_size = bhist.GetClusterSize(0);
      double max_size2 = bhist.GetClusterSize(1);
      double ratio = max_size2 / max_size;
      mout << mdu->GetSimulationTime();
      mout << " " << max_size;
      mout << " " << max_size2;
      mout << " " << ratio;
      mout << " " << mdu->Temperature();
      mout << " " << mdu->Pressure();
      mout << " #relax" << std::endl;
    }
  }
  bhist.Clear();

  const int TOTAL_LOOP = param.GetIntegerDef("TotalLoop", 1000);
  const double MAX_CLUSTER_SIZE = param.GetDoubleDef("MaxClusterSize", 500.0);
  std::string base_dir = param.GetStringDef("BaseDir", ".");

  for (int i = 0; i < TOTAL_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mdu->MakePairList();
      //mdu->SaveDensityBinary(save_grid_size);
      //mdu->SaveDumpFileSequential();
      //bhist.SaveBubbleDistribution(mdu,base_dir.c_str());
      bhist.Analyse(mdu->GetVariables());
      double max_size = bhist.GetClusterSize(0);
      double max_size2 = bhist.GetClusterSize(1);
      double ratio = max_size2 / max_size;
      mout << mdu->GetSimulationTime();
      mout << " " << max_size;
      mout << " " << max_size2;
      mout << " " << ratio;
      mout << " " << mdu->Temperature();
      mout << " " << mdu->Pressure();
      mout << " " << mdu->TotalEnergy();
      mout << " #observe" << std::endl;
      if (max_size > MAX_CLUSTER_SIZE) {
        double w_time = mdu->GetSimulationTime() - start_time;
        mout << "# Waiting time = " << w_time << std::endl;
        break;
      }
    }
  }
  if (param.GetBooleanDef("ReadSeed", true)) {
    mout.SaveToFile(output_filename);
  }
  mdu->StopMeasurement();
  mdu->ShowMUPS(TOTAL_LOOP);
}
//----------------------------------------------------------------------
void
Cavitation::Thermalize(MDUnit *mdu, Parameter &param) {
  const int T_LOOP = param.GetIntegerDef("ThermalizeLoop", 10000);
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  BubbleHist bhist(mdu->GetRank(), mdu->GetProcessNumber(), mdu->GetRect(), mdu->GetSystemSize());
  for (int i = 0; i < T_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mdu->MakePairList();
      bhist.Analyse(mdu->GetVariables());
      double max_size = bhist.GetClusterSize(0);
      double max_size2 = bhist.GetClusterSize(1);
      double ratio = max_size2 / max_size;
      mout << mdu->GetSimulationTime();
      mout << " " << max_size;
      mout << " " << max_size2;
      mout << " " << ratio;
      mout << " " << mdu->Temperature();
      mout << " " << mdu->Pressure();
      mout << " #thermalize" << std::endl;
    }
  }
}
//----------------------------------------------------------------------
void
Cavitation::MakeConfiguration(MDUnit *mdu, Parameter &param) {
  double *L = mdu->GetSystemSize();
  const double density = param.GetDoubleDef("Density", 0.7);;
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L[X] / s);
  int sy = static_cast<int>(L[Y] / s);
  int sz = static_cast<int>(L[Z] / s);
  double x[D];
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s + hs;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s + hs;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s;
        mdu->AddParticle(x);
      }
    }
  }
}
//----------------------------------------------------------------------
void
Cavitation::MakeConfigurationExtract(MDUnit *mdu, Parameter &param) {
  double *L = mdu->GetSystemSize();
  const double density = param.GetDoubleDef("Density", 0.7);;
  double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L[X] / s) + 1;
  int sy = static_cast<int>(L[Y] / s) + 1;
  int sz = static_cast<int>(L[Z] / s) + 1;
  s = L[X] / static_cast<double>(sx);
  double x[D];
  int n_full = sx * sy * sz * 4;
  int n_exact = static_cast<int>(L[X] * L[Y] * L[Z] * density);
  int n_delete = n_full - n_exact;
  std::vector <int> vec_delete;
  MT::SetSeed(1);
  while (vec_delete.size() < n_delete) {
    int v = static_cast<int>(n_full * MT::GetDouble());
    if (std::find(vec_delete.begin(), vec_delete.end(), v) == vec_delete.end()) {
      vec_delete.push_back(v);
    }
  }
  std::sort(vec_delete.begin(), vec_delete.end());
  int n = 0;
  int n_pos = 0;
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s + hs;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s + hs;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;
      }
    }
  }
}
//----------------------------------------------------------------------
void
Cavitation::MakeBubble(MDUnit *mdu, Parameter &param) {
  double *L = mdu->GetSystemSize();
  const double density = param.GetDoubleDef("Density", 0.7);;
  double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  const double cx = 7.0;
  const double cy = 7.0;
  const double cz = 7.0;

  const double cx2 = 47.0;
  const double cy2 = 47.0;
  const double cz2 = 47.0;

  const double R = 6.0;
  const double R2 = 5.0;
  int sx = static_cast<int>(L[X] / s);
  int sy = static_cast<int>(L[Y] / s);
  int sz = static_cast<int>(L[Z] / s);
  s = L[X] / sx;
  double x[D];
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s;
        if ((x[X] - cx) * (x[X] - cx) + (x[Y] - cy) * (x[Y] - cy) + (x[Z] - cz) * (x[Z] - cz) > R * R &&
            (x[X] - cx2) * (x[X] - cx2) + (x[Y] - cy2) * (x[Y] - cy2) + (x[Z] - cz2) * (x[Z] - cz2) > R2 * R2)
          mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s + hs;
        if ((x[X] - cx) * (x[X] - cx) + (x[Y] - cy) * (x[Y] - cy) + (x[Z] - cz) * (x[Z] - cz) > R * R &&
            (x[X] - cx2) * (x[X] - cx2) + (x[Y] - cy2) * (x[Y] - cy2) + (x[Z] - cz2) * (x[Z] - cz2) > R2 * R2)
          mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s + hs;
        if ((x[X] - cx) * (x[X] - cx) + (x[Y] - cy) * (x[Y] - cy) + (x[Z] - cz) * (x[Z] - cz) > R * R &&
            (x[X] - cx2) * (x[X] - cx2) + (x[Y] - cy2) * (x[Y] - cy2) + (x[Z] - cz2) * (x[Z] - cz2) > R2 * R2)
          mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s;
        if ((x[X] - cx) * (x[X] - cx) + (x[Y] - cy) * (x[Y] - cy) + (x[Z] - cz) * (x[Z] - cz) > R * R &&
            (x[X] - cx2) * (x[X] - cx2) + (x[Y] - cy2) * (x[Y] - cy2) + (x[Z] - cz2) * (x[Z] - cz2) > R2 * R2)
          mdu->AddParticle(x);
      }
    }
  }
}
//----------------------------------------------------------------------
