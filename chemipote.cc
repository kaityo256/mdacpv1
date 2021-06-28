//----------------------------------------------------------------------
// Project for measuring chemical potential
//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "chemipote.h"
#include "makeconf.h"
#include "communicator.h"
#include "mt.h"
//----------------------------------------------------------------------
ChemiPote chemipote;
//----------------------------------------------------------------------
static const double DZ = 0.25;
//----------------------------------------------------------------------
class DensityProfile {
private:
  double * density_z;
  int number_of_bins;
  int number_of_observation;
public:
  DensityProfile(MDUnit *mdu) {
    double *L = mdu->GetSystemSize();
    number_of_bins = static_cast<int>(L[Z] / DZ);
    density_z = new double[number_of_bins];
    for (int i = 0; i < number_of_bins; i++) {
      density_z[i] = 0.0;
    }
    number_of_observation = 0;
  };
  ~DensityProfile(void) {
    delete [] density_z;
  };
  void Observe(MDUnit *mdu);
  void SaveToFile(const char* filename, MDUnit *mdu);
};
//----------------------------------------------------------------------
void
DensityProfile::Observe(MDUnit *mdu) {
  double *L = mdu->GetSystemSize();
  const double dz = L[Z] / static_cast<double>(number_of_bins);
  Variables *vars = mdu->GetVariables();
  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  for (int i = 0; i < pn; i++) {
    int index = static_cast<int>(q[i][Z] / dz);
    if (index < 0) {
      index += number_of_bins;
    } else if (index >= number_of_bins) {
      index -= number_of_bins;
    }
    density_z[index] += 1.0;
  }
  number_of_observation++;
}
//----------------------------------------------------------------------
void
DensityProfile::SaveToFile(const char* filename, MDUnit *mdu) {
  double *result_z = new double[number_of_bins];
  Communicator::AllReduceDoubleBuffer(density_z, number_of_bins, result_z);
  const double *L = mdu->GetSystemSize();
  const double inv = 1.0 / (L[X] * L[Y] * DZ * static_cast<double>(number_of_observation));
  const int rank = mdu->GetRank();
  if (0 == rank) {
    std::ofstream ofs(filename);
    for (int i = 0; i < number_of_bins; i++) {
      ofs << (i + 0.5)*DZ << " " << result_z[i]*inv << std::endl;
    }
  }
  delete [] result_z;
}
//----------------------------------------------------------------------
void
ChemiPote::Run(MDUnit *mdu, Parameter &param) {
  mdu->SetPeriodic(true);
  std::string dumpfile = param.GetStringDef("DumpFile", "dumpfile.dmp");
  if (param.GetBooleanDef("Restart", false)) {
    mdu->LoadDumpFile(dumpfile.c_str());
    mdu->CheckPairList();
    mdu->ShowSystemInfo();
    mout << "#Restart from t = " << mdu->GetSimulationTime() << std::endl;

  } else {
    MakeConfHalfFCC(mdu, param);
    mdu->SetInitialVelocity(param.GetDoubleDef("InitialVelocity", 1.0));
  }

  if (param.GetBooleanDef("Thermalize", false)) {
    Thermalize(mdu, param);
  }
  if (param.GetBooleanDef("Observe", false)) {
    Observe(mdu, param);
  }

  if (param.GetBooleanDef("Restart", false)) {
    mout.AppendToFile(param.GetStringDef("OutputFile", "mout.out"));
  } else {
    mout.SaveToFile(param.GetStringDef("OutputFile", "mout.out"));
  }
  //mdu->SaveDumpFile(dumpfile.c_str());
  //mdu->SaveConfigurationSequential();
}
//----------------------------------------------------------------------
void
ChemiPote::Observe(MDUnit *mdu, Parameter &param) {
  mdu->SetControlTemperature(false);
  mdu->StartMeasurement();
  double pressure_sum = 0.0;
  double gamma_sum = 0.0;
  int observe_count = 0;
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  const int TOTAL_LOOP = param.GetIntegerDef("TotalLoop", 20000);
  mdu->SetControlTemperature(false);
  DensityProfile dprofile(mdu);
  for (int i = 0; i < TOTAL_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mdu->CheckPairList();
      //mdu->SaveAsCdviewSequential();
      const double p = mdu->Pressure();
      const double g = Gamma(mdu);
      dprofile.Observe(mdu);
      pressure_sum += p;
      gamma_sum += g;
      observe_count++;
      mout << mdu->GetSimulationTime();
      mout << " " << mdu->Temperature();
      mout << " " << p;
      mout << " " << g;
      mout << " " << mdu->PotentialEnergy();
      mout << " #Observe" << std::endl;
    }
  }
  pressure_sum /= static_cast<double>(observe_count);
  gamma_sum /= static_cast<double>(observe_count);
  mout << "# Gamma = " << gamma_sum << std::endl;
  mout << "# Pressure = " << pressure_sum << std::endl;
  dprofile.SaveToFile(param.GetStringDef("FileName", "test.dat").c_str(), mdu);
  mdu->StopMeasurement();
  mdu->ShowMUPS(TOTAL_LOOP);
}
//----------------------------------------------------------------------
void
ChemiPote::Thermalize(MDUnit *mdu, Parameter &param) {
  const int THERMALIZE_LOOP = param.GetIntegerDef("ThermalizeLoop", 10000);
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  const double temperature = param.GetDoubleDef("AimedTemperature", 0.9);
  mdu->SetControlTemperature(true);
  mdu->ShowSystemInfo();
  for (int i = 0; i < THERMALIZE_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdu->GetSimulationTime();
      mout << " " << mdu->Temperature();
      mout << " " << mdu->Pressure();
      mout << " " << Gamma(mdu);
      mout << " " << mdu->PotentialEnergy();
      mout << " #Thermalize" << std::endl;
    }
  }
}
//----------------------------------------------------------------------
void
ChemiPote::MakeConfHalf(MDUnit *mdu, Parameter &param) {
  double q[D];
  double *L = mdu->GetSystemSize();
  double density = param.GetDoubleDef("LiquidDensity", 0.7);;
  double s = 1.0 / pow(density, 1.0 / 3.0);
  int sx = static_cast<int>(L[X] / s);
  int sy = static_cast<int>(L[Y] / s);
  int sz = static_cast<int>(L[Z] * 0.5 / s);
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        q[X] = (ix + 0.5) * s;
        q[Y] = (iy + 0.5) * s;
        q[Z] = (iz + 0.5) * s;
        mdu->AddParticle(q);
      }
    }
  }
  density = param.GetDoubleDef("GasDensity", 0.1);;
  s = 1.0 / pow(density, 1.0 / 3.0);
  sx = static_cast<int>(L[X] / s);
  sy = static_cast<int>(L[Y] / s);
  sz = static_cast<int>(L[Z] * 0.5 / s);
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        q[X] = (ix + 0.5) * s;
        q[Y] = (iy + 0.5) * s;
        q[Z] = (iz + 0.5) * s + L[Z] * 0.5;
        mdu->AddParticle(q);
      }
    }
  }
}
//----------------------------------------------------------------------
void
ChemiPote::MakeConfHalfFCC(MDUnit *mdu, Parameter &param) {
  double *L = mdu->GetSystemSize();
  double density = param.GetDoubleDef("LiquidDensity", 0.7);;
  double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  double hs = s * 0.5;
  int sx = static_cast<int>(L[X] / s);
  int sy = static_cast<int>(L[Y] / s);
  int sz = static_cast<int>(L[Z] * 0.5 / s);
  double x[D];
  const double e = 0.0000001;
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = static_cast<double>(ix) * s + e;
        x[Y] = static_cast<double>(iy) * s + e;
        x[Z] = static_cast<double>(iz) * s + e;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + e;
        x[Y] = static_cast<double>(iy) * s + hs + e;
        x[Z] = static_cast<double>(iz) * s + hs + e;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs + e;
        x[Y] = static_cast<double>(iy) * s + e;
        x[Z] = static_cast<double>(iz) * s + hs + e;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs + e;
        x[Y] = static_cast<double>(iy) * s + hs + e;
        x[Z] = static_cast<double>(iz) * s + e;
        mdu->AddParticle(x);
      }
    }
  }
  density = param.GetDoubleDef("GasDensity", 0.7);;
  s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  hs = s * 0.5;
  sx = static_cast<int>(L[X] / s);
  sy = static_cast<int>(L[Y] / s);
  sz = static_cast<int>(L[Z] * 0.5 / s);

  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = static_cast<double>(ix) * s + e;
        x[Y] = static_cast<double>(iy) * s + e;
        x[Z] = static_cast<double>(iz) * s + e + L[Z] * 0.5;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + e;
        x[Y] = static_cast<double>(iy) * s + hs + e;
        x[Z] = static_cast<double>(iz) * s + hs + e + L[Z] * 0.5;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs + e;
        x[Y] = static_cast<double>(iy) * s + e;
        x[Z] = static_cast<double>(iz) * s + hs + e + L[Z] * 0.5;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs + e;
        x[Y] = static_cast<double>(iy) * s + hs + e;
        x[Z] = static_cast<double>(iz) * s + e + L[Z] * 0.5;
        mdu->AddParticle(x);
      }
    }
  }
}
//----------------------------------------------------------------------
double
ChemiPote::Gamma(MDUnit *mdu) {
  Variables *vars = mdu->GetVariables();
  double *L = mdu->GetSystemSize();
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  double (*q)[D] = vars->q;
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);
  const int pn = vars->GetParticleNumber();

  double gamma = 0.0;
  int *key_index = mdu->GetKeyParticles();
  int *partner_index = mdu->GetPartnerParticles();

  const int number_of_pairs = mdu->GetPairNumber();
  for (int k = 0; k < number_of_pairs; k++) {
    int i = key_index[k];
    int j = partner_index[k];
    const double dx = q[j][X] - q[i][X];
    const double dy = q[j][Y] - q[i][Y];
    const double dz = q[j][Z] - q[i][Z];
    const double r2 = (dx * dx + dy * dy + dz * dz) * sigma_inv2;
    if (r2 > CL2)continue;
    const double r6 = r2 * r2 * r2;
    double g = -((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * (dz * dz - 0.5 * dx * dx - 0.5 * dy * dy);
    if (j >= pn) {
      g *= 0.5;
    }
    gamma += g;
  }
  gamma = Communicator::AllReduceDouble(gamma);
  gamma = gamma * 0.5 / L[X] / L[Y];
  return gamma;
}
//----------------------------------------------------------------------
