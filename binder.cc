//----------------------------------------------------------------------
// Project for measuring Binder parameter
//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "binder.h"
#include "makeconf.h"
#include "communicator.h"
#include "mt.h"
//----------------------------------------------------------------------
Binder binder;
//----------------------------------------------------------------------
class BinderCalculatorMPI {
private:
  int msize;
  int msize_local;
  double dsize;
  double *dprof;
  double *dprof_local;
  int observe_count;
  double value;
public:
  BinderCalculatorMPI(const int m, MDUnit *mdu) {
    msize = m;
    double *L = mdu->GetSystemSize();
    dsize = L[X] / static_cast<double>(msize);
    dprof = new double[msize * msize * msize];
    observe_count = 0;
    MDRect r = mdu->GetRect();
    msize_local = static_cast<int>(r.GetWidth(X) / dsize);
    dprof_local = new double[msize_local * msize_local * msize_local];
  };
  ~BinderCalculatorMPI(void) {
    delete [] dprof;
    delete [] dprof_local;
  };
  double Calculate(MDUnit *mdu);
  double GetValue(void) {return value / static_cast<double>(observe_count);};
  double GetSize(void) {return dsize;};
};
//----------------------------------------------------------------------
double
BinderCalculatorMPI::Calculate(MDUnit *mdu) {
  observe_count++;
  mdu->CheckPairList();
  Variables *vars = mdu->GetVariables();
  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  MDRect myrect = mdu->GetRect();
  double *L = mdu->GetSystemSize();
  const int number_of_localgrids = msize_local * msize_local * msize_local;

  for (int i = 0; i < number_of_localgrids; i++) {
    dprof_local[i] = 0.0;
  }

  const double d_inv = 1.0 / dsize;
  const double v_inv = 1.0 / (dsize * dsize * dsize);
  for (int i = 0; i < pn; i++) {
    int ix = static_cast<int>((q[i][X] - myrect.s[X]) * d_inv);
    int iy = static_cast<int>((q[i][Y] - myrect.s[Y]) * d_inv);
    int iz = static_cast<int>((q[i][Z] - myrect.s[Z]) * d_inv);
    int index = ix + iy * msize_local + iz * msize_local * msize_local;
    if ( q[i][X] < myrect.s[X] || q[i][X] > myrect.e[X] ||
         q[i][Y] < myrect.s[Y] || q[i][Y] > myrect.e[Y] ||
         q[i][Z] < myrect.s[Z] || q[i][Z] > myrect.e[Z] ||
         index < 0 || index >= number_of_localgrids) {
      show_error("Invalid Index.");
      mout << mdu->GetSimulationTime() << std::endl;
      mout << ix << " " << iy << " " << iz << std::endl;
      mout << q[i][X] << " " << q[i][Y] << " " << q[i][Z] << std::endl;
      exit(1);
    }
    dprof_local[index] = dprof_local[index] + v_inv;
  }

  Communicator::AllGatherDouble(dprof_local, number_of_localgrids, dprof);
  const int number_of_grids = msize * msize * msize;

  double m0 = 0.0;
  for (int i = 0; i < number_of_grids; i++) {
    m0 += dprof[i];
  }
  m0 = m0 / (number_of_grids);

  double m2 = 0.0;
  double m4 = 0.0;

  for (int i = 0; i < number_of_grids; i++) {
    double phi = dprof[i] - m0;
    m2 += phi * phi;
    m4 += phi * phi * phi * phi;
  }
  m2 /= (number_of_grids);
  m4 /= (number_of_grids);
  value = value + m4 / (m2 * m2);
  return m4 / (m2 * m2);
}
//----------------------------------------------------------------------
class BinderCalculator {
private:
  int msize;
  double dsize;
  double *dprof;
  int observe_count;
  double value;
public:
  BinderCalculator(const int m, MDUnit *mdu) {
    msize = m;
    double *L = mdu->GetSystemSize();
    dsize = L[X] / static_cast<double>(msize);
    dprof = new double[msize * msize * msize];
    observe_count = 0;
  };
  ~BinderCalculator(void) {
    delete [] dprof;
  };
  double Calculate(MDUnit *mdu);
  double GetValue(void) {return value / static_cast<double>(observe_count);};
  double GetSize(void) {return dsize;};
};
//----------------------------------------------------------------------
double
BinderCalculator::Calculate(MDUnit *mdu) {

  observe_count++;
  mdu->CheckPairList();
  Variables *vars = mdu->GetVariables();
  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double *L = mdu->GetSystemSize();
  const int number_of_grids = msize * msize * msize;
  for (int i = 0; i < number_of_grids; i++) {
    dprof[i] = 0.0;
  }
  const double d_inv = 1.0 / dsize;
  const double v_inv = 1.0 / (dsize * dsize * dsize);
  for (int i = 0; i < pn; i++) {
    int ix = static_cast<int>(q[i][X] * d_inv);
    int iy = static_cast<int>(q[i][Y] * d_inv);
    int iz = static_cast<int>(q[i][Z] * d_inv);
    int index = ix + iy * msize + iz * msize * msize;
    if ( q[i][X] < 0 || q[i][X] > L[X] ||
         q[i][Y] < 0 || q[i][Y] > L[Y] ||
         q[i][Z] < 0 || q[i][Z] > L[Z]) {
      show_error("Invalid Index.");
      mout << mdu->GetSimulationTime() << std::endl;
      mout << ix << " " << iy << " " << iz << std::endl;
      mout << q[i][X] << " " << q[i][Y] << " " << q[i][Z] << std::endl;
      exit(1);
    }
    dprof[index] = dprof[index] + v_inv;
  }

  double m0 = 0.0;
  for (int i = 0; i < number_of_grids; i++) {
    m0 += dprof[i];
  }
  m0 = m0 / (number_of_grids);

  double m2 = 0.0;
  double m4 = 0.0;

  for (int i = 0; i < number_of_grids; i++) {
    double phi = dprof[i] - m0;
    m2 += phi * phi;
    m4 += phi * phi * phi * phi;
  }
  m2 /= (number_of_grids);
  m4 /= (number_of_grids);
  value = value + m4 / (m2 * m2);
  return m4 / (m2 * m2);
}
//----------------------------------------------------------------------
void
Binder::Run(MDUnit *mdu, Parameter &param) {
  mdu->SetPeriodic(true);
  std::string dumpfile = param.GetStringDef("DumpFile", "dumpfile.dmp");
  if (param.GetBooleanDef("Restart", false)) {
    mdu->LoadDumpFile(dumpfile.c_str());
    mdu->CheckPairList();
    mdu->ShowSystemInfo();
    mout << "#Restart from t = " << mdu->GetSimulationTime() << std::endl;
  } else {
    MakeConfigurationPeriodic(mdu, param);
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
  if (param.GetBooleanDef("SaveDumpFile", false)) {
    mdu->SaveDumpFile(dumpfile.c_str());
  }
}
//----------------------------------------------------------------------
void
Binder::Observe(MDUnit *mdu, Parameter &param) {
  mdu->SetControlTemperature(false);
  mdu->StartMeasurement();
  double pressure_sum = 0.0;
  int observe_count = 0;
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  const int TOTAL_LOOP = param.GetIntegerDef("TotalLoop", 20000);
  const double temperature = param.GetDoubleDef("AimedTemperature", 0.9);
  mdu->SetControlTemperature(false);
  BinderCalculatorMPI b4(4, mdu);
  for (int i = 0; i < TOTAL_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mdu->CheckPairList();
      const double p = mdu->Pressure();
      pressure_sum += p;
      observe_count++;
      mout << mdu->GetSimulationTime();
      mout << " " << mdu->Temperature();
      mout << " " << p;
      mout << " " << mdu->PotentialEnergy();
      mout << " " << b4.Calculate(mdu);
      mout << " #Observe" << std::endl;
    }
  }
  pressure_sum /= static_cast<double>(observe_count);
  mout << "# Pressure = " << pressure_sum << std::endl;
  mout << "# " << temperature << " " << b4.GetValue() << " # Binder" << std::endl;
  mdu->StopMeasurement();
  mdu->ShowMUPS(TOTAL_LOOP);
}
//----------------------------------------------------------------------
void
Binder::Thermalize(MDUnit *mdu, Parameter &param) {
  const int THERMALIZE_LOOP = param.GetIntegerDef("ThermalizeLoop", 10000);
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  const double temperature = param.GetDoubleDef("AimedTemperature", 0.9);
  mdu->SetControlTemperature(true);
  mdu->ShowSystemInfo();
  BinderCalculatorMPI b4(4, mdu);
  for (int i = 0; i < THERMALIZE_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdu->GetSimulationTime();
      mout << " " << mdu->Temperature();
      mout << " " << mdu->Pressure();
      mout << " " << mdu->PotentialEnergy();
      mout << " " << b4.Calculate(mdu);
      mout << " #Thermalize" << std::endl;
    }
  }
}
//----------------------------------------------------------------------
