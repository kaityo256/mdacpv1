//----------------------------------------------------------------------
#include <iostream>
#include <string>
#include <fstream>
#include "benchmark.h"
#include "makeconf.h"
//----------------------------------------------------------------------
Benchmark bench;
//----------------------------------------------------------------------
void
Benchmark::Run(MDUnit *mdu, Parameter &param) {

  std::string dumpfile = param.GetStringDef("DumpFileName", "dumpfile.dmp");

  //MT::SetSeed(mdu->GetRank());
  MT::SetSeed(1);

  if (param.GetBooleanDef("Restart", false)) {
    LoadDumpFileWithDuplication(dumpfile.c_str(), mdu);
    mdu->CheckPairList();
    mdu->ShowSystemInfo();
    mout << "#Restart from t = " << mdu->GetSimulationTime() << std::endl;
  } else {
    MakeConfigurationFCC(mdu, param);
    mdu->ShowSystemInfo();
    const double v0 = param.GetDoubleDef("InitialVelocity", 1.0);
    mdu->SetInitialVelocity(v0);
    Thermalize(mdu, param);
  }
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  const int LOOP = param.GetIntegerDef("TotalLoop", 1000);
  mdu->StartMeasurement();
  for (int i = 0; i < LOOP; i++) {
    mdu->CalculateBenchmark();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdu->GetSimulationTime();
      mout << " " << mdu->Temperature();
      mout << " " << mdu->Pressure();
      mout << " " << mdu->TotalEnergy();
      mout << "# observe" << std::endl;
    }
  }
  mdu->StopMeasurement();
  mdu->ShowMUPS(LOOP);
  if (param.GetBooleanDef("SaveDumpFile", false)) {
    mdu->SaveDumpFile(dumpfile.c_str());
  }
}
//----------------------------------------------------------------------
void
Benchmark::Thermalize(MDUnit *mdu, Parameter &param) {
  const int T_LOOP = param.GetIntegerDef("ThermalizeLoop", 150);
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  for (int i = 0; i < T_LOOP; i++) {
    mdu->CalculateBenchmark();
  }
}
//----------------------------------------------------------------------
void
Benchmark::LoadDumpFileWithDuplication(const char *filename, MDUnit *mdu) {
  double x[D], v[D];
  std::ifstream ifs(filename, std::ios::binary);
  int pn;
  double st;
  double zeta;
  Variables *vars = mdu->GetVariables();
  MDRect myrect = mdu->GetRect();
  ifs.read((char*)&pn, sizeof(pn));
  ifs.read((char*)&st, sizeof(st));
  ifs.read((char*)&zeta, sizeof(zeta));
  vars->SetParticleNumber(0);
  vars->SimulationTime = st;
  vars->Zeta = zeta;
  mout << "#Time = " << st << std::endl;
  mout << "#pn = " << pn << std::endl;
  for (int i = 0; i < pn; i++) {
    for (int d = 0; d < D; d++) {
      ifs.read((char*)&x[d], sizeof(double));
      ifs.read((char*)&v[d], sizeof(double));
      x[d] += myrect.s[d];
    }
    mdu->AddParticle(x, v);
  }
}
//----------------------------------------------------------------------
