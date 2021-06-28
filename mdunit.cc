//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <vector>
#include "mdunit.h"
#include "fcalculator.h"
#include "communicator.h"
#include "measure.h"
#include "stopwatch.h"
//----------------------------------------------------------------------
MDUnit::MDUnit(MPIInfo *mi, Parameter &param) {
  minfo = mi;
  vars = new Variables();
  obs = new Observer();
  plist = new PairList();
  measure = new Measure();
  cmanager = new CommunicationManager();
  minfo->GetGridPosition(grid_position);
  int grid_size[D];
  minfo->GetGridSize(grid_size);
  sinfo = new SimulationInfo(param, grid_size);
  double s[D], e[D];
  for (int d = 0; d < D; d++) {
    double ul = sinfo->L[d] / static_cast<double>(grid_size[d]);
    s[d] = ul * static_cast<double>(grid_position[d]);
    e[d] = s[d] + ul;
  }
  myrect = MDRect(s, e);
  mesh = new MeshList(sinfo, vars, myrect);
}
//----------------------------------------------------------------------
MDUnit::~MDUnit(void) {
  delete measure;
  delete plist;
  delete obs;
  delete mesh;
  delete vars;
  delete cmanager;
  delete sinfo;
}
//----------------------------------------------------------------------
void
MDUnit::ShowSystemInfo(void) {
  const unsigned long int tn = obs->GetTotalParticleNumber(vars);
  const double sigma = vars->GetSigma();
  sinfo->ShowAll(tn, sigma);
}
//----------------------------------------------------------------------
void
MDUnit::CheckPairList(void) {
  if (!plist->IsFresh()) {
    MakePairList();
  }
}
//----------------------------------------------------------------------
void
MDUnit::MakePairList(void) {
  SendParticles();
  vars->AdjustPeriodicBoundary(sinfo);
  mesh->Sort(vars, sinfo, myrect);
  MakeNeighborInformation();
  plist->Init(vars, sinfo);
  mesh->MakeList(vars, sinfo, myrect);
}
//----------------------------------------------------------------------
void
MDUnit::ChangeScale(double s) {
  vars->SetSigma(s);
  mesh->ChangeScale(s, sinfo, myrect);
  MakePairList();
  mout << "# Change Scale to " << s << std::endl;
}
//----------------------------------------------------------------------
void
MDUnit::Calculate(void) {
  if (IsPairListExpired()) {
    MakePairList();
  }

  ForceCalculator::UpdatePositionHalf(vars, sinfo);
  SendNeighborInformation();

  if (sinfo->ControlTemperature) {
    CalculateHeatbath();
  } else {
    ForceCalculator::CalculateForce(vars, mesh, sinfo, minfo);
  }

  ForceCalculator::UpdatePositionHalf(vars, sinfo);

  vars->SimulationTime += sinfo->TimeStep;
  plist->SetFresh(false);
}
//----------------------------------------------------------------------
void
MDUnit::CalculateHeatbath(void) {
  if (sinfo->HeatbathType == HT_NOSEHOOVER) {
    ForceCalculator::HeatbathZeta(vars->Zeta, Temperature(), sinfo);
    ForceCalculator::HeatbathMomenta(vars, sinfo);
    ForceCalculator::CalculateForce(vars, mesh, sinfo, minfo);
    ForceCalculator::HeatbathMomenta(vars, sinfo);
    ForceCalculator::HeatbathZeta(vars->Zeta, Temperature(), sinfo);
  } else if (sinfo->HeatbathType == HT_LANGEVIN) {
    ForceCalculator::CalculateForce(vars, mesh, sinfo, minfo);
    ForceCalculator::Langevin(vars, sinfo);
  }
}
//----------------------------------------------------------------------
void
MDUnit::CalculateBenchmark(void) {
  static StopWatch swAll(minfo->GetRank(), "all");
  static StopWatch swForce(minfo->GetRank(), "force");
  static StopWatch swComm(minfo->GetRank(), "comm");
  static StopWatch swPair(minfo->GetRank(), "pair");
  swAll.Start(GetSimulationTime());
  if (IsPairListExpired()) {
    swPair.Start(GetSimulationTime());
    MakePairList();
    //mout << "# " << GetSimulationTime() << " #Expired!" << std::endl;
    swPair.Stop();
  }
  ForceCalculator::UpdatePositionHalf(vars, sinfo);

  swComm.Start(GetSimulationTime());
  SendNeighborInformation();
  swComm.Stop();

  swForce.Start(GetSimulationTime());
  ForceCalculator::CalculateForceBenchmark(vars, mesh, sinfo);
  swForce.Stop();
  ForceCalculator::UpdatePositionHalf(vars, sinfo);
  vars->SimulationTime += sinfo->TimeStep;
  plist->SetFresh(false);
  swAll.Stop();
}
//----------------------------------------------------------------------
bool
MDUnit::IsPairListExpired(void) {
  bool expired = plist->IsPairListExpired(vars, mesh, sinfo);
  return Communicator::AllReduceBoolean(expired);
}
//----------------------------------------------------------------------
void
MDUnit::AddParticle(double x[D]) {
  double v[D];
  v[X] = 0.0;
  v[Y] = 0.0;
  v[Z] = 0.0;
  AddParticle(x, v);
}
//----------------------------------------------------------------------
void
MDUnit::AddParticle(double x[D], double v[D]) {
  if (myrect.IsInside(x)) {
    vars->AddParticle(x, v);
  }
}
//----------------------------------------------------------------------
void
MDUnit::SaveConfigurationSequential(void) {
  static int index = 0;
  char filename[256];
  sprintf(filename, "%s/conf%03d.cnf", sinfo->BaseDir.c_str(), index);
  index++;
  SaveConfiguration(filename);
}
//----------------------------------------------------------------------
void
MDUnit::SaveConfiguration(const char *filename) {
  CheckPairList();
  const unsigned long int tn = obs->GetTotalParticleNumber(vars);
  if (0 == minfo->GetRank()) {
    std::ofstream fs(filename);
    if (!fs) {
      std::cout << "Could not open " << filename << std::endl;
    }
    fs << "# N=" << tn;
    fs << ",LX=" << sinfo->L[X];
    fs << ",LY=" << sinfo->L[Y];
    fs << ",LZ=" << sinfo->L[Z];
    fs << ",Time=" << vars->SimulationTime;
    fs << std::endl;
    fs.close();
  }
  for (int r = 0; r < minfo->GetProcessNumber(); r++) {
    Communicator::Barrier();
    if (minfo->GetRank() != r)continue;
    std::ofstream fs(filename, std::ios::app);
    vars->SaveConfiguration(fs);
    fs.close();
  }
}
//----------------------------------------------------------------------
void
MDUnit::SaveAsCdviewSequential(void) {
  static int index = 0;
  char filename[256];
  sprintf(filename, "%s/conf%03d.cd", sinfo->BaseDir.c_str(), index);
  index++;
  SaveAsCdview(filename);
}
//----------------------------------------------------------------------
void
MDUnit::SaveAsCdview(const char *filename) {
  CheckPairList();
  const unsigned long int tn = obs->GetTotalParticleNumber(vars);
  if (0 == minfo->GetRank()) {
    std::ofstream fs(filename);
    if (!fs) {
      std::cout << "Could not open " << filename << std::endl;
    }
    fs << "# N=" << tn;
    fs << ",LX=" << sinfo->L[X];
    fs << ",LY=" << sinfo->L[Y];
    fs << ",LZ=" << sinfo->L[Z];
    fs << ",Time=" << vars->SimulationTime;
    fs << std::endl;
    fs.close();
  }
  int num = 0;
  for (int r = 0; r < minfo->GetProcessNumber(); r++) {
    Communicator::Barrier();
    double (*q)[D] = vars->q;
    double (*p)[D] = vars->p;
    if (minfo->GetRank() == r) {
      std::ofstream fs(filename, std::ios::app);
      for (int i = 0; i < vars->GetParticleNumber(); i++) {
        fs << num << " 0 " << q[i][X] << " " << q[i][Y] << " " << q[i][Z];
        fs << " " << p[i][X] << " " << p[i][Y] << " " << p[i][Z] << "\n";
        num++;
      }
      fs.close();
    }
    Communicator::BroadcastInteger(num, r);
  }
}
//----------------------------------------------------------------------
void
MDUnit::LoadCdview(const char *filename) {
  std::ifstream ifs(filename);
  if (!ifs) {
    mout << "# Could not open " << filename << std::endl;
    return;
  }
  std::string line;
  while (!ifs.eof()) {
    if (ifs.peek() == '#') {
      getline(ifs, line);
      continue;
    }
    getline(ifs, line);
    std::stringstream ss(line);
    int i, j;
    double x[D], v[D];
    ss >> i;
    ss >> j;
    ss >> x[X];
    ss >> x[Y];
    ss >> x[Z];
    ss >> v[X];
    ss >> v[Y];
    ss >> v[Z];
    AddParticle(x, v);
  }
  mout << "# Loaded  " << filename << std::endl;
}
//----------------------------------------------------------------------
void
MDUnit::SaveDensityBinary(double gridsize) {
  CheckPairList();
  static int index = 0;
  char filename[256];
  sprintf(filename, "%s/conf%03d.density", sinfo->BaseDir.c_str(), index);
  index++;
  const unsigned long int tn = obs->GetTotalParticleNumber(vars);

  const int gx = static_cast<int>(sinfo->L[X] / gridsize);
  const int gy = static_cast<int>(sinfo->L[Y] / gridsize);
  const int gz = static_cast<int>(sinfo->L[Z] / gridsize);

  if (0 == minfo->GetRank()) {
    std::ofstream fs(filename, std::ios::out | std::ios::binary);
    fs.write((char *)&tn, sizeof(int));
    fs.write((char *)&sinfo->L[X], sizeof(double));
    fs.write((char *)&sinfo->L[Y], sizeof(double));
    fs.write((char *)&sinfo->L[Z], sizeof(double));
    fs.write((char *)&gx, sizeof(int));
    fs.write((char *)&gy, sizeof(int));
    fs.write((char *)&gz, sizeof(int));
    fs.write((char *)&vars->SimulationTime, sizeof(double));
    fs.close();
  }

  for (int r = 0; r < minfo->GetProcessNumber(); r++) {
    Communicator::Barrier();
    if (minfo->GetRank() != r)continue;
    std::ofstream fs(filename, std::ios::app | std::ios::binary);
    vars->SaveDensityBinary(fs, gridsize, sinfo, myrect);
    fs.close();
  }
}
//----------------------------------------------------------------------
void
MDUnit::SaveDensityBinaryEach(double gridsize) {
  CheckPairList();
  static int index = 0;
  char filename[256];
  sprintf(filename, "%s/conf%03d_%04d.density", sinfo->BaseDir.c_str(), index, minfo->GetRank());
  index++;
  const int tn = obs->GetTotalParticleNumber(vars);

  const int gx = static_cast<int>(sinfo->L[X] / gridsize);
  const int gy = static_cast<int>(sinfo->L[Y] / gridsize);
  const int gz = static_cast<int>(sinfo->L[Z] / gridsize);

  if (0 == minfo->GetRank()) {
    std::ofstream fs(filename, std::ios::out | std::ios::binary);
    fs.write((char *)&tn, sizeof(int));
    fs.write((char *)&sinfo->L[X], sizeof(double));
    fs.write((char *)&sinfo->L[Y], sizeof(double));
    fs.write((char *)&sinfo->L[Z], sizeof(double));
    fs.write((char *)&gx, sizeof(int));
    fs.write((char *)&gy, sizeof(int));
    fs.write((char *)&gz, sizeof(int));
    fs.write((char *)&vars->SimulationTime, sizeof(double));
    fs.close();
  }
  std::ofstream fs(filename, std::ios::binary);
  vars->SaveDensityBinary(fs, gridsize, sinfo, myrect);
  fs.close();
}
//----------------------------------------------------------------------
void
MDUnit::SaveDensity(double gridsize) {
  CheckPairList();
  static int index = 0;
  char filename[256];
  sprintf(filename, "conf/conf%03d.density", index);
  index++;
  const int tn = obs->GetTotalParticleNumber(vars);

  const int gx = static_cast<int>(sinfo->L[X] / gridsize);
  const int gy = static_cast<int>(sinfo->L[Y] / gridsize);
  const int gz = static_cast<int>(sinfo->L[Z] / gridsize);

  if (0 == minfo->GetRank()) {
    std::ofstream fs(filename);
    if (!fs) {
      std::cout << "Could not open " << filename << std::endl;
    }
    fs << "# N=" << tn;
    fs << ",LX=" << sinfo->L[X];
    fs << ",LY=" << sinfo->L[Y];
    fs << ",LZ=" << sinfo->L[Z];
    fs << ",GX=" << gx;
    fs << ",GY=" << gy;
    fs << ",GZ=" << gz;
    fs << ",Time=" << vars->SimulationTime;
    fs << std::endl;
    fs.close();
  }
  for (int r = 0; r < minfo->GetProcessNumber(); r++) {
    Communicator::Barrier();
    if (minfo->GetRank() != r)continue;
    std::ofstream fs(filename, std::ios::app);
    vars->SaveDensity(fs, gridsize, sinfo, myrect);
    fs.close();
  }
}
//----------------------------------------------------------------------
void
MDUnit::SaveDumpFileSequential(void) {
  static int index = 0;
  char filename[256];
  sprintf(filename, "%s/conf%03d.dmp", sinfo->BaseDir.c_str(), index);
  index++;
  SaveDumpFile(filename);
}
//----------------------------------------------------------------------
void
MDUnit::SaveDumpFile(const char *filename) {
  CheckPairList();
  const int pn = obs->GetTotalParticleNumber(vars);
  if (0 == minfo->GetRank()) {
    std::ofstream fs(filename, std::ios::binary);
    double st = vars->SimulationTime;
    double zeta = vars->Zeta;
    fs.write((const char*)&pn, sizeof(pn));
    fs.write((const char*)&st, sizeof(st));
    fs.write((const char*)&zeta, sizeof(zeta));
    fs.close();
  }
  for (int r = 0; r < minfo->GetProcessNumber(); r++) {
    Communicator::Barrier();
    if (minfo->GetRank() != r)continue;
    std::ofstream fs(filename, std::ios::app | std::ios::binary);
    vars->SaveToStream(fs);
    fs.close();
  }
}
//----------------------------------------------------------------------
void
MDUnit::LoadDumpFile(const char *filename) {
  double x[D], v[D];
  std::ifstream ifs(filename, std::ios::binary);
  int pn;
  double st;
  double zeta;
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
    }
    AddParticle(x, v);
  }
}
//----------------------------------------------------------------------
void
MDUnit::StartMeasurement(void) {
  mesh->ClearNumberOfConstructions();
  measure->Start();
}
//----------------------------------------------------------------------
void
MDUnit::StopMeasurement(void) {
  measure->Stop();
}
//----------------------------------------------------------------------
void
MDUnit::ShowMUPS(const int LOOP) {
  const unsigned long int tn = obs->GetTotalParticleNumber(vars);
  measure->ShowMUPS(tn, LOOP);
  mout << "# Number of Pair-list Construction: " << mesh->GetNumberOfConstructions() << std::endl;
}
//----------------------------------------------------------------------
void
MDUnit::DumpPairList(void) {
  const int pn = vars->GetParticleNumber();
  const int *sorted_list = mesh->GetSortedList();
  for (int i = 0; i < pn; i++) {
    const int np = mesh->GetPartnerNumber(i);
    const int kp = mesh->GetKeyPointer(i);
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      printf("%d %d\n", i, j);
    }
  }
}
//----------------------------------------------------------------------
