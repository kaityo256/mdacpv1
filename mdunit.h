//----------------------------------------------------------------------
#ifndef mdunit_h
#define mdunit_h
//----------------------------------------------------------------------
#include "mpiinfo.h"
#include "mdconfig.h"
#include "variables.h"
#include "meshlist.h"
#include "cmanager.h"
#include "observer.h"
#include "pairlist.h"
#include "measure.h"
#include "parameter.h"
#include "simulationinfo.h"
#include "mt.h"
//----------------------------------------------------------------------
class MDUnit {
private:
  SimulationInfo *sinfo;
  Variables *vars;
  MeshList *mesh;
  Observer *obs;
  PairList *plist;
  Measure *measure;
  CommunicationManager *cmanager;
  MDRect myrect;
  MPIInfo *minfo;
  int grid_position[D];

  bool IsPairListExpired(void);
  void SendParticles(void) {cmanager->SendParticles(vars, sinfo, myrect, minfo);};
  void MakeNeighborInformation(void) {cmanager->MakeNeighborInformation(vars, sinfo, myrect, minfo);};
  void SendNeighborInformation(void) {cmanager->SendNeighborInformation(vars, sinfo, minfo);};
  void CalculateHeatbath(void);

public:
  MDUnit(MPIInfo *minfo, Parameter &param);
  ~MDUnit(void);

  void Calculate(void);
  void CalculateBenchmark(void);
  void MakePairList(void);
  void CheckPairList(void);
  void AddParticle(double x[D], double v[D]);
  void AddParticle(double x[D]);

  // File I/O
  void SaveConfiguration(const char *filename);
  void SaveConfigurationSequential(void);
  void SaveDensity(double gridsize);
  void SaveDensityBinary(double gridsize);
  void SaveDensityBinaryEach(double gridsize);
  void SaveAsCdviewSequential(void);
  void SaveAsCdview(const char *filename);
  void LoadCdview(const char *filename);
  void SaveDumpFileSequential(void);
  void SaveDumpFile(const char *filename);
  void LoadDumpFile(const char *filename);

  MDRect GetRect(void) {return myrect;};

  double Density(void) {return obs->Density(vars, sinfo);};
  double Temperature(void) {return obs->Temperature(vars);};
  double KineticEnergy(void) {return obs->KineticEnergy(vars);};
  double PotentialEnergy(void) {return obs->PotentialEnergy(vars, mesh, sinfo, minfo);};
  double TotalEnergy(void) {return obs->TotalEnergy(vars, mesh, sinfo, minfo);};
  double Pressure(void) {return obs->Pressure(vars, mesh, sinfo);};

  Variables* GetVariables(void) {return vars;};
  void SetInitialVelocity(double v0) {vars->SetInitialVelocity(v0);};
  double GetZeta(void) {return vars->Zeta;};

  void SetSeed(int s) {MT::SetSeed(s);};


  void ChangeScale(double s);
  int GetRank(void) {return minfo->GetRank();};
  int GetProcessNumber(void) {return minfo->GetProcessNumber();};
  double* GetSystemSize(void) {return sinfo->L;};
  void SetPeriodic(bool b) {sinfo->IsPeriodic = b;};
  bool IsPeriodic(void) {return sinfo->IsPeriodic;};
  void SetAimedTemperature(double t) {sinfo->AimedTemperature = t;};
  void ShowSystemInfo(void);
  double GetSimulationTime(void) {return vars->SimulationTime;};
  void SetControlTemperature(bool b) {sinfo->ControlTemperature = b;};
  void SetTimeStep(double t) {sinfo->TimeStep = t;};
  double GetTimeStep(void) {return sinfo->TimeStep;};
  void SetSimulationTime(double t) {vars->SimulationTime = t;};
  int GetParticleNumber(void) {return vars->GetParticleNumber();};

  void StartMeasurement(void);
  void StopMeasurement(void);
  void ShowMUPS(const int LOOP);

  // Mesh Information
  int* GetKeyParticles(void) {return mesh->GetKeyParticles();};
  int* GetPartnerParticles(void) {return mesh->GetPartnerParticles();};
  int GetPairNumber(void) {return mesh->GetPairNumber();};

  //For Debug
  void Shuffle(void) {vars->Shuffle();};
  void DumpPairList(void);

};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
