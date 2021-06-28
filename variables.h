//----------------------------------------------------------------------
#ifndef variables_h
#define variables_h
//----------------------------------------------------------------------
#include <iostream>
#include "mdconfig.h"
#include "simulationinfo.h"
#include "mpiinfo.h"
#include "mdrect.h"
//----------------------------------------------------------------------
class Variables {
private:
  int particle_number;
  int total_particle_number;
  double C0, C2;
  double sigma;
public:
  double q[N][D];
  double p[N][D];
  double Zeta;
  double SimulationTime;
  double GetC0(void) {return C0;};
  double GetC2(void) {return C2;};
  double GetSigma(void) {return sigma;};
  void SetSigma(double s) {sigma = s;};
  int GetParticleNumber(void) {return particle_number;};
  int GetTotalParticleNumber(void) {return total_particle_number;};

  Variables(void);
  void AddParticle(double x[D], double v[D]);
  void SaveToStream(std::ostream &fs);
  void LoadFromStream(std::istream &fs);
  void SaveConfiguration(MPIInfo *minfo);
  void SaveConfiguration(std::ostream &fs);
  void SaveDensity(std::ostream &fs, double gridsize, SimulationInfo *sinfo, MDRect myrect);
  void SaveDensityBinary(std::ostream &fs, double gridsize, SimulationInfo *sinfo, MDRect myrect);
  void SetParticleNumber(int pn) {particle_number = pn;};
  void SetTotalParticleNumber(int pn) {total_particle_number = pn;};
  void AdjustPeriodicBoundary(SimulationInfo *sinfo);
  void SetInitialVelocity(double V0);
  // For Debug
  void Shuffle(void);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
