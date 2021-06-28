//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include "observer.h"
#include "communicator.h"
//----------------------------------------------------------------------
unsigned long int
Observer::GetTotalParticleNumber(Variables *vars) {
  unsigned long int pn = vars->GetParticleNumber();
  return Communicator::AllReduceUnsignedLongInteger(pn);
}
//----------------------------------------------------------------------
double
Observer::Density(Variables *vars, SimulationInfo *sinfo) {
  const double volume = sinfo->L[X] * sinfo->L[Y] * sinfo->L[Z];
  const double pn = static_cast<double>(GetTotalParticleNumber(vars));
  return pn / volume;
}
//----------------------------------------------------------------------
double
Observer::KineticEnergy(Variables *vars) {
  double e = 0;
  const int pn = vars->GetParticleNumber();
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);
  double (*p)[D] = vars->p;
  for (int i = 0; i < pn; i++) {
    for (int d = 0; d < D; d++) {
      e += 0.5 * p[i][d] * p[i][d];
    }
  }
  int tn = GetTotalParticleNumber(vars);
  e = Communicator::AllReduceDouble(e);
  e /= static_cast<double>(tn);
  e *= sigma_inv2;
  return e;
}
//----------------------------------------------------------------------
double
Observer::Temperature(Variables *vars) {
  return KineticEnergy(vars) / 1.5;
}
//----------------------------------------------------------------------
double
Observer::PotentialEnergy(Variables *vars, MeshList *mesh, SimulationInfo *sinfo, MPIInfo *minfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  double (*q)[D] = vars->q;
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);
  const double C2 = vars->GetC2();
  const double C0 = vars->GetC0();
  const int pn = vars->GetParticleNumber();

  double energy = 0.0;

  const int s = mesh->GetPairNumber();
  int *key_particles = mesh->GetKeyParticles();
  int *partner_particles = mesh->GetPartnerParticles();
  for (int k = 0; k < s; k++) {
    int i = key_particles[k];
    int j = partner_particles[k];
    double dx = q[i][X] - q[j][X];
    double dy = q[i][Y] - q[j][Y];
    double dz = q[i][Z] - q[j][Z];
    const double r2 = (dx * dx + dy * dy + dz * dz) * sigma_inv2;
    if (r2 > CL2) continue;
    double e = 4.0 * (1.0 / (r2 * r2 * r2 * r2 * r2 * r2) - 1.0 / (r2 * r2 * r2) + C2 * r2 + C0);
    if (i >= pn || j >= pn ) {
      e *= 0.5;
    }
    energy += e;
  }

  if (!sinfo->IsPeriodic) {
    energy += PotentialEnergyOnEdge(vars, sinfo, minfo);
  }

  int tn = GetTotalParticleNumber(vars);
  energy = Communicator::AllReduceDouble(energy);
  energy /= static_cast<double>(tn);
  return energy;
}
//----------------------------------------------------------------------
double
Observer::PotentialEnergyOnEdge(Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo) {
  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  const double C0 = vars->GetC0();
  const double C2 = vars->GetC2();
  int dir1[D] = {D_LEFT, D_BACK, D_DOWN};
  int dir2[D] = {D_RIGHT, D_FORWARD, D_UP};

  double energy = 0.0;
  for (int d = 0; d < D; d++) {
    if (minfo->IsEdge(dir1[d])) {
      for (int i = 0; i < pn; i++) {
        if (q[i][d] > CUTOFF_LENGTH) continue;
        double dx = q[i][d];
        double r2 = dx * dx;
        double e = 4.0 * (1.0 / (r2 * r2 * r2 * r2 * r2 * r2) - 1.0 / (r2 * r2 * r2) + C2 * r2 + C0);
        energy += e;
      }
    }
    if (minfo->IsEdge(dir2[d])) {
      for (int i = 0; i < pn; i++) {
        if (sinfo->L[d] - q[i][d] > CUTOFF_LENGTH) continue;
        double dx = sinfo->L[d] - q[i][d];
        double r2 = dx * dx;
        double e = 4.0 * (1.0 / (r2 * r2 * r2 * r2 * r2 * r2) - 1.0 / (r2 * r2 * r2) + C2 * r2 + C0);
        energy += e;
      }
    }
  }
  return energy;
}
//----------------------------------------------------------------------
double
Observer::TotalEnergy(Variables *vars, MeshList *mesh, SimulationInfo *sinfo, MPIInfo *minfo) {
  return KineticEnergy(vars) + PotentialEnergy(vars, mesh, sinfo, minfo);
}
//----------------------------------------------------------------------
double
Observer::Pressure(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  double phi = 0;
  double (*q)[D] = vars->q;
  const int pn = vars->GetParticleNumber();
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);
  const double C_2 = vars->GetC2();
  const double C2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const int s = mesh->GetPairNumber();
  int *key_particles = mesh->GetKeyParticles();
  int *partner_particles = mesh->GetPartnerParticles();
  for (int k = 0; k < s; k++) {
    int i = key_particles[k];
    int j = partner_particles[k];
    double dx = q[i][X] - q[j][X];
    double dy = q[i][Y] - q[j][Y];
    double dz = q[i][Z] - q[j][Z];
    const double r2 = (dx * dx + dy * dy + dz * dz) * sigma_inv2;
    if (r2 > C2) continue;
    const double r6 = r2 * r2 * r2;
    double pp = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C_2 * 8.0) * r2 * sigma_inv2;
    if (i >= pn && j >= pn) {
      continue;
    }
    if (i >= pn || j >= pn ) {
      pp *= 0.5;
    }
    phi += pp;
  }
  phi = Communicator::AllReduceDouble(phi);
  double V = sinfo->L[X] * sinfo->L[Y] * sinfo->L[Z];
  const double T = Temperature(vars);
  const int tn = GetTotalParticleNumber(vars);
  return (T * tn - phi / 3.0) / V;
}
//----------------------------------------------------------------------
