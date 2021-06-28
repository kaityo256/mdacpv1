//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include "fcalculator.h"
#include "mt.h"
//----------------------------------------------------------------------
void
ForceCalculator::UpdatePositionHalf(Variables *vars, SimulationInfo *sinfo) {
  const double dt2 = sinfo->TimeStep * 0.5;
  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  for (int i = 0; i < pn; i++) {
    q[i][X] += p[i][X] * dt2;
    q[i][Y] += p[i][Y] * dt2;
    q[i][Z] += p[i][Z] * dt2;
  }
}
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForce(Variables *vars, MeshList *mesh, SimulationInfo *sinfo, MPIInfo *minfo) {

#ifdef POWER
  CalculateForceUnroll(vars, mesh, sinfo);
#else
  CalculateForceNext(vars, mesh, sinfo);
#endif
  if (!sinfo->IsPeriodic) {
    CalculateForceOnEdge(vars, sinfo, minfo);
  }
}
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
#ifdef POWER
  CalculateForceUnrollBenchmark(vars, mesh, sinfo);
#else
  CalculateForceNextBenchmark(vars, mesh, sinfo);
  //CalculateForceSortedBenchmark(vars,mesh,sinfo);
#endif
}
//----------------------------------------------------------------------
// Calculate Force for general cases
// Unrolling by double
// Subscript a and b denotes even and odd indexes in the unrolled loop.
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceUnroll(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double C2_8 = C2 * 8.0;
  const double dt = sinfo->TimeStep;


  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    for (int k = kp; k < (kp + np - 1); k += 2) {
      const int j_a = sorted_list[k];
      const int j_b = sorted_list[k + 1];
      const double dx_a = q[j_a][X] - qx_key;
      const double dy_a = q[j_a][Y] - qy_key;
      const double dz_a = q[j_a][Z] - qz_key;

      const double dx_b = q[j_b][X] - qx_key;
      const double dy_b = q[j_b][Y] - qy_key;
      const double dz_b = q[j_b][Z] - qz_key;

      const double r2_a = (dx_a * dx_a + dy_a * dy_a + dz_a * dz_a) * sigma_inv2;
      const double r2_b = (dx_b * dx_b + dy_b * dy_b + dz_b * dz_b) * sigma_inv2;

      const double r6_a = r2_a * r2_a * r2_a;
      const double r6_b = r2_b * r2_b * r2_b;

      const double r14_a = r6_a * r6_a * r2_a;
      const double r14_b = r6_b * r6_b * r2_b;

      const double r14_inv_ab = 1.0 / (r14_a * r14_b);
      double df_a = (24.0 * r6_a - 48.0) * r14_inv_ab * r14_b;
      double df_b = (24.0 * r6_b - 48.0) * r14_inv_ab * r14_a;

      df_a = (df_a + C2_8) * dt * sigma_inv2;
      df_b = (df_b + C2_8) * dt * sigma_inv2;

      if (r2_a > CL2) {
        df_a = 0.0;
      }

      if (r2_b > CL2) {
        df_b = 0.0;
      }

      pfx += df_a * dx_a;
      pfy += df_a * dy_a;
      pfz += df_a * dz_a;

      pfx += df_b * dx_b;
      pfy += df_b * dy_b;
      pfz += df_b * dz_b;

      p[j_a][X] -= df_a * dx_a;
      p[j_a][Y] -= df_a * dy_a;
      p[j_a][Z] -= df_a * dz_a;

      p[j_b][X] -= df_b * dx_b;
      p[j_b][Y] -= df_b * dy_b;
      p[j_b][Z] -= df_b * dz_b;
    }
    if (np % 2 == 1) {
      int j = sorted_list[kp + np - 1];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz) * sigma_inv2;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt * sigma_inv2;
      if (r2 > CL2) {
        df = 0.0;
      }
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
// CalculateForce For Benchmark (Optimized for Power )
// Assuming that sigma is always unity.
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceUnrollBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double C2_8 = C2 * 8.0;
  const double dt = sinfo->TimeStep;

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    for (int k = kp; k < (kp + np - 1); k += 2) {
      const int j_a = sorted_list[k];
      const int j_b = sorted_list[k + 1];
      const double dx_a = q[j_a][X] - qx_key;
      const double dy_a = q[j_a][Y] - qy_key;
      const double dz_a = q[j_a][Z] - qz_key;

      const double dx_b = q[j_b][X] - qx_key;
      const double dy_b = q[j_b][Y] - qy_key;
      const double dz_b = q[j_b][Z] - qz_key;

      const double r2_a = (dx_a * dx_a + dy_a * dy_a + dz_a * dz_a);
      const double r2_b = (dx_b * dx_b + dy_b * dy_b + dz_b * dz_b);

      const double r6_a = r2_a * r2_a * r2_a;
      const double r6_b = r2_b * r2_b * r2_b;

      const double r14_a = r6_a * r6_a * r2_a;
      const double r14_b = r6_b * r6_b * r2_b;

      const double r14_inv_ab = 1.0 / (r14_a * r14_b);
      double df_a = (24.0 * r6_a - 48.0) * r14_inv_ab * r14_b;
      double df_b = (24.0 * r6_b - 48.0) * r14_inv_ab * r14_a;

      df_a = (df_a + C2_8) * dt;
      df_b = (df_b + C2_8) * dt;

      if (r2_a > CL2) {
        df_a = 0.0;
      }

      if (r2_b > CL2) {
        df_b = 0.0;
      }

      pfx += df_a * dx_a;
      pfy += df_a * dy_a;
      pfz += df_a * dz_a;

      pfx += df_b * dx_b;
      pfy += df_b * dy_b;
      pfz += df_b * dz_b;

      p[j_a][X] -= df_a * dx_a;
      p[j_a][Y] -= df_a * dy_a;
      p[j_a][Z] -= df_a * dz_a;

      p[j_b][X] -= df_b * dx_b;
      p[j_b][Y] -= df_b * dy_b;
      p[j_b][Z] -= df_b * dz_b;
    }
    if (np % 2 == 1) {
      int j = sorted_list[kp + np - 1];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      if (r2 > CL2) {
        df = 0.0;
      }
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
// Calculate Force using sorted list.
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceSorted(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {

  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz) * sigma_inv2;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt * sigma_inv2;
      if (r2 > CL2) {
        df = 0.0;
      }
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
// Calculate Force using sorted list.
// Assuming sigma = 1.0
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceSortedBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;

  double (*q2)[D] = vars->q;
  double (*p2)[D] = vars->p;
  typedef double mdvar_t[N][D];
  static mdvar_t q, p;

  const int pn = vars->GetParticleNumber();
  for (int i = 0; i < N; i++) {
    q[i][X] = q2[i][X];
    q[i][Y] = q2[i][Y];
    q[i][Z] = q2[i][Z];
    p[i][X] = p2[i][X];
    p[i][Y] = p2[i][Y];
    p[i][Z] = p2[i][Z];
  }

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      if (r2 > CL2) {
        df = 0.0;
      }
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
  for (int i = 0; i < N; i++) {
    q2[i][X] = q[i][X];
    q2[i][Y] = q[i][Y];
    q2[i][Z] = q[i][Z];
    p2[i][X] = p[i][X];
    p2[i][Y] = p[i][Y];
    p2[i][Z] = p[i][Z];
  }

}
//----------------------------------------------------------------------
// Calculate Force for Benchmark (Optimized for Intel)
// Calculate Next Pair on previous loop
// Assuming that sigma is always unity.
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceNextBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {

  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2() * 8.0;
  const double dt = sinfo->TimeStep;

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    int ja = sorted_list[kp];
    double dxa = q[ja][X] - qx_key;
    double dya = q[ja][Y] - qy_key;
    double dza = q[ja][Z] - qz_key;
    double df = 0.0;
    double dxb = 0.0, dyb = 0.0, dzb = 0.0;
    int jb = 0;

    const int np = mesh->GetPartnerNumber(i);
    for (int k = kp; k < np + kp; k++) {

      const double dx = dxa;
      const double dy = dya;
      const double dz = dza;
      double r2 = (dx * dx + dy * dy + dz * dz);
      const int j = ja;
      ja = sorted_list[k + 1];
      dxa = q[ja][X] - qx_key;
      dya = q[ja][Y] - qy_key;
      dza = q[ja][Z] - qz_key;
      if (r2 > CL2)continue;
      pfx += df * dxb;
      pfy += df * dyb;
      pfz += df * dzb;
      p[jb][X] -= df * dxb;
      p[jb][Y] -= df * dyb;
      p[jb][Z] -= df * dzb;
      const double r6 = r2 * r2 * r2;
      df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2) * dt;
      jb = j;
      dxb = dx;
      dyb = dy;
      dzb = dz;
    }
    p[jb][X] -= df * dxb;
    p[jb][Y] -= df * dyb;
    p[jb][Z] -= df * dzb;
    p[i][X] += pfx + df * dxb;
    p[i][Y] += pfy + df * dyb;
    p[i][Z] += pfz + df * dzb;
  }
}
//----------------------------------------------------------------------
// Calculate Force using sorted list.
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceNext(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {

  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    int ja = sorted_list[kp];
    double dxa = q[ja][X] - qx_key;
    double dya = q[ja][Y] - qy_key;
    double dza = q[ja][Z] - qz_key;
    double df = 0.0;
    double dxb = 0.0, dyb = 0.0, dzb = 0.0;
    int jb = ja;

    for (int k = kp; k < np + kp; k++) {

      const double dx = dxa;
      const double dy = dya;
      const double dz = dza;
      double r2 = (dx * dx + dy * dy + dz * dz) * sigma_inv2;
      const int j = ja;
      ja = sorted_list[k + 1];
      dxa = q[ja][X] - qx_key;
      dya = q[ja][Y] - qy_key;
      dza = q[ja][Z] - qz_key;
      if (r2 > CL2)continue;
      pfx += df * dxb;
      pfy += df * dyb;
      pfz += df * dzb;
      p[jb][X] -= df * dxb;
      p[jb][Y] -= df * dyb;
      p[jb][Z] -= df * dzb;
      double r6 = r2 * r2 * r2;
      df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt * sigma_inv2;
      jb = j;
      dxb = dx;
      dyb = dy;
      dzb = dz;
    }
    p[jb][X] -= df * dxb;
    p[jb][Y] -= df * dyb;
    p[jb][Z] -= df * dzb;
    p[i][X] += pfx + df * dxb;
    p[i][Y] += pfy + df * dyb;
    p[i][Z] += pfz + df * dzb;
  }

}
//----------------------------------------------------------------------
// Calculate Force without optimization
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForcePair(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const double sigma = vars->GetSigma();
  const double sigma_inv2 = 1.0 / (sigma * sigma);

  int *key_index = mesh->GetKeyParticles();
  int *partner_index = mesh->GetPartnerParticles();

  const int number_of_pairs = mesh->GetPairNumber();
  for (int k = 0; k < number_of_pairs; k++) {
    int i = key_index[k];
    int j = partner_index[k];
    double dx = q[j][X] - q[i][X];
    double dy = q[j][Y] - q[i][Y];
    double dz = q[j][Z] - q[i][Z];
    double r2 = (dx * dx + dy * dy + dz * dz) * sigma_inv2;
    if (r2 > CL2)continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt * sigma_inv2;
    p[i][X] += df * dx;
    p[i][Y] += df * dy;
    p[i][Z] += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceOnEdge(Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo) {

  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;

  const double C2 = vars->GetC2();

  int dir1[D] = {D_LEFT, D_BACK, D_DOWN};
  int dir2[D] = {D_RIGHT, D_FORWARD, D_UP};

  const double dt = sinfo->TimeStep;
  for (int d = 0; d < D; d++) {
    if (minfo->IsEdge(dir1[d])) {
      for (int i = 0; i < pn; i++) {
        if (q[i][d] > CUTOFF_LENGTH) continue;
        double dx = q[i][d];
        double r2 = dx * dx;
        double df = -4.0 * (12.0 - 6.0 * r2 * r2 * r2) / (r2 * r2 * r2 * r2 * r2 * r2 * r2) + C2 * 8.0;
        df *= dt;
        p[i][d] -= df * dx;
      }
    }
    if (minfo->IsEdge(dir2[d])) {
      for (int i = 0; i < pn; i++) {
        if (sinfo->L[d] - q[i][d] > CUTOFF_LENGTH) continue;
        double dx = sinfo->L[d] - q[i][d];
        double r2 = dx * dx;
        double df = -4.0 * (12.0 - 6.0 * r2 * r2 * r2) / (r2 * r2 * r2 * r2 * r2 * r2 * r2) + C2 * 8.0;
        df *= dt;
        p[i][d] += df * dx;
      }
    }
  }
}
//----------------------------------------------------------------------
// For Heatbath
//----------------------------------------------------------------------
void
ForceCalculator::HeatbathZeta(double &zeta, double current_temperature, SimulationInfo *sinfo) {
  const double dt2 = sinfo->TimeStep * 0.5;
  const double tau = sinfo->HeatbathTau;
  double t1 = (current_temperature - sinfo->AimedTemperature) / (tau * tau);
  zeta += t1 * dt2;
}
//----------------------------------------------------------------------
void
ForceCalculator::HeatbathMomenta(Variables *vars, SimulationInfo *sinfo) {
  const double dt2 = sinfo->TimeStep * 0.5;
  const int pn = vars->GetParticleNumber();
  double (*p)[D] = vars->p;

  const double exp1 = exp(-dt2 * vars->Zeta);
  for (int i = 0; i < pn; i++) {
    for (int d = 0; d < D; d++) {
      p[i][d] *= exp1;
    }
  }
}
//----------------------------------------------------------------------
void
ForceCalculator::Langevin(Variables *vars, SimulationInfo *sinfo) {
  const double dt = sinfo->TimeStep;
  const int pn = vars->GetParticleNumber();
  double (*p)[D] = vars->p;
  const double hb_gamma = sinfo->HeatbathGamma;
  const double T = sinfo->AimedTemperature;
  const double hb_D = sqrt(2.0 * hb_gamma * T / dt);
  for (int i = 0; i < pn; i++) {
    for (int d = 0; d < D; d++) {
      const double r = MT::GetGauss() * hb_D;
      p[i][d] += (-hb_gamma * p[i][d] + r) * dt;
    }
  }
}
//----------------------------------------------------------------------
