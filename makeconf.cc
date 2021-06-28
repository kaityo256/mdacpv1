//----------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <algorithm>
#include "makeconf.h"
//----------------------------------------------------------------------
// Simple Cubic for Free
//----------------------------------------------------------------------
void
MakeConfigurationFree(MDUnit *mdu, Parameter &param) {
  double density = 0.5;
  double s = 1.0 / pow(density, 1.0 / 3.0);
  double *L = mdu->GetSystemSize();
  double wx = L[X] - CUTOFF_LENGTH * 2.0;
  double wy = L[Y] - CUTOFF_LENGTH * 2.0;
  double wz = L[Z] - CUTOFF_LENGTH * 2.0;
  const int sx = static_cast<int>(wx / s);
  const int sy = static_cast<int>(wy / s);
  const int sz = static_cast<int>(wz / s);
  double x[D];
  for (int ix = 0; ix < sx; ix++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int iz = 0; iz < sz; iz++) {
        x[X] = ix * s + CUTOFF_LENGTH;
        x[Y] = iy * s + CUTOFF_LENGTH;
        x[Z] = iz * s + CUTOFF_LENGTH;
        mdu->AddParticle(x);
      }
    }
  }
}
//----------------------------------------------------------------------
// Simple Cubic for Periodic
//----------------------------------------------------------------------
void
MakeConfigurationPeriodic(MDUnit *mdu, Parameter &param) {
  double q[D];
  double *L = mdu->GetSystemSize();
  double density = param.GetDoubleDef("Density", 0.5);;
  double s = 1.0 / pow(density, 1.0 / 3.0);
  int sx = static_cast<int>(L[X] / s);
  int sy = static_cast<int>(L[Y] / s);
  int sz = static_cast<int>(L[Z] / s);
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
}
//----------------------------------------------------------------------
void
MakeConfigurationFCC(MDUnit *mdu, Parameter &param) {
  double *L = mdu->GetSystemSize();
  const double density = param.GetDoubleDef("Density", 0.5);;
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L[X] / s);
  int sy = static_cast<int>(L[Y] / s);
  int sz = static_cast<int>(L[Z] / s);
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
}
//----------------------------------------------------------------------
void
MakeConfigurationFCC2(MDUnit *mdu, Parameter &param) {
  double *L = mdu->GetSystemSize();
  const double density = param.GetDoubleDef("Density", 0.5);;
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
// Make Configuration FCC with accurate density
// Prepare the system with higher density, and remove some particles
// in order to adjust density
//----------------------------------------------------------------------
void
MakeConfigurationExtract(MDUnit *mdu, Parameter &param) {
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
//----------------------------------------------------------------------
