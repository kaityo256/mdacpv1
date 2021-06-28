//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include "variables.h"
#include "mt.h"

//----------------------------------------------------------------------
Variables::Variables(void) {
  particle_number = 0;
  total_particle_number = 0;
  SimulationTime = 0.0;
  sigma = 1.0;
  Zeta = 0.0;
  const double s2 = 1.0 / (CUTOFF_LENGTH * CUTOFF_LENGTH);
  const double s6 = s2 * s2 * s2;
  const double s8 = s6 * s2;
  const double s12 = s6 * s6;
  const double s14 = s12 * s2;
  C2 = 6.0 * s14 - 3.0 * s8;
  C0 = -s12 + s6 - C2 / s2;
}
//----------------------------------------------------------------------
void
Variables::AddParticle(double x[D], double v[D]) {
  q[particle_number][X] = x[X];
  q[particle_number][Y] = x[Y];
  q[particle_number][Z] = x[Z];
  p[particle_number][X] = v[X];
  p[particle_number][Y] = v[Y];
  p[particle_number][Z] = v[Z];
  particle_number++;
  total_particle_number = particle_number;
  assert(particle_number < N);
}
//----------------------------------------------------------------------
void
Variables::SaveConfiguration(MPIInfo *minfo) {
  char filename[256];
  sprintf(filename, "conf%02d.dat", minfo->GetRank());
  std::ofstream fs(filename);
  for (int i = 0; i < total_particle_number; i++) {
    fs << q[i][X];
    fs << " " << q[i][Y];
    fs << " " << q[i][Z];
    fs << " " << p[i][X];
    fs << " " << p[i][Y];
    fs << " " << p[i][Z];
    fs << std::endl;
  }
}
//----------------------------------------------------------------------
void
Variables::SaveDensityBinary(std::ostream &fs, double gridsize, SimulationInfo *sinfo, MDRect myrect) {
  const int pn = GetParticleNumber();
  const int gn = static_cast<int>(myrect.GetWidth(X) / gridsize);
  const int gridnum = gn * gn * gn;
  double *grid = new double[gridnum];
  for (int i = 0; i < gridnum; i++) {
    grid[i] = 0.0;
  }
  const double gsinv = 1.0 / gridsize;
  double *s = myrect.GetStartPosition();
  const double vinv = 1.0 / (gridsize * gridsize * gridsize);
  for (int i = 0; i < pn; i++) {
    int ix = static_cast<int>((q[i][X] - s[X]) * gsinv);
    int iy = static_cast<int>((q[i][Y] - s[Y]) * gsinv);
    int iz = static_cast<int>((q[i][Z] - s[Z]) * gsinv);
    int index = ix + iy * gn + iz * gn * gn;
    grid[index] += vinv;
  }
  const int sx = static_cast<int>(s[X] / gridsize);
  const int sy = static_cast<int>(s[Y] / gridsize);
  const int sz = static_cast<int>(s[Z] / gridsize);
  const int glx = static_cast<int>(sinfo->L[X] / gridsize);
  const int gly = static_cast<int>(sinfo->L[Y] / gridsize);

  const unsigned int bufsize = gn * gn * gn * (sizeof(int) + sizeof(double));
  char *buf = new char[bufsize];
  int pos = 0;
  for (int iz = 0; iz < gn; iz++) {
    for (int iy = 0; iy < gn; iy++) {
      for (int ix = 0; ix < gn; ix++) {
        int index = ix + iy * gn + iz * gn * gn;
        int gindex = (ix + sx) + (iy + sy) * glx + (iz + sz) * glx * gly;
        memcpy((&buf[pos]), &gindex, sizeof(int));
        pos += sizeof(int);
        memcpy((&buf[pos]), &grid[index], sizeof(double));
        pos += sizeof(double);
      }
    }
  }
  fs.write(buf, bufsize);
  delete [] buf;
  delete [] grid;
}
//----------------------------------------------------------------------
void
Variables::SaveDensity(std::ostream &fs, double gridsize, SimulationInfo *sinfo, MDRect myrect) {
  const int pn = GetParticleNumber();
  const int gn = static_cast<int>(myrect.GetWidth(X) / gridsize);
  const int gridnum = gn * gn * gn;
  double *grid = new double[gridnum];
  for (int i = 0; i < gridnum; i++) {
    grid[i] = 0.0;
  }
  const double gsinv = 1.0 / gridsize;
  double *s = myrect.GetStartPosition();
  const double vinv = 1.0 / (gridsize * gridsize * gridsize);
  for (int i = 0; i < pn; i++) {
    int ix = static_cast<int>((q[i][X] - s[X]) * gsinv);
    int iy = static_cast<int>((q[i][Y] - s[Y]) * gsinv);
    int iz = static_cast<int>((q[i][Z] - s[Z]) * gsinv);
    int index = ix + iy * gn + iz * gn * gn;
    grid[index] += vinv;
  }
  const int sx = static_cast<int>(s[X] / gridsize);
  const int sy = static_cast<int>(s[Y] / gridsize);
  const int sz = static_cast<int>(s[Z] / gridsize);
  const int glx = static_cast<int>(sinfo->L[X] / gridsize);
  const int gly = static_cast<int>(sinfo->L[Y] / gridsize);

  for (int iz = 0; iz < gn; iz++) {
    for (int iy = 0; iy < gn; iy++) {
      for (int ix = 0; ix < gn; ix++) {
        int index = ix + iy * gn + iz * gn * gn;
        int gindex = (ix + sx) + (iy + sy) * glx + (iz + sz) * glx * gly;
        fs << gindex << " " << grid[index] << "\n";
      }
    }
  }
  delete [] grid;
}
//----------------------------------------------------------------------
void
Variables::SaveConfiguration(std::ostream &fs) {
  for (int i = 0; i < particle_number; i++) {
    fs << q[i][X];
    fs << " " << q[i][Y];
    fs << " " << q[i][Z];
    fs << " " << p[i][X];
    fs << " " << p[i][Y];
    fs << " " << p[i][Z];
    fs << "\n";
  }
}
//----------------------------------------------------------------------
void
Variables::AdjustPeriodicBoundary(SimulationInfo *sinfo) {
  for (int i = 0; i < particle_number; i++) {
    for (int d = 0; d < D; d++) {
      if (q[i][d] < 0.0) {
        q[i][d] += sinfo->L[d];
      } else if (q[i][d] > sinfo->L[d]) {
        q[i][d] -= sinfo->L[d];
      }
    }
  }

  for (int i = 0; i < particle_number; i++) {
    for (int d = 0; d < D; d++) {
      if (q[i][d] < 0.0 || q[i][d] > sinfo->L[d]) {
        show_error("Invalid position");
        printf("%f: %d %d %f\n", SimulationTime, i, d, q[i][d]);
        exit(1);
      }
    }
  }
}
//----------------------------------------------------------------------
void
Variables::SetInitialVelocity(const double V0) {
  for (int i = 0; i < particle_number; i++) {
    double z = MT::GetDouble() * 2.0 - 1.0;
    double phi = 2.0 * MT::GetDouble() * M_PI;
    p[i][X] = V0 * sqrt(1 - z * z) * cos(phi);
    p[i][Y] = V0 * sqrt(1 - z * z) * sin(phi);
    p[i][Z] = V0 * z;
  }
  double avx = 0;
  double avy = 0;
  double avz = 0;
  for (int i = 0; i < particle_number; i++) {
    avx += p[i][X];
    avy += p[i][Y];
    avz += p[i][Z];
  }
  avx /= static_cast<double>(particle_number);
  avy /= static_cast<double>(particle_number);
  avz /= static_cast<double>(particle_number);
  for (int i = 0; i < particle_number; i++) {
    p[i][X] -= avx;
    p[i][Y] -= avy;
    p[i][Z] -= avz;
  }
}
//----------------------------------------------------------------------
void
Variables::SaveToStream(std::ostream &fs) {
  for (int i = 0; i < particle_number; i++) {
    for (int d = 0; d < D; d++) {
      fs.write((const char*)&q[i][d], sizeof(double));
      fs.write((const char*)&p[i][d], sizeof(double));
    }
  }
}
//----------------------------------------------------------------------
void
Variables::LoadFromStream(std::istream &fs) {
  fs.read((char*)&particle_number, sizeof(particle_number));
  fs.read((char*)&SimulationTime, sizeof(SimulationTime));
  fs.read((char*)&Zeta, sizeof(Zeta));
  for (int i = 0; i < particle_number; i++) {
    for (int d = 0; d < D; d++) {
      fs.read((char*)&q[i][d], sizeof(double));
      fs.read((char*)&p[i][d], sizeof(double));
    }
  }
}
//----------------------------------------------------------------------
// For Debug
//----------------------------------------------------------------------
void
Variables::Shuffle(void) {
  const int pn = particle_number;
  int *a = new int[pn];
  for (int i = 0; i < pn; i++) {
    a[i] = i;
  }
  std::random_shuffle(&a[0], &a[pn]);
  double (*qb)[D] = new double[pn][D];
  double (*pb)[D] = new double[pn][D];

  for (int i = 0; i < pn; i++) {
    for (int d = 0; d < D; d++) {
      qb[i][d] = q[i][d];
      pb[i][d] = p[i][d];
    }
  }

  for (int i = 0; i < pn; i++) {
    int j = a[i];
    for (int d = 0; d < D; d++) {
      q[j][d] = qb[i][d];
      p[j][d] = pb[i][d];
    }
  }
  delete [] qb;
  delete [] pb;
  delete [] a;
}
//----------------------------------------------------------------------
