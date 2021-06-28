#include <iostream>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <functional>
#include <algorithm>
#include "bubblehist.h"
#include "communicator.h"
//----------------------------------------------------------------------
//const double density_threshold = 0.4;
const double density_threshold = 0.2;
//----------------------------------------------------------------------
BubbleHist::BubbleHist(int r, int p, MDRect rect, double *ssize) {
  rank = r;
  proc_num = p;
  myrect = rect;
  L = ssize;
  num_observation = 0;
  const double GRIDSIZE = 2.0;
  double wx = myrect.GetWidth(X);
  double wy = myrect.GetWidth(Y);
  double wz = myrect.GetWidth(Z);

  local_grid_x = static_cast<int>(wx / GRIDSIZE);
  local_grid_y = static_cast<int>(wy / GRIDSIZE);
  local_grid_z = static_cast<int>(wz / GRIDSIZE);
  grid_size = wx / static_cast<double>(local_grid_x);

  gnx = static_cast<int>(L[X] / grid_size + 0.5);
  gny = static_cast<int>(L[Y] / grid_size + 0.5);
  gnz = static_cast<int>(L[Z] / grid_size + 0.5);

  local_grid_num = local_grid_x * local_grid_y * local_grid_z;
  grid_num = gnx * gny * gnz;
  density_grid = new double[local_grid_num];
  clustering_grid = new bool[grid_num];
  cluster_number = new int[grid_num];
  cluster_size = new int[grid_num];
  size_distribution = new double[grid_num];
  index_list = new int[grid_num];
  for (int i = 0; i < grid_num; i++) {
    size_distribution[i] = 0.0;
  }
}
//----------------------------------------------------------------------
BubbleHist::~BubbleHist(void) {
  delete [] density_grid;
  delete [] clustering_grid;
  delete [] cluster_number;
  delete [] cluster_size;
  delete [] size_distribution;
  delete [] index_list;
}
//----------------------------------------------------------------------
void
BubbleHist::Clear(void) {
  num_observation = 0;
  for (int i = 0; i < grid_num; i++) {
    size_distribution[i] = 0.0;
  }
}
//----------------------------------------------------------------------
int
BubbleHist::GetIndex(int ix, int iy, int iz) {
  if ( ix < 0) {
    ix += gnx;
  } else if (ix >= gnx) {
    ix -= gnx;
  }
  if ( iy < 0) {
    iy += gny;
  } else if (iy >= gny) {
    iy -= gny;
  }
  if ( iz < 0) {
    iz += gnz;
  } else if (iz >= gnz) {
    iz -= gnz;
  }
  return ix + iy * gnx + iz * gnx * gny;
}
//----------------------------------------------------------------------
int
BubbleHist::GetClusterIndex(int index) {
  while (index != cluster_number[index]) {
    index = cluster_number[index];
  }
  return index;
}
//----------------------------------------------------------------------
int
BubbleHist::GetClusterIndex(int ix, int iy, int iz) {
  int index = GetIndex(ix, iy, iz);
  while (index != cluster_number[index]) {
    index = cluster_number[index];
  }
  return index;
}
//----------------------------------------------------------------------
void
BubbleHist::Check(int ix1, int iy1, int iz1, int ix2, int iy2, int iz2) {
  int i1 = GetIndex(ix1, iy1, iz1);
  int i2 = GetIndex(ix2, iy2, iz2);
  if (!clustering_grid[i1] || !clustering_grid[i2]) {
    return;
  }
  int c1 = GetClusterIndex(ix1, iy1, iz1);
  int c2 = GetClusterIndex(ix2, iy2, iz2);
  if (c1 < c2) {
    cluster_number[c2] = c1;
  } else {
    cluster_number[c1] = c2;
  }
}
//----------------------------------------------------------------------
void
BubbleHist::Clustering(void) {
  for (int i = 0; i < grid_num; i++) {
    cluster_number[i] = i;
  }
  for (int iz = 0; iz < gnz; iz++) {
    for (int iy = 0; iy < gny; iy++) {
      for (int ix = 0; ix < gnx; ix++) {
        Check(ix, iy, iz, ix + 1, iy, iz);
        Check(ix, iy, iz, ix, iy + 1, iz);
        Check(ix, iy, iz, ix, iy, iz + 1);
      }
    }
  }
}
//----------------------------------------------------------------------
double
BubbleHist::GetClusterSize(int index) {
  const double v = grid_size * grid_size * grid_size;
  return v * cluster_size[index];
}
//----------------------------------------------------------------------
void
BubbleHist::MakeDensity(Variables *vars) {
  for (int i = 0; i < local_grid_num; i++) {
    density_grid[i] = 0.0;
  }
  const int pn = vars->GetParticleNumber();
  double *s = myrect.GetStartPosition();
  double (*q)[D] = vars->q;

  const double gsinv = 1.0 / grid_size;
  const double vinv = gsinv * gsinv * gsinv;
  for (int i = 0; i < pn; i++) {
    double qx = q[i][X];
    double qy = q[i][Y];
    double qz = q[i][Z];
    if (qx < 0) {
      qx += L[X];
    } else if (qx > L[X]) {
      qx -= L[X];
    }
    if (qy < 0) {
      qy += L[Y];
    } else if (qy > L[Y]) {
      qy -= L[Y];
    }
    if (qz < 0) {
      qz += L[Z];
    } else if (qz > L[Z]) {
      qz -= L[Z];
    }

    int ix = static_cast<int>((qx - s[X]) * gsinv);
    int iy = static_cast<int>((qy - s[Y]) * gsinv);
    int iz = static_cast<int>((qz - s[Z]) * gsinv);
    int local_index = ix + iy * local_grid_x + iz * local_grid_x * local_grid_y;
    assert(local_index >= 0);
    assert(local_index < local_grid_num);
    density_grid[local_index] += vinv;
  }
  for (int i = 0; i < local_grid_num; i++) {
    int ii = i;
    int ix = ii % local_grid_x + static_cast<int>(s[X] * gsinv);
    ii /= local_grid_x;
    int iy = ii % local_grid_y + static_cast<int>(s[Y] * gsinv);
    int iz = ii / local_grid_y + static_cast<int>(s[Z] * gsinv);
    int index = ix + gnx * iy + gnx * gny * iz;
    index_list[i] = index;
  }
}
//----------------------------------------------------------------------
void
BubbleHist::Analyse(Variables *vars) {
  num_observation++;
  for (int i = 0; i < proc_num; i++) {
    if (rank == i) {
      MakeDensity(vars);
    }
    Communicator::Barrier();
  }

  double *recv_density = new double[grid_num];
  int *recv_index = new int[grid_num];

  Communicator::AllGatherInteger(index_list, local_grid_num, recv_index);
  Communicator::AllGatherDouble(density_grid, local_grid_num, recv_density);

  for (int i = 0; i < grid_num; i++) {
    int index = recv_index[i];
    clustering_grid[index] = (recv_density[i] < density_threshold);
  }

  delete [] recv_density;
  delete [] recv_index;

  Clustering();

  for (int i = 0; i < grid_num; i++) {
    cluster_size[i] = 0;
  }
  for (int i = 0; i < grid_num; i++) {
    int index = GetClusterIndex(i);
    cluster_size[index]++;
  }
  std::sort(&cluster_size[0], &cluster_size[grid_num], std::greater<int>());
  const double V = L[X] * L[Y] * L[Z];
  const double Vinv = 1.0 / V;
  for (int i = 0; i < grid_num; i++) {
    int n = cluster_size[i];
    size_distribution[n] += Vinv;
  }
}
//----------------------------------------------------------------------
void
BubbleHist::SaveHistogram(std::string filename) {
  std::ofstream fs(filename.c_str());
  double v = grid_size * grid_size * grid_size;
  double maxsize = GetClusterSize(0);
  const double V = L[X] * L[Y] * L[Z];
  const double Vinv = 1.0 / V;
  fs << "#";
  fs << ",system_volume=" << V;
  fs << ",grid_volume=" << v;
  fs << ",num_observation=" << num_observation;
  fs << std::endl;

  double ob_inv = 1.0 / static_cast<double>(num_observation);
  for (int i = 0; i < grid_num; i++) {
    size_distribution[i] *= ob_inv;
  }

  for (int i = 0; i < grid_num; i++) {
    if (size_distribution[i] > Vinv * ob_inv * 0.5) {
      fs << v*i;
      fs << " " << size_distribution[i];
      fs << " " << log(size_distribution[i]);
      fs << std::endl;
    }
    if (v * i > maxsize) {
      break;
    }
  }
}
//----------------------------------------------------------------------
void
BubbleHist::SaveBubbleDistribution(MDUnit *mdu, const char *base_dir) {
  mdu->CheckPairList();
  static int index = 0;
  char filename[256];
  sprintf(filename, "%s/bubble%03d.dist", base_dir, index);
  index++;
  Clear();
  Analyse(mdu->GetVariables());
  if (0 == rank) {
    double grid_volume = grid_size * grid_size * grid_size;
    int *bubble_num = new int[grid_num];
    const double volume = L[X] * L[Y] * L[Z];
    std::ofstream fs(filename);
    fs << "#";
    fs << ",system_volume=" << volume;
    fs << ",grid_volume=" << grid_volume;
    fs << ",time=" << mdu->GetSimulationTime();
    fs << std::endl;
    for (int i = 0; i < grid_num; i++) {
      bubble_num[i] = 0;
    }
    for (int i = 0; i < grid_num; i++) {
      int n = cluster_size[i];
      bubble_num[n]++;
    }
    int lastsize = 0;
    for (int i = 0; i < grid_num; i++) {
      int n = cluster_size[i];
      if (n == lastsize) {
        continue;
      }
      lastsize = n;
      fs << n*grid_volume << " " << bubble_num[n] << std::endl;
      if (n < 2) {
        break;
      }
    }
    delete [] bubble_num;
  }
}
//----------------------------------------------------------------------
