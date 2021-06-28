//----------------------------------------------------------------------
#include <iostream>
#include <mpi.h>
#include <iomanip>
#include "mpistream.h"
#include "mpiinfo.h"
#include "mdunit.h"
//----------------------------------------------------------------------
MPIInfo::MPIInfo(int *argc, char*** argv) {
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  mout.SetRank(rank);
  mout << "# " << MDACP_VERSION << std::endl;
  mout << "# Number of Processes: " << num_procs << std::endl;
  grid = new int[num_procs];
  qx = new int[num_procs];
  qy = new int[num_procs];
  qz = new int[num_procs];
  for (int i = 0; i < num_procs; i++) {
    grid[i] = i;
  }
  MakeGrid();
  RegisterPosition();
}
//----------------------------------------------------------------------
MPIInfo::~MPIInfo(void) {
  delete [] qx;
  delete [] qy;
  delete [] qz;
  delete [] grid;
  MPI_Finalize();
}
//----------------------------------------------------------------------
bool
MPIInfo::SetGlobalGridSize(Parameter &param) {
  if (!param.Contains("GridX")) {
    return true;
  }
  if (!param.Contains("GridY")) {
    show_error("#GridX is set, but GridY is not set");
    return false;
  }
  if (!param.Contains("GridZ")) {
    show_warning("#GridX is set, but GridZ is not set");
    return false;
  }
  int ix = param.GetInteger("GridX");
  int iy = param.GetInteger("GridY");
  int iz = param.GetInteger("GridZ");
  if (num_procs != ix * iy * iz) {
    show_error("Invalid Grid Numbers");
    mout << "NumProcs = " << num_procs << std::endl;
    mout << "GridX=" << ix;
    mout << " GridY=" << iy;
    mout << " GridZ=" << iz;
    mout << " Total=" << ix*iy*iz << std::endl;
    return false;
  }
  g[X] = ix;
  g[Y] = iy;
  g[Z] = iz;
  return true;
}
//----------------------------------------------------------------------
bool
MPIInfo::SetLocalGridSize(Parameter &param) {
  if (!param.Contains("LocalGridX")) {
    return true;
  }
  if (!param.Contains("LocalGridY")) {
    show_error("#LocalGridX is set, but LocalGridY is not set");
    return false;
  }
  if (!param.Contains("LocalGridZ")) {
    show_warning("#LocalGridX is set, but LocalGridZ is not set");
    return false;
  }
  int sx = param.GetInteger("LocalGridX");
  int sy = param.GetInteger("LocalGridY");
  int sz = param.GetInteger("LocalGridZ");
  if (g[X] % sx != 0) {
    show_error("#LocalGridX is not divisor of GridX");
    return false;
  }
  if (g[Y] % sy != 0) {
    show_error("#LocalGridY is not divisor of GridY");
    return false;
  }
  if (g[Z] % sz != 0) {
    show_error("#LocalGridZ is not divisor of GridZ");
    return false;
  }

  int lx = g[X] / sx;
  int ly = g[Y] / sy;

  for (int i = 0; i < num_procs; i++) {
    int i_global = i / (sx * sy * sz);
    int i_local =  i % (sx * sy * sz);
    int x_global = i_global % lx;
    i_global /= lx;
    int y_global = i_global % ly;
    int z_global = i_global / ly;

    int x_local = i_local % sx;
    i_local /= sx;
    int y_local = i_local % sy;
    int z_local = i_local / sy;

    int ix = x_global * sx + x_local;
    int iy = y_global * sy + y_local;
    int iz = z_global * sz + z_local;
    int index = ix + iy * g[X] + iz * g[X] * g[Y];
    grid[index] = i;
  }
  return true;
}
//----------------------------------------------------------------------
bool
MPIInfo::SetParameter(Parameter &param) {
  bool isValid = true;
  if (!SetGlobalGridSize(param)) {
    isValid = false;
  }
  if (!SetLocalGridSize(param)) {
    isValid = false;
  }
  RegisterPosition();

  if (num_procs != g[X]*g[Y]*g[Z]) {
    show_error("# A total number of cells is invalid.");
    isValid = false;
  }

  return isValid;
}
//----------------------------------------------------------------------
void
MPIInfo::Index2Pos(int index, int &ix, int &iy, int &iz) {
  ix = index % g[X];
  index -= ix;
  index /= g[X];
  iy = index % g[Y];
  iz = (index - iy) / g[Y];
}
//----------------------------------------------------------------------
void
MPIInfo::MakeGrid(void) {
  int ix = 1;
  int iy = 1;
  int iz = 1;
  int np = num_procs;
  while (np > 1) {
    if (np > 1) {
      ix *= 2;
      np /= 2;
    }
    if (np > 1) {
      iy *= 2;
      np /= 2;
    }
    if (np > 1) {
      iz *= 2;
      np /= 2;
    }
  }
  g[X] = ix;
  g[Y] = iy;
  g[Z] = iz;
}
//----------------------------------------------------------------------
void
MPIInfo::ShowGrid(void) {
  for (int iz = 0; iz < g[Z]; iz++) {
    for (int iy = 0; iy < g[Y]; iy++) {
      mout << "# ";
      for (int ix = 0; ix < g[X]; ix++) {
        int index = ix + iy * g[X] + iz * g[X] * g[Y];
        mout << std::setw(3);
        mout << grid[index];
      }
      mout << std::endl;
    }
    mout << "# " << std::endl;
  }
}
//----------------------------------------------------------------------
void
MPIInfo::RegisterPosition(void) {
  for (int i = 0; i < num_procs; i++) {
    int index = grid[i];
    int x, y, z;
    Index2Pos(i, x, y, z);
    qx[index] = x;
    qy[index] = y;
    qz[index] = z;
  }
}
//----------------------------------------------------------------------
int
MPIInfo::GetNeighborID(int ix, int iy, int iz, int dir) {
  int dx = dir % 3 - 1;
  dir /= 3;
  int dy = dir % 3 - 1;
  dir /= 3;
  int dz = dir % 3 - 1;
  ix += dx;
  iy += dy;
  iz += dz;
  if (ix < 0) ix = g[X] - 1;
  else if (ix >= g[X]) ix = 0;
  if (iy < 0) iy = g[Y] - 1;
  else if (iy >= g[Y]) iy = 0;
  if (iz < 0) iz = g[Z] - 1;
  else if (iz >= g[Z]) iz = 0;
  int id = ix + iy * g[X] + iz * g[X] * g[Y];
  return grid[id];
}
//----------------------------------------------------------------------
bool
MPIInfo::IsOverBoundary(int ix, int iy, int iz, int dir) {
  int dx = dir % 3 - 1;
  dir /= 3;
  int dy = dir % 3 - 1;
  dir /= 3;
  int dz = dir % 3 - 1;
  ix += dx;
  iy += dy;
  iz += dz;
  if (ix < 0) return true;
  else if (ix >= g[X]) return true;
  if (iy < 0) return true;
  else if (iy >= g[Y]) return true;
  if (iz < 0) return true;
  else if (iz >= g[Z]) return true;

  return false;
}
//----------------------------------------------------------------------
int
MPIInfo::GetOppositeNeighborID(int ix, int iy, int iz, int dir) {
  int dx = dir % 3 - 1;
  dir /= 3;
  int dy = dir % 3 - 1;
  dir /= 3;
  int dz = dir % 3 - 1;
  ix -= dx;
  iy -= dy;
  iz -= dz;
  if (ix < 0) ix = g[X] - 1;
  else if (ix >= g[X]) ix = 0;
  if (iy < 0) iy = g[Y] - 1;
  else if (iy >= g[Y]) iy = 0;
  if (iz < 0) iz = g[Z] - 1;
  else if (iz >= g[Z]) iz = 0;
  int id = ix + iy * g[X] + iz * g[X] * g[Y];
  return grid[id];
}
//----------------------------------------------------------------------
void
MPIInfo::GetGridPosition(int q[D]) {
  q[X] = qx[rank];
  q[Y] = qy[rank];
  q[Z] = qz[rank];
}
//----------------------------------------------------------------------
int
MPIInfo::GetParity(int dir) {
  switch (dir) {
  case D_LEFT:
  case D_RIGHT:
    return (qx[rank] % 2);

  case D_FORWARD:
  case D_BACK:
    return (qy[rank] % 2);

  case D_UP:
  case D_DOWN:
    return (qz[rank] % 2);
  }
  show_error("At GetParity");
  return -1;
}
//----------------------------------------------------------------------
bool
MPIInfo::IsEdge(int dir) {
  int dx = dir % 3 - 1;
  dir /= 3;
  int dy = dir % 3 - 1;
  dir /= 3;
  int dz = dir % 3 - 1;
  dx += qx[rank];
  dy += qy[rank];
  dz += qz[rank];
  if (dx < 0) return true;
  if (dy < 0) return true;
  if (dz < 0) return true;
  if (dx >= g[X]) return true;
  if (dy >= g[Y]) return true;
  if (dz >= g[Z]) return true;
  return false;
}
//----------------------------------------------------------------------

