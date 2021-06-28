//----------------------------------------------------------------------
#ifndef mpiinfo_h
#define mpiinfo_h
#include "parameter.h"
#include "mdconfig.h"
//----------------------------------------------------------------------
class MPIInfo {
private:
  int rank;
  int num_procs;
  int g[D];
  int *grid;
  int *qx, *qy, *qz;
  void MakeGrid(void);
  void Index2Pos(int index, int &ix, int &iy, int &iz);
  void RegisterPosition(void);
  bool SetGlobalGridSize(Parameter &param);
  bool SetLocalGridSize(Parameter &param);
public:
  MPIInfo(int *argc, char  ***argv);
  ~MPIInfo(void);
  void ShowGrid(void);
  bool SetParameter(Parameter &param);
  int GetRank(void) {return rank;};
  int GetProcessNumber(void) {return num_procs;};
  bool IsOverBoundary(int dir) {return IsOverBoundary(qx[rank], qy[rank], qz[rank], dir);};
  bool IsOverBoundary(int ix, int iy, int iz, int dir);
  int GetNeighborID(int dir) {return GetNeighborID(qx[rank], qy[rank], qz[rank], dir);};
  int GetNeighborID(int ix, int iy, int iz, int dir);
  int GetOppositeNeighborID(int dir) {return GetOppositeNeighborID(qx[rank], qy[rank], qz[rank], dir);};
  int GetOppositeNeighborID(int ix, int iy, int iz, int dir);
  void GetGridPosition(int q[D]);
  void GetGridSize(int g2[D]) {g2[X] = g[X]; g2[Y] = g[Y]; g2[Z] = g[Z];};
  int GetGridSize(int d) {return g[d];};
  int GetParity(int dir);
  bool IsEdge(int dir);
};
//----------------------------------------------------------------------
#endif
