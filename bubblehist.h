#ifndef bubblehist_h
#define bubblehist_h
//----------------------------------------------------------------------
#include <string.h>
#include "mdconfig.h"
#include "variables.h"
#include "mdunit.h"
//----------------------------------------------------------------------
class BubbleHist {
private:
  int num_observation;
  double grid_size;
  double *density_grid;
  bool *clustering_grid;
  int *cluster_number;
  int *cluster_size;
  int *index_list;
  int gnx, gny, gnz;
  int local_grid_x, local_grid_y, local_grid_z;
  int local_grid_num;
  int grid_num;
  double *size_distribution;
  double *L;
  int rank;
  int proc_num;
  MDRect myrect;
  void Clustering(void);
  int GetClusterIndex(int ix);
  int GetClusterIndex(int ix, int iy, int iz);
  int GetIndex(int ix, int iy, int iz);
  void Check(int, int, int, int, int, int);
  void MakeDensity(Variables *vars);
public:
  BubbleHist(int r, int pn, MDRect myrect, double *ssize);
  ~BubbleHist(void);
  void Clear();
  void Analyse(Variables *vars);
  double GetClusterSize(int index);
  double GetLargestClusterSize(void) {return GetClusterSize(0);};
  void SaveHistogram(std::string filename);
  void SaveBubbleDistribution(MDUnit *mdu, const char *base_dir = ".");
};
//----------------------------------------------------------------------
#endif

