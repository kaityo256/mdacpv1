//----------------------------------------------------------------------
#ifndef fcalculator_h
#define fcalculator_h
//----------------------------------------------------------------------
#include "mdconfig.h"
#include "variables.h"
#include "meshlist.h"
#include "mpiinfo.h"
#include "simulationinfo.h"
//----------------------------------------------------------------------
class ForceCalculator {
private:
  static void CalculateForceUnroll(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void CalculateForceUnrollBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void CalculateForceSorted(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void CalculateForceSortedBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void CalculateForceNext(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void CalculateForceNextBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void CalculateForcePair(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void CalculateForceOnEdge(Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);

public:
  static void UpdatePositionHalf(Variables *vars, SimulationInfo *sinfo);
  static void CalculateForce(Variables *vars, MeshList *mesh, SimulationInfo *sinfo, MPIInfo *minfo);
  static void CalculateForceBenchmark(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
  static void HeatbathZeta(double &zeta, double ct, SimulationInfo *sinfo);
  static void HeatbathMomenta(Variables *vars, SimulationInfo *sinfo);
  static void Langevin(Variables *vars, SimulationInfo *sinfo);

};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
