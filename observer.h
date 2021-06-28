//----------------------------------------------------------------------
#ifndef observer_h
#define observer_h
#include "mdconfig.h"
#include "meshlist.h"
#include "variables.h"
#include "mpiinfo.h"
//----------------------------------------------------------------------
class Observer {
private:
  double PotentialEnergyOnEdge(Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);
public:
  double Density(Variables *vars, SimulationInfo *sinfo);
  double Temperature(Variables *vars);
  double KineticEnergy(Variables *vars);
  double PotentialEnergy(Variables *vars, MeshList *mesh, SimulationInfo *sinfo, MPIInfo *minfo);
  double TotalEnergy(Variables *vars, MeshList *mesh, SimulationInfo *sinfo, MPIInfo *minfo);
  unsigned long int GetTotalParticleNumber(Variables *vars);
  double Pressure(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

