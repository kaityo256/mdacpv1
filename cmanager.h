//----------------------------------------------------------------------
#ifndef cmanager_h
#define cmanager_h
//----------------------------------------------------------------------
#include "mdconfig.h"
#include "simulationinfo.h"
#include "variables.h"
#include "mpiinfo.h"
//----------------------------------------------------------------------
class CommunicationManager {
private:
  int send_number[MAX_DIR];
  int recv_number[MAX_DIR];
  int send_index[MAX_DIR][N];
  double *send_buffer[MAX_DIR];
  void MakeSendBuffer(int dir, double *s_buffer, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);
  void ReceiveParticle(int dir, double *r_buffer, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);
//  void SendNeighborPairwise(int dir, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);
  void MakeNeighborInformationSub(int dir, Variables *vars, SimulationInfo *sinfo, MDRect myrect, MPIInfo *minfo);
public:
  CommunicationManager(void);
  void SendParticles(Variables *vars, SimulationInfo *sinfo, MDRect myrect, MPIInfo *minfo);
  void SendParticlesSub(int dir, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);
  void MakeNeighborInformation(Variables *vars, SimulationInfo *sinfo, MDRect myrect, MPIInfo *minfo);

  void SendNeighborInformation(Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);
  void SendNeighborInformationSub(int dir, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
