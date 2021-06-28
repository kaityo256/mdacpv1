//----------------------------------------------------------------------
#include <iostream>
#include "cmanager.h"
#include "communicator.h"

//----------------------------------------------------------------------
CommunicationManager::CommunicationManager(void) {
#ifdef PAIR
  mout << "# Communication: Two-step oneway (PAIR)" << std::endl;
#elif PAIRWISE
  mout << "# Communication: Pairwise (PAIRWISE)" << std::endl;
#else
  mout << "# Communication: Simple Oneway (Oneway)" << std::endl;
#endif
}
//----------------------------------------------------------------------
void
CommunicationManager::SendParticles(Variables *vars, SimulationInfo *sinfo, MDRect myrect, MPIInfo *minfo) {
  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  static bool send_flag[N];

  for (int i = 0; i < pn; i++) {
    send_flag[i] = false;
  }

  for (int i = 0; i < MAX_DIR; i++) {
    send_number[i] = 0;
  }

  double x[D];
  for (int i = 0; i < pn; i++) {
    x[X] = q[i][X];
    x[Y] = q[i][Y];
    x[Z] = q[i][Z];
    int dir = myrect.GetDirection(x);
    int dest_rank = minfo->GetNeighborID(dir);
    if (D_NULL != dir && dest_rank != minfo->GetRank()) {
      send_index[dir][send_number[dir]] = i;
      send_number[dir]++;
      send_flag[i] = true;
    }
  }
  for (int i = 0; i < MAX_DIR; i++) {
    send_buffer[i] = new double[send_number[i]*D * 2];
  }

  for (int dir = 0; dir < MAX_DIR; dir++) {
    for (int j = 0; j < send_number[dir]; j++) {
      const int i = send_index[dir][j];
      send_buffer[dir][j * D * 2    ] = q[i][X];
      send_buffer[dir][j * D * 2 + 1] = p[i][X];
      send_buffer[dir][j * D * 2 + 2] = q[i][Y];
      send_buffer[dir][j * D * 2 + 3] = p[i][Y];
      send_buffer[dir][j * D * 2 + 4] = q[i][Z];
      send_buffer[dir][j * D * 2 + 5] = p[i][Z];
    }
  }

  //Pack
  int n = 0;
  for (int i = 0; i < pn; i++) {
    if (!send_flag[i]) {
      q[n][X] = q[i][X];
      q[n][Y] = q[i][Y];
      q[n][Z] = q[i][Z];
      p[n][X] = p[i][X];
      p[n][Y] = p[i][Y];
      p[n][Z] = p[i][Z];
      n++;
    }
  }
  vars->SetParticleNumber(n);

  //Send
  for (int dir = 0; dir < MAX_DIR; dir++) {
    if (D_NULL == dir)continue;
    SendParticlesSub(dir, vars, sinfo, minfo);
  }

  for (int i = 0; i < MAX_DIR; i++) {
    delete [] send_buffer[i];
  }
}
//----------------------------------------------------------------------
void
CommunicationManager::SendParticlesSub(int dir, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo) {
  int dest_rank = minfo->GetNeighborID(dir);
  int src_rank = minfo->GetOppositeNeighborID(dir);
  if (dest_rank == minfo->GetRank())return;
  if (src_rank == minfo->GetRank())return;

  Communicator::SendRecvInteger(send_number[dir], dest_rank, recv_number[dir], src_rank);

  double *recv_buffer = new double[recv_number[dir]*D * 2];
  Communicator::SendRecvDouble(send_buffer[dir], send_number[dir] * 2 * D, dest_rank,
                               recv_buffer, recv_number[dir] * 2 * D, src_rank);
  int n = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  for (int i = 0; i < recv_number[dir]; i++) {
    q[n][X] = recv_buffer[i * D * 2];
    p[n][X] = recv_buffer[i * D * 2 + 1];
    q[n][Y] = recv_buffer[i * D * 2 + 2];
    p[n][Y] = recv_buffer[i * D * 2 + 3];
    q[n][Z] = recv_buffer[i * D * 2 + 4];
    p[n][Z] = recv_buffer[i * D * 2 + 5];
    n++;
  }
  vars->SetParticleNumber(n);
  delete [] recv_buffer;
}
//----------------------------------------------------------------------
void
CommunicationManager::MakeNeighborInformation(Variables *vars, SimulationInfo *sinfo, MDRect myrect, MPIInfo *minfo) {

  vars->SetTotalParticleNumber(vars->GetParticleNumber());
  MakeNeighborInformationSub(D_LEFT, vars, sinfo, myrect, minfo);
  MakeNeighborInformationSub(D_RIGHT, vars, sinfo, myrect, minfo);
  MakeNeighborInformationSub(D_BACK, vars, sinfo, myrect, minfo);
  MakeNeighborInformationSub(D_FORWARD, vars, sinfo, myrect, minfo);
  MakeNeighborInformationSub(D_DOWN, vars, sinfo, myrect, minfo);
  MakeNeighborInformationSub(D_UP, vars, sinfo, myrect, minfo);
}
//----------------------------------------------------------------------
void
CommunicationManager::SendNeighborInformation(Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo) {
  vars->SetTotalParticleNumber(vars->GetParticleNumber());
  SendNeighborInformationSub(D_LEFT, vars, sinfo, minfo);
  SendNeighborInformationSub(D_RIGHT, vars, sinfo, minfo);
  SendNeighborInformationSub(D_BACK, vars, sinfo, minfo);
  SendNeighborInformationSub(D_FORWARD, vars, sinfo, minfo);
  SendNeighborInformationSub(D_DOWN, vars, sinfo, minfo);
  SendNeighborInformationSub(D_UP, vars, sinfo, minfo);
}
//----------------------------------------------------------------------
void
CommunicationManager::MakeNeighborInformationSub(int dir, Variables *vars, SimulationInfo *sinfo, MDRect myrect, MPIInfo *minfo) {

  int send_dir = dir;
  int recv_dir = dir;
#ifdef PAIR
  if (minfo->GetGridSize(Direction::GetCoordinate(dir)) > 1) {
    if (minfo->GetParity(dir)) {
      send_dir = dir;
      recv_dir = Direction::Flip(dir);
    } else {
      send_dir = Direction::Flip(dir);
      recv_dir = dir;
    }
  }
#endif

  double (*q)[D] = vars->q;
  const int pn = vars->GetTotalParticleNumber();
  int n = 0;
  double x[D];
  for (int i = 0; i < pn; i++) {
    x[X] = q[i][X];
    x[Y] = q[i][Y];
    x[Z] = q[i][Z];
    if (myrect.IsInsideEdge(send_dir, x, sinfo)) {
      send_index[send_dir][n] = i;
      n++;
    }
  }
  send_number[send_dir] = n;
  int dest_rank = minfo->GetNeighborID(send_dir);
  int src_rank = minfo->GetOppositeNeighborID(recv_dir);
  Communicator::SendRecvInteger(send_number[send_dir], dest_rank, recv_number[recv_dir], src_rank);
  SendNeighborInformationSub(dir, vars, sinfo, minfo);
}
//----------------------------------------------------------------------
void
CommunicationManager::MakeSendBuffer(int dir, double * s_buffer, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo) {
  double (*q)[D] = vars->q;
  for (int i = 0; i < send_number[dir]; i++) {
    int j = send_index[dir][i];
    s_buffer[i * D]   = q[j][X];
    s_buffer[i * D + 1] = q[j][Y];
    s_buffer[i * D + 2] = q[j][Z];
    if (minfo->IsOverBoundary(dir)) {
      if (D_LEFT == dir) s_buffer[i * D] += sinfo->L[X];
      if (D_RIGHT == dir) s_buffer[i * D] -= sinfo->L[X];
      if (D_BACK == dir) s_buffer[i * D + 1] += sinfo->L[Y];
      if (D_FORWARD == dir) s_buffer[i * D + 1] -= sinfo->L[Y];
      if (D_DOWN == dir) s_buffer[i * D + 2] += sinfo->L[Z];
      if (D_UP == dir) s_buffer[i * D + 2] -= sinfo->L[Z];
    }
  }
}
//----------------------------------------------------------------------
void
CommunicationManager::ReceiveParticle(int dir, double * r_buffer, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo) {
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;

  int n = vars->GetTotalParticleNumber();
  for (int i = 0; i < recv_number[dir]; i++) {
    q[n][X] = r_buffer[i * D];
    q[n][Y] = r_buffer[i * D + 1];
    q[n][Z] = r_buffer[i * D + 2];
    p[n][X] = 0.0;
    p[n][Y] = 0.0;
    p[n][Z] = 0.0;
    n++;
  }
  vars->SetTotalParticleNumber(n);
}
//----------------------------------------------------------------------
void
CommunicationManager::SendNeighborInformationSub(int dir, Variables *vars, SimulationInfo *sinfo, MPIInfo *minfo) {
  int send_dir = dir;
  int recv_dir = dir;
#ifdef PAIR
  if (minfo->GetGridSize(Direction::GetCoordinate(dir)) > 1) {
    if (minfo->GetParity(dir)) {
      send_dir = dir;
      recv_dir = Direction::Flip(dir);
    } else {
      send_dir = Direction::Flip(dir);
      recv_dir = dir;
    }
  }
#endif

  int dest_rank = minfo->GetNeighborID(send_dir);
  int src_rank = minfo->GetOppositeNeighborID(recv_dir);

  double *s_buffer = new double[send_number[send_dir]*D];
  double *r_buffer = new double[recv_number[recv_dir]*D];

  MakeSendBuffer(send_dir, s_buffer, vars, sinfo, minfo);

#ifdef PAIRWISE
  if (minfo->GetParity(dir)) {
    Communicator::SendDouble(s_buffer, send_number[dir]*D, dest_rank);
  } else {
    Communicator::RecvDouble(r_buffer, recv_number[dir]*D, src_rank);
  }
  Communicator::Barrier();
  if (minfo->GetParity(dir)) {
    Communicator::RecvDouble(r_buffer, recv_number[dir]*D, src_rank);
  } else {
    Communicator::SendDouble(s_buffer, send_number[dir]*D, dest_rank);
  }
#else
  Communicator::SendRecvDouble(s_buffer, send_number[send_dir]*D, dest_rank,
                               r_buffer, recv_number[recv_dir]*D, src_rank);
#endif

  ReceiveParticle(recv_dir, r_buffer, vars, sinfo, minfo);
  delete [] s_buffer;
  delete [] r_buffer;
}
//----------------------------------------------------------------------
