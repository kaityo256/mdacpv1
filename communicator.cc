//----------------------------------------------------------------------
// MPI Communication Class
//----------------------------------------------------------------------
#include "communicator.h"
//----------------------------------------------------------------------
void
Communicator::Barrier(void) {
  MPI_Barrier(MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::SendInteger(int &number, int dest_rank) {
  MPI_Send(&number, 1, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::RecvInteger(int &number, int src_rank) {
  MPI_Status st;
  MPI_Recv(&number, 1, MPI_INT, src_rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
void
Communicator::SendRecvInteger(int &send_number, int dest_rank, int &recv_number, int src_rank) {
  MPI_Status st;
  MPI_Sendrecv(&send_number, 1, MPI_INT, dest_rank, 0, &recv_number, 1, MPI_INT, src_rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
void
Communicator::SendRecvDouble(void *sendbuf, int send_number, int dest_rank,
                             void *recvbuf, int recv_number, int src_rank) {
  MPI_Status st;
  MPI_Sendrecv(sendbuf, send_number, MPI_DOUBLE, dest_rank, 0, recvbuf, recv_number, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
void
Communicator::SendDouble(double *buffer, int n, int dest_rank) {
  MPI_Send(buffer, n, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::RecvDouble(double *buffer, int n, int src_rank) {
  MPI_Status st;
  MPI_Recv(buffer, n, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
void
Communicator::SendInteger(int *buffer, int n, int dest_rank) {
  MPI_Send(buffer, n, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::RecvInteger(int *buffer, int n, int src_rank) {
  MPI_Status st;
  MPI_Recv(buffer, n, MPI_INT, src_rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
void
Communicator::ISendInteger(int &number, int dest_rank, int tag, MPI_Request &req) {
  MPI_Isend(&number, 1, MPI_INT, dest_rank, tag, MPI_COMM_WORLD, &req);
}
//----------------------------------------------------------------------
void
Communicator::IRecvInteger(int &number, int src_rank, int tag, MPI_Request &req) {
  MPI_Irecv(&number, 1, MPI_INT, src_rank, tag, MPI_COMM_WORLD, &req);
}
//----------------------------------------------------------------------
void
Communicator::SendDouble(void *buf, int size, int rank) {
  MPI_Send(buf, size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::RecvDouble(void *buf, int size, int rank) {
  MPI_Status st;
  MPI_Recv(buf, size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
void
Communicator::Wait(MPI_Request &req) {
  MPI_Status stat;
  MPI_Wait(&req, &stat);
}
//----------------------------------------------------------------------
bool
Communicator::AllReduceBoolean(bool flag) {
  int n = 0;
  if (flag) {
    n = 1;
  }
  int sum = 0;
  MPI_Allreduce(&n, & sum, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  if (sum > 0) {
    return true;
  } else {
    return false;
  }
}
//----------------------------------------------------------------------
int
Communicator::AllReduceInteger(int value) {
  int sum  = 0;
  MPI_Allreduce(&value, & sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return sum;
}
//----------------------------------------------------------------------
unsigned long int
Communicator::AllReduceUnsignedLongInteger(unsigned long int value) {
  unsigned long int  sum  = 0;
  MPI_Allreduce(&value, & sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  return sum;
}
//----------------------------------------------------------------------
double
Communicator::FindMaxDouble(double value) {
  double max = 0;
  MPI_Allreduce(&value, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return max;
}
//----------------------------------------------------------------------
double
Communicator::AllReduceDouble(double value) {
  double sum = 0;
  MPI_Allreduce(&value, & sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sum;
}
//----------------------------------------------------------------------
void
Communicator::AllReduceDoubleBuffer(double *sendbuf, int size, double *recvbuf) {
  MPI_Allreduce(sendbuf, recvbuf, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::AllGatherInteger(int *sendbuf, int number, int *recvbuf) {
  MPI_Allgather(sendbuf, number, MPI_INT, recvbuf, number, MPI_INT, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::AllGatherDouble(double *sendbuf, int number, double *recvbuf) {
  MPI_Allgather(sendbuf, number, MPI_DOUBLE, recvbuf, number, MPI_DOUBLE, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
double
Communicator::GetTime(void) {
  return MPI_Wtime();
}
//----------------------------------------------------------------------
void
Communicator::BroadcastInteger(int &value, int root) {
  MPI_Bcast(&value, 1, MPI_INT, root, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
