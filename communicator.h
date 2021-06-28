//----------------------------------------------------------------------
// MPI Communication Class
//----------------------------------------------------------------------

#ifndef communicator_h
#define communicator_h
#include <mpi.h>
#include "mdconfig.h"

class Communicator {
public:
  static void Barrier(void);
  static void SendInteger(int &number, int dest_rank);
  static void RecvInteger(int &number, int src_rank);
  static void SendDouble(double *buffer, int number, int dest_rank);
  static void RecvDouble(double *buffer, int number, int src_rank);
  static void SendInteger(int *buffer, int number, int dest_rank);
  static void RecvInteger(int *buffer, int number, int src_rank);
  static void SendRecvInteger(int &sendnumber, int dest_rank, int &recvnumber, int src_rank);
  static void SendRecvDouble(void *sendbuf, int sendnumber, int dest_rank, void *recvbuf, int recvnumber, int src_rank);

  static void ISendInteger(int &number, int rank, int tag, MPI_Request &req);
  static void IRecvInteger(int &number, int rank, int tag, MPI_Request &req);
  static void SendDouble(void *buf, int size, int rank);
  static void RecvDouble(void *buf, int size, int rank);
  static void Wait(MPI_Request &req);
  static double GetTime(void);

  static double FindMaxDouble(double value);
  static bool AllReduceBoolean(bool flag);
  static double AllReduceDouble(double value);
  static void AllReduceDoubleBuffer(double *sendbuf, int size, double *recvbuf);
  static int AllReduceInteger(int value);
  static unsigned long int AllReduceUnsignedLongInteger(unsigned long int value);
  static void AllGatherInteger(int *sendbuf, int number, int *recvbuf);
  static void AllGatherDouble(double *sendbuf, int number, double *recvbuf);
  static void BroadcastInteger(int &value, int root);

};
#endif
