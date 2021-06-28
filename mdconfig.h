//---------------------------------------------------------------------
#ifndef mdconfig_h
#define mdconfig_h
//---------------------------------------------------------------------
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "parameter.h"
#include "mpistream.h"

const int D = 3;
const int X = 0, Y = 1, Z = 2;
const int N = 1000000;
const int PAIRLIST_SIZE = N * 50;

const double CUTOFF_LENGTH = 3.0;
//const double CUTOFF_LENGTH = 2.5;

extern const char *MDACP_VERSION;
//---------------------------------------------------------------------------
#define show_error(MSG) { mout << "# Error at " << __FILE__ <<":" << __LINE__ << std::endl;mout << MSG << std::endl;}
#define show_warning(MSG) { mout << "# Warning at " << __FILE__ <<":" << __LINE__ << std::endl;mout << MSG << std::endl;}
//---------------------------------------------------------------------------
const int MAX_DIR = 27;

enum DIRECTION {D_LEFT_BACK_DOWN, D_BACK_DOWN, D_RIGHT_BACK_DOWN, D_LEFT_DOWN, D_DOWN, D_RIGHT_DOWN, D_LEFT_FORWARD_DOWN, D_FORWARD_DOWN, D_RIGHT_FORWARD_DOWN, D_LEFT_BACK, D_BACK, D_RIGHT_BACK, D_LEFT, D_NULL, D_RIGHT, D_LEFT_FORWARD, D_FORWARD, D_RIGHT_FORWARD, D_LEFT_BACK_UP, D_BACK_UP, D_RIGHT_BACK_UP, D_LEFT_UP, D_UP, D_RIGHT_UP, D_LEFT_FORWARD_UP, D_FORWARD_UP, D_RIGHT_FORWARD_UP};
//---------------------------------------------------------------------------
enum HEATBATH_TYPE {HT_NOSEHOOVER, HT_LANGEVIN};

//---------------------------------------------------------------------------
class Direction {
private:
  static const char *name_str[MAX_DIR];
public:
  static const char * Name(int dir) {
    return name_str[dir];
  };
  static int Flip(int dir);
  static int GetCoordinate(int dir);
};
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------
