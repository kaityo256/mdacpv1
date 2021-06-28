#include <iostream>
#include "mdrect.h"
//----------------------------------------------------------------------
int
MDRect::GetDirection(double q[D]) {
  int dir = 0;
  if (q[X] > e[X]) dir += 2;
  else if (q[X] >= s[X]) dir += 1;
  if (q[Y] > e[Y]) dir += 6;
  else if (q[Y] >= s[Y]) dir += 3;
  if (q[Z] > e[Z]) dir += 18;
  else if (q[Z] >= s[Z]) dir += 9;
  return dir;
}
//---------------------------------------------------------------------
bool
MDRect::IsInsideEdge(int dir, double q[D], SimulationInfo *sinfo) {
  const double SL = sinfo->SearchLength;
  switch (dir) {
  case D_LEFT:
    return (q[X] < s[X] + SL && q[X] > s[X]);

  case D_RIGHT:
    return (q[X] > e[X] - SL && q[X] < e[X]);

  case D_BACK:
    return (q[Y] < s[Y] + SL && q[Y] > s[Y]);

  case D_FORWARD:
    return (q[Y] > e[Y] - SL && q[Y] < e[Y]);

  case D_DOWN:
    return (q[Z] < s[Z] + SL && q[Z] > s[Z]);

  case D_UP:
    return (q[Z] > e[Z] - SL && q[Z] < e[Z]);

  }
  return false;
}
//---------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const MDRect& rect) {
  return (os << "(" << rect.s[X] << "," << rect.s[Y] << "," << rect.s[Z] << ")-("
          << rect.e[X] << "," << rect.e[Y] << "," << rect.e[Z] << ")");
}
//---------------------------------------------------------------------

