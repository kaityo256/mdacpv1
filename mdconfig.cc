#include "mdconfig.h"
//----------------------------------------------------------------------
const char *MDACP_VERSION = "MDACP Ver 1.02";
//----------------------------------------------------------------------
const char * Direction::name_str[MAX_DIR] = {"left-back-down", "back-down", "right-back-down", "left-down", "down", "right-down", "left-forward-down", "forward-down", "right-forward-down", "left-back", "back", "right-back", "left", "null", "right", "left-forward", "forward", "right-forward", "left-back-up", "back-up", "right-back-up", "left-up", "up", "right-up", "left-forward-up", "forward-up", "right-forward-up"};
//----------------------------------------------------------------------
int
Direction::Flip(int dir) {
  switch (dir) {
  case D_LEFT:
    return D_RIGHT;
  case D_RIGHT:
    return D_LEFT;
  case D_UP:
    return D_DOWN;
  case D_DOWN:
    return D_UP;
  case D_FORWARD:
    return D_BACK;
  case D_BACK:
    return D_FORWARD;
  }
  return D_NULL;
}
//----------------------------------------------------------------------
int
Direction::GetCoordinate(int dir) {
  switch (dir) {
  case D_LEFT:
  case D_RIGHT:
    return X;
  case D_FORWARD:
  case D_BACK:
    return Y;
  case D_UP:
  case D_DOWN:
    return Z;
  }
  return D_NULL;
}
//----------------------------------------------------------------------
