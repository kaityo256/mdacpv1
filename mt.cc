#include "mt.h"
#include <algorithm>

std::mt19937 MT::mt;

void MT::SetSeed(int seed) {
  mt.seed(seed);
}

double MT::GetDouble(void) {
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  return ud(mt);
}

double MT::GetGauss(void) {
  std::normal_distribution<double> nd(0.0, 1.0);
  return nd(mt);
}