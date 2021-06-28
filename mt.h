#pragma once
#include <random>

class MT {
private:
  static std::mt19937 mt;

public:
  static void SetSeed(int seed);
  static double GetDouble(void);
  static double GetGauss(void);
};
