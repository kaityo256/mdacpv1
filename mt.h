//---------------------------------------------------------------------------
#ifndef mt_h
#define mt_h
//---------------------------------------------------------------------------
class MT {
public:
  static void SetSeed(int seed);
  static double GetDouble(void);
  static double GetGauss(void);
  static void Save(char * filename);
  static void Load(char * filename);
};
//---------------------------------------------------------------------------
#endif
