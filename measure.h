//----------------------------------------------------------------------
#ifndef measure_h
#define measure_h
class Measure {
private:
  double start_time;
  double sum_time;
  int measure_times;
  void Clear(void);
public:
  Measure(void);
  void Start(void);
  int Stop(void);
  double Average(void);
  void ShowMUPS(const unsigned long int pn, const int LOOP);
};
#endif
