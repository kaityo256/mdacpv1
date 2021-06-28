#include <iostream>
#include <sys/time.h>
#include "measure.h"
#include "mpistream.h"
//----------------------------------------------------------------------
Measure::Measure(void) {
  measure_times = 0;
  start_time = 0.0;
  sum_time = 0.0;
  measure_times = 0;
}
//----------------------------------------------------------------------
void
Measure::Clear(void) {
  measure_times = 0;
  sum_time = 0.0;
}
//----------------------------------------------------------------------
void
Measure::Start(void) {
  timeval tv;
  gettimeofday(&tv, NULL);
  start_time = tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

//----------------------------------------------------------------------
int
Measure::Stop(void) {
  timeval tv;
  gettimeofday(&tv, NULL);
  double stop_time = tv.tv_sec + (double)tv.tv_usec * 1e-6;
  sum_time += (stop_time - start_time);
  measure_times++;
  return measure_times;
}
//----------------------------------------------------------------------
double
Measure::Average(void) {
  double ave =  sum_time / (double)(measure_times);
  Clear();
  return ave;
}
//----------------------------------------------------------------------
void
Measure::ShowMUPS(const unsigned long int pn, const int LOOP) {
  double sec = Average();
  //double sec = sum_time;
  double mups = static_cast<double>(LOOP);
  mups = mups * static_cast<double>(pn) / sec / 1000000.0;
  mout << "# N=" << pn << " ";
  mout << sec << " [SEC] ";
  mout << mups << " [MUPS]" << std::endl;
}
//----------------------------------------------------------------------
