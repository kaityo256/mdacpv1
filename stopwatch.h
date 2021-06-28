#ifndef stopwatch_h
//----------------------------------------------------------------------
#include <vector>
//----------------------------------------------------------------------
class StopWatch {
private:
  std::vector<double> data;
  std::vector<double> stime;
  int count;
  double current_time;
  const char *basename;
  int id;
public:
  StopWatch(int rank, const char* bname) {
    basename = bname;
    count = 0;
    id = rank;
  };
  ~StopWatch(void) {
    //SaveToFile();
  }
  void Start(double st = 0.0) {
    stime.push_back(st);
    current_time = Communicator::GetTime();
  };
  void Stop(void) {
    data.push_back(Communicator::GetTime() - current_time);
  };

  void SaveToFile(void) {
    char filename[256];
    sprintf(filename, "%s%04d.dat", basename, id);
    std::ofstream ofs(filename);
    for (int i = 0; i < data.size(); i++) {
      ofs << i << " " << data[i] << " " << stime[i] << std::endl;
    }
  };
};
//----------------------------------------------------------------------
#endif
