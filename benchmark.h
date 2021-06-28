//----------------------------------------------------------------------
#ifndef benchmark_h
#define benchmark_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class Benchmark : public Project {
private:
  void Thermalize(MDUnit *mdu, Parameter &param);
  void LoadDumpFileWithDuplication(const char *filename, MDUnit *mdu);
public:
  Benchmark(void) {
    ProjectManager::GetInstance().AddProject("Benchmark", this);
  };
  void Run(MDUnit *mdu, Parameter &param);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

