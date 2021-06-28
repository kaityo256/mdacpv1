//----------------------------------------------------------------------
#ifndef chemipote_h
#define chemipote_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class ChemiPote : public Project {
private:
  void Observe(MDUnit *mdu, Parameter &param);
  void Thermalize(MDUnit *mdu, Parameter &param);
  void MakeConfHalf(MDUnit *mdu, Parameter &param);
  void MakeConfHalfFCC(MDUnit *mdu, Parameter &param);
  double Gamma(MDUnit *mdu);
public:
  void Run(MDUnit *mdu, Parameter &param);
  ChemiPote(void) {
    ProjectManager::GetInstance().AddProject("ChemiPote", this);
  };
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

