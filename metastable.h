//----------------------------------------------------------------------
#ifndef metastable_h
#define metastable_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class MetaStable : public Project {
private:
  void Thermalize(MDUnit *mdu, Parameter &param);
public:
  void Run(MDUnit *mdu, Parameter &param);
  MetaStable(void) {
    ProjectManager::GetInstance().AddProject("MetaStable", this);
  };
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

