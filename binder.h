//----------------------------------------------------------------------
#ifndef binder_h
#define binder_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class Binder : public Project {
private:
  void Observe(MDUnit *mdu, Parameter &param);
  void Thermalize(MDUnit *mdu, Parameter &param);
public:
  void Run(MDUnit *mdu, Parameter &param);
  Binder(void) {
    ProjectManager::GetInstance().AddProject("Binder", this);
  };
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

