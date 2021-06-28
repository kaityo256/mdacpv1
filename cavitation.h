//----------------------------------------------------------------------
#ifndef cavitation_h
#define cavitation_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class Cavitation : public Project {
private:
  void MakeConfiguration(MDUnit *mdu, Parameter &param);
  void MakeConfigurationExtract(MDUnit *mdu, Parameter &param);
  void MakeBubble(MDUnit *mdu, Parameter &param);
  void Thermalize(MDUnit *mdu, Parameter &param);
public:
  void Run(MDUnit *mdu, Parameter &param);
  Cavitation(void) {
    ProjectManager::GetInstance().AddProject("Cavitation", this);
  };
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

