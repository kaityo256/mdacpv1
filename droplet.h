//----------------------------------------------------------------------
#ifndef collision_h
#define collision_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class Droplet : public Project {
private:
  void AddBall(double c[D], double r, double density, MDUnit *mdu);
  void MakeInitialCondition(MDUnit *mdu, Parameter &param);
public:
  Droplet(void) {
    ProjectManager::GetInstance().AddProject("Droplet", this);
  };
  void Run(MDUnit *mdu, Parameter &param);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

