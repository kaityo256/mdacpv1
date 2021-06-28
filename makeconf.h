//----------------------------------------------------------------------
#ifndef makeconf_h
#define makeconf_h
//----------------------------------------------------------------------
#include "mdconfig.h"
#include "variables.h"
#include "mdunit.h"
//----------------------------------------------------------------------
void MakeConfigurationFree(MDUnit *mdu, Parameter &param);
void MakeConfigurationPeriodic(MDUnit *mdu, Parameter &param);
void MakeConfigurationFCC(MDUnit *mdu, Parameter &param);
void MakeConfigurationFCC2(MDUnit *mdu, Parameter &param);
void MakeConfigurationExtract(MDUnit *mdu, Parameter &param);
#endif
//----------------------------------------------------------------------

