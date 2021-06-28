#include <iostream>
#include "mpiinfo.h"
#include "projectmanager.h"
//----------------------------------------------------------------------
int
main(int argc, char **argv) {
  setvbuf(stdout, NULL, _IOLBF, 0);
  MPIInfo minfo(&argc, &argv);
  std::string inputfile;

  if (argc > 1) {
    inputfile = argv[1];
  } else {
    mout << "# Input file is not specified. input.cfg is used." << std::endl;
    inputfile = "input.cfg";
  }

  Parameter param(inputfile.c_str());

  if (!param.IsValid()) {
    return -1;
  }

  if (!minfo.SetParameter(param)) {
    return -1;
  }

  ProjectManager::GetInstance().ExecuteProject(param, &minfo);
}
//----------------------------------------------------------------------

