//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "metastable.h"
#include "makeconf.h"
#include "bubblehist.h"
#include "communicator.h"
#include "mt.h"
//----------------------------------------------------------------------
MetaStable metastable;
//----------------------------------------------------------------------
void
MetaStable::Run(MDUnit *mdu, Parameter &param) {
  double density = param.GetDoubleDef("Density", 0.70);
  const double density_step = param.GetDoubleDef("DensityStep", 0.01);
  const double temperature = param.GetDoubleDef("AimedTemperature", 0.9);
  double *L = mdu->GetSystemSize();
  char output_filename[256];
  int step = 0;
  if (param.GetBooleanDef("ReadStep", true)) {
    std::ifstream ifs("step.dat");
    ifs >> step;
    ifs.close();
    Communicator::Barrier();
    if (0 == mdu->GetRank()) {
      std::ofstream ofs("step.dat");
      ofs << step + 1 << std::endl;
      ofs.close();
    }
    density = density - density_step * step;
    param.SetDouble("Density", density);
  }
  mout << "# Step = " << step << std::endl;
  sprintf(output_filename, "S%03d.dat", step);
  mdu->SetPeriodic(true);
  MakeConfigurationExtract(mdu, param);
  mdu->SetInitialVelocity(param.GetDoubleDef("InitialVelocity", 1.0));
  mdu->SetControlTemperature(true);
  mdu->ShowSystemInfo();
  mout << "#time,temperature,pressure,bubblesize" << std::endl;

  mdu->StartMeasurement();

  BubbleHist bhist(mdu->GetRank(), mdu->GetProcessNumber(), mdu->GetRect(), mdu->GetSystemSize());
  bhist.Clear();

  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  const int THERMALIZE_LOOP = param.GetIntegerDef("ThermalizeLoop", 10000);
  const int TOTAL_LOOP = param.GetIntegerDef("TotalLoop", 20000);
  const double MAX_CLUSTER_SIZE = param.GetDoubleDef("MaxClusterSize", 500.0);

  double pressure_average = 0.0;
  int num_observed = 0;
  for (int i = 0; i < TOTAL_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mdu->MakePairList();
      bhist.Analyse(mdu->GetVariables());
      double bubble_size = bhist.GetLargestClusterSize();
      double p = mdu->Pressure();
      mout << mdu->GetSimulationTime();
      mout << " " << mdu->Temperature();
      mout << " " << p;
      mout << " " << bubble_size;
      mout << std::endl;
      if (bubble_size > MAX_CLUSTER_SIZE) {
        mout << "# Bubble appears at " << mdu->GetSimulationTime() << std::endl;
        break;
      }
      if ( i > THERMALIZE_LOOP) {
        pressure_average += p;
        num_observed++;
      }
    }
  }

  pressure_average /= static_cast<double>(num_observed);
  mout << "# " << mdu->Density() << " " << pressure_average;
  mout << " " << temperature << " # (density,pressure,temperature) " << std::endl;
  mout.SaveToFile(output_filename);
  mdu->StopMeasurement();
  mdu->ShowMUPS(TOTAL_LOOP);
}
//----------------------------------------------------------------------
