
#include <iostream>
#include "droplet.h"
#include "makeconf.h"
//----------------------------------------------------------------------
Droplet droplet;
//----------------------------------------------------------------------
void
Droplet::Run(MDUnit *mdu, Parameter &param) {
  MakeInitialCondition(mdu, param);
  const double v0 = param.GetDoubleDef("InitialVelocity", 1.0);
  mdu->SetInitialVelocity(v0);
  mdu->SetControlTemperature(true);
  mdu->ShowSystemInfo();
  const int OBSERVE_LOOP = param.GetIntegerDef("ObserveLoop", 100);
  const int TOTAL_LOOP = param.GetIntegerDef("TotalLoop", 1000);
  const bool save_cdview = param.GetBooleanDef("SaveCdviewFile", false);
  mdu->StartMeasurement();
  mout << "#Energy = " << mdu->TotalEnergy() << std::endl;
  for (int i = 0; i < TOTAL_LOOP; i++) {
    mdu->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdu->GetSimulationTime();
      mout << " " << mdu->Temperature();
      mout << "# observe" << std::endl;
      if (save_cdview) {
        mdu->SaveAsCdviewSequential();
      }
    }
  }
  mdu->StopMeasurement();
  mdu->ShowMUPS(TOTAL_LOOP);
  mout << "#Energy = " << mdu->TotalEnergy() << std::endl;
}
//----------------------------------------------------------------------
void
Droplet::MakeInitialCondition(MDUnit *mdu, Parameter &param) {
  double *L = mdu->GetSystemSize();
  double radius = param.GetDoubleDef("DropletRadius", 20);
  double density = param.GetDoubleDef("Density", 0.5);
  double c[D];
  c[X] = L[X] * 0.5;
  c[Y] = L[Y] * 0.5;
  c[Z] = L[Z] * 0.5;
  AddBall(c, radius, density, mdu);
}
//----------------------------------------------------------------------
void
Droplet::AddBall(double c[D], double R, double density, MDUnit *mdu) {
  double *L = mdu->GetSystemSize();
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L[X] / s);
  int sy = static_cast<int>(L[Y] / s);
  int sz = static_cast<int>(L[Z] / s);
  double x[D];

  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s;
        if ((x[X] - c[X]) * (x[X] - c[X]) + (x[Y] - c[Y]) * (x[Y] - c[Y]) + (x[Z] - c[Z]) * (x[Z] - c[Z]) < R * R) {
          mdu->AddParticle(x);
        }

        x[X] = static_cast<double>(ix) * s;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s + hs;
        if ((x[X] - c[X]) * (x[X] - c[X]) + (x[Y] - c[Y]) * (x[Y] - c[Y]) + (x[Z] - c[Z]) * (x[Z] - c[Z]) < R * R) {
          mdu->AddParticle(x);
        }

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s;
        x[Z] = static_cast<double>(iz) * s + hs;
        if ((x[X] - c[X]) * (x[X] - c[X]) + (x[Y] - c[Y]) * (x[Y] - c[Y]) + (x[Z] - c[Z]) * (x[Z] - c[Z]) < R * R) {
          mdu->AddParticle(x);
        }

        x[X] = static_cast<double>(ix) * s + hs;
        x[Y] = static_cast<double>(iy) * s + hs;
        x[Z] = static_cast<double>(iz) * s;
        if ((x[X] - c[X]) * (x[X] - c[X]) + (x[Y] - c[Y]) * (x[Y] - c[Y]) + (x[Z] - c[Z]) * (x[Z] - c[Z]) < R * R) {
          mdu->AddParticle(x);
        }
      }
    }
  }
}
//----------------------------------------------------------------------
