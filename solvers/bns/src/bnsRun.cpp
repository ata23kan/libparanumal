/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "bns.hpp"

void bns_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         c,
                         nu,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_q);

  initialPositionKernel(mesh.Nelements,
                        meshN1.o_x,
                        meshN1.o_y,
                        meshN1.o_z,
                        o_VX);

  // Hold the initial positions for explicit deformation
  o_VX0.copyFrom(o_VX, meshN1.Np*meshN1.Nelements*mesh.dim, 0, properties_t("async", true));

  if (mesh.NpmlElements)
    pmlInitialConditionKernel(mesh.NpmlElements,
                             c,
                             nu,
                             startTime,
                             mesh.o_x,
                             mesh.o_y,
                             mesh.o_z,
                             o_pmlq);

  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat vmax = MaxWaveSpeed();

  dfloat dtAdv  = hmin/(vmax*(mesh.N+1.)*(mesh.N+1.));
  dfloat dtVisc = 1.0/tauInv;

  dfloat dt = (semiAnalytic) ? cfl*dtAdv : cfl*std::min(dtAdv, dtVisc);
  /*
    Artificial warping of time step size for multirate testing
    */
#if 0
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3"))
    dt /= (1<<(mesh.mrNlevels-1));
#endif
  timeStepper.SetTimeStep(dt);
  // printf("Initial time step: %e\n", timeStepper.GetTimeStep());

  // timeStepper.RunWithAle(*this, o_q, o_VX, startTime, finalTime);
  timeStepper.RunWithAlePml(*this, o_q, o_VX, o_pmlq, startTime, finalTime);

  printf("Total Volume Kernel time: %.4e\n", time_volume);
  printf("Total Pml Volume Kernel time: %.4e\n", time_volumePml);
  printf("Total Surface Kernel time: %.4e\n", time_surface);
  printf("Total Pml Surface Kernel time: %.4e\n", time_surfacePml);
  printf("Total Cubature Kernel time: %.4e\n", time_cubature);
  printf("Total Pml Cubature Kernel time: %.4e\n", time_cubaturePml);
  printf("Total MoveMesh Kernel time: %.4e\n", time_meshDeform);
  printf("Total Vgeo Kernel time: %.4e\n", time_vgeo);
  printf("Total Sgeo Kernel time: %.4e\n", time_sgeo);
  printf("Total UpdateX Kernel time: %.4e\n", time_updateX);

  // output norm of final solution
  {
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
    mesh.MassMatrixApply(o_q, o_Mq);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}


