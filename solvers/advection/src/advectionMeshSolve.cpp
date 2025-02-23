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

#include "advection.hpp"

// Find the mesh velocity and update the geometric factors
void advection_t::MeshSolve(const dfloat T, const dfloat aleT){
// void advection_t::MeshSolve(const dfloat T, const dfloat dt){

	dlong Ntotal = (meshN1.Nelements+meshN1.totalHaloPairs)*meshN1.Np*mdsNfields;

	deviceMemory<dfloat> o_rhs = platform.reserve<dfloat>(Ntotal);
	deviceMemory<dfloat> o_xL = platform.reserve<dfloat>(Ntotal);

  // set x to zero
  platform.linAlg().set(meshN1.Nelements*meshN1.Np*mdsNfields, (dfloat)0.0, o_xL);

  aleRhsKernel(mesh.Nelements,
               meshN1.o_wJ,
               meshN1.o_ggeo,
               meshN1.o_vgeo,
               meshN1.o_S,
               meshN1.o_Se,
               meshN1.o_vmapM,
               lambda,
               mu,
               T,
               meshN1.o_x,
               meshN1.o_y,
               meshN1.o_z,
               mdsSolver.o_mapB,
               o_rhs);

  int maxIter = 5000;
  int verbose = 0;

  deviceMemory<dfloat> o_Grhs = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);
  deviceMemory<dfloat> o_Gx   = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);
  mdsSolver.ogsMasked.Gather(o_Grhs, o_rhs, mdsNfields, ogs::Add, ogs::Trans);
  mdsSolver.ogsMasked.Gather(o_Gx, o_xL, mdsNfields, ogs::Add, ogs::NoTrans);
  Niter = mdsSolver.Solve(mdsLinearSolver, o_Gx, o_Grhs, mdsTOL, maxIter, verbose);
  mdsSolver.ogsMasked.Scatter(o_xL, o_Gx, mdsNfields, ogs::NoTrans);
  o_Grhs.free(); o_Gx.free();

  // merge arrays back and enter BCs
  deviceMemory<dfloat> o_mQ = platform.reserve<dfloat>(mdsNfields*meshN1.Np*meshN1.Nelements);

  aleBCKernel(meshN1.Nelements,
              meshN1.o_x,
              meshN1.o_y,
              meshN1.o_z,
              T,
              mdsSolver.o_mapB,
              o_mQ,
              o_xL); 

  memory<dfloat> xuL(mesh.Np*mesh.Nelements);
  memory<dfloat> test(mesh.Np*mesh.Nelements);
  // o_mQ.copyTo(xuL);

  // for(int e=0; e<mesh.Nelements;e++){
  //   int id=e*meshN1.Np;
  //   int iid = e*meshN1.Np*mdsNfields;
  //   printf("Element: %d\n", e);
  //   printf("%f %f\t%f %f\n", meshN1.x[id+0], xuL[iid+0+0*meshN1.Np], meshN1.y[id+0], xuL[iid+0+1*meshN1.Np]);
  //   printf("%f %f\t%f %f\n", meshN1.x[id+1], xuL[iid+1+0*meshN1.Np], meshN1.y[id+1], xuL[iid+1+1*meshN1.Np]);
  //   printf("%f %f\t%f %f\n", meshN1.x[id+2], xuL[iid+2+0*meshN1.Np], meshN1.y[id+2], xuL[iid+2+1*meshN1.Np]);
  //   printf("\n");
  // }


  updateVgeoKernel(mesh.Nelements,
                   meshN1.o_x,
                   meshN1.o_y,
                   meshN1.o_z,
                   mesh.o_wJ,
                   mesh.o_vgeo,
                   mesh.o_ggeo);
                   // o_mQ

  updateSgeoKernel(mesh.Nelements,
                   meshN1.o_x,
                   meshN1.o_y,
                   meshN1.o_z,
                   mesh.o_sgeo);
                   // o_mQ

  deviceMemory<dfloat>o_test = platform.reserve<dfloat>(Ntotal);
  platform.linAlg().set(Ntotal, (dfloat)1.0, o_test);

  interpolationKernel(mesh.Nelements,
                      o_IM,
                      meshN1.o_x,
                      meshN1.o_y,
                      meshN1.o_z,
                      mesh.o_x,
                      mesh.o_y,
                      mesh.o_z,
                      o_meshVelx,
                      o_meshVely);

  // std::exit(EXIT_SUCCESS);

  // o_test.copyTo(test);

  // o_meshVely.copyTo(test);
  platform.linAlg().scale(mesh.Np*mesh.Nelements, 1/aleT, o_meshVelx);
  platform.linAlg().scale(mesh.Np*mesh.Nelements, 1/aleT, o_meshVely);
  // for(int id=0;id<Ntotal;id++){
  o_meshVely.copyTo(xuL);

  // for(int e=0; e<mesh.Nelements; e++){
  //   for(int n=0; n<mesh.Np; ++n){
  //     int id = e*mesh.Np + n;
  //     printf("%f\n", xuL[id]);
  //   }
  // }
  // std::exit(EXIT_SUCCESS);

  // printf("Time: %f\n", T);
  // printf("Ale Time: %f\n\n", aleT);
  // std::string name;
  // settings.getSetting("OUTPUT FILE NAME", name);
  // char fname[BUFSIZ];
  // sprintf(fname, "%s_newMesh_%04d_%04f.vtu", name.c_str(), meshN1.rank, T);
  // mdsSolver.PlotNewMesh2(xuL, fname);

}