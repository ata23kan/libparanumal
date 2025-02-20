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
void advection_t::MeshSolve(const dfloat T){
// void advection_t::MeshSolve(const dfloat T, const dfloat dt){

	dlong Ntotal = (meshN1.Nelements+meshN1.totalHaloPairs)*meshN1.Np*mdsNfields;

	deviceMemory<dfloat> o_rhs = platform.reserve<dfloat>(Ntotal);
	deviceMemory<dfloat> o_xL = platform.reserve<dfloat>(Ntotal);

  // set x to zero
  platform.linAlg().set(meshN1.Nelements*meshN1.Np*mdsNfields, (dfloat)0.0, o_xL);

  aleRhsKernel(mesh.Nelements,
               meshN1.o_wJ,
               meshN1.o_ggeo,
               meshN1.o_sgeo,
               meshN1.o_vgeo,
               meshN1.o_D,
               meshN1.o_S,
               meshN1.o_Se,
               meshN1.o_MM,
               meshN1.o_vmapM,
               meshN1.o_sM,
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
  memory<dfloat> xuL(Ntotal);

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

  // o_xL.copyTo(xuL);
  // int zerocnt = 0;
  // for(int i=0;i<Ntotal;i++){
  //   // printf("%f\n", xuL[i]);
  //   if (xuL[i]==0) zerocnt++;
  // }
  // printf("zerocnt: %d\n", zerocnt);

  updateGgeoKernel(mesh.Nelements,
                   mesh.o_wJ,
                   mesh.o_vgeo,
                   mesh.o_ggeo,
                   o_mQ);

  updateSgeoKernel(mesh.Nelements,
                   mesh.o_sgeo,
                   o_mQ);

  // std::exit(EXIT_SUCCESS);
  // o_mQ.copyTo(xuL);
  // std::string name;
  // settings.getSetting("OUTPUT FILE NAME", name);
  // char fname[BUFSIZ];
  // sprintf(fname, "%s_newMesh_%04d_%04f.vtu", name.c_str(), meshN1.rank, T);
  // mdsSolver.PlotNewMesh2(xuL, fname);

}