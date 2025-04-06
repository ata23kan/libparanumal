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
  platform.linAlg().set(meshN1.Nelements*meshN1.Np*mdsNfields, (dfloat)0.0, o_rhs);

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

  aleBCKernel(meshN1.Nelements,
              meshN1.o_x,
              meshN1.o_y,
              meshN1.o_z,
              T,
              mdsSolver.o_mapB,
              o_dx,
              o_dy,
              o_xL); 

  platform.linAlg().axpy(meshN1.Nelements*meshN1.Np, 1.0, o_dx, 1.0, meshN1.o_x);  // update the vertex positions
  platform.linAlg().axpy(meshN1.Nelements*meshN1.Np, 1.0, o_dy, 1.0, meshN1.o_y);  // update the vertex positions

  platform.linAlg().set(meshN1.Nelements*meshN1.Np, (dfloat)0.0, o_meshVelx);  // make sure the mesh velocity is zero
  platform.linAlg().set(meshN1.Nelements*meshN1.Np, (dfloat)0.0, o_meshVely);  // make sure the mesh velocity is zero

  dfloat invaleT = 1/aleT;
  velInterpolationKernel(mesh.Nelements,
                      invaleT,
                      o_IM,
                      o_dx,
                      o_dy,
                      mesh.o_x,
                      mesh.o_y,
                      mesh.o_z,
                      o_meshVelx,
                      o_meshVely);


}