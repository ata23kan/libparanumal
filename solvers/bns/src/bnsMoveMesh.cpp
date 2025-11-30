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
#define PI 3.14159265

void bns_t::MoveMesh(deviceMemory<dfloat>& o_Vx, deviceMemory<dfloat>& o_rhsX, const dfloat T){

  // switch(testCase){

  // default:
  // case 1:{
  //   const dlong nx = 1;
  //   const dlong ny = 1;
  //   const dlong nt = 1;
  //   const dlong Lx = 20;  // BOX DIMX
  //   const dlong Ly = 20;  // BOX DIMY
  //   const dfloat t0 = sqrt(50);
  //   const dfloat Ax = 0.5;
  //   const dfloat Ay = 0.5;

  //   const dfloat omega = 2 * PI * nt / t0;
  //   // const dfloat S     = sin(omega * T);
  //   const dfloat kx    = 2 * PI * nx / Lx;
  //   const dfloat ky    = 2 * PI * ny / Ly;

  //   // Explicit deformation
  //   explicitDeformationKernel(mesh.NnonPmlElements,
  //                             mesh.o_nonPmlElements,
  //                             T,
  //                             Ax,
  //                             Ay,
  //                             kx,
  //                             ky,
  //                             omega,
  //                             o_VX0,
  //                             o_rhsX,
  //                             o_Vx);
  // } // end case 1 (BOX)

  // case 2:{
  //   // Plunging airfoil
  //   const dfloat xc = 0.5, yc = 0.0;;
  //   const dfloat r1 = 0.20, r2 = 4.0;
  //   const dfloat H0 = 0.2; // Plunging amplitude
  //   const dfloat f  = 0.01;  // Plunging frequency
  //   const dfloat omega = 2*PI*f;
  //   const dfloat ex = 0.0, ey = 1.0; // Plunging directions

  //   const dfloat inv_dr = 1.0 / (r2 - r1);
  //   const dfloat h = H0 * sin(omega*T);
  //   const dfloat dhdt = H0 * omega * cos(omega*T);

  //   // Explicit deformation
  //   explicitDeformationKernel(mesh.NnonPmlElements,
  //                             mesh.o_nonPmlElements,
  //                             T,
  //                             xc,
  //                             yc,
  //                             ex,
  //                             ey,
  //                             h,
  //                             dhdt,
  //                             r1,
  //                             r2,
  //                             inv_dr,
  //                             o_VX0,
  //                             o_rhsX,
  //                             o_Vx);


  // } // end case 2 (PLUNGINGAIRFOIL)

  // case 3:{

    // AA: This may only solve for nonPml elements
    dlong Ntotal = (meshN1.Nelements+meshN1.totalHaloPairs)*meshN1.Np*mdsNfields;

    deviceMemory<dfloat> o_rhsV = platform.reserve<dfloat>(Ntotal);
    deviceMemory<dfloat> o_vL   = platform.reserve<dfloat>(Ntotal);

    // set x to zero
    platform.linAlg().set(meshN1.Nelements*meshN1.Np*mdsNfields, (dfloat)0.0, o_vL);
    platform.linAlg().set(meshN1.Nelements*meshN1.Np*mdsNfields, (dfloat)0.0, o_rhsV);

    aleRhsKernel(mesh.Nelements,
                 meshN1.o_wJ,
                 meshN1.o_ggeo,
                 meshN1.o_vgeo,
                 meshN1.o_S,
                 meshN1.o_Se,
                 meshN1.o_vmapM,
                 mdsLambda,
                 mdsMu,
                 T,
                 meshN1.o_x,
                 meshN1.o_y,
                 meshN1.o_z,
                 mdsSolver.o_mapB,
                 o_rhsV);

    int maxIter = 5000;
    int verbose = 0;

    // Create gather arrays
    deviceMemory<dfloat> o_GrhsV = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);
    deviceMemory<dfloat> o_Gv    = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);

    // Gather - Solve - Scatter
    mdsSolver.ogsMasked.Gather(o_GrhsV, o_rhsV, mdsNfields, ogs::Add, ogs::Trans);
    mdsSolver.ogsMasked.Gather(o_Gv, o_vL, mdsNfields, ogs::Add, ogs::NoTrans);
    Niter = mdsSolver.Solve(mdsLinearSolver, o_Gv, o_GrhsV, mdsTOL, maxIter, verbose);
    // printf("Total Number of Iterations: %d\n", Niter);
    mdsSolver.ogsMasked.Scatter(o_vL, o_Gv, mdsNfields, ogs::NoTrans);
    o_GrhsV.free(); o_Gv.free();

    aleBCKernel(meshN1.Nelements,
                meshN1.o_x,
                meshN1.o_y,
                meshN1.o_z,
                T,
                mdsSolver.o_mapB,
                o_rhsX,
                o_vL); 

  // } // end case 3 (Solve Mesh)

  // }
  // Interpolate vertex velocities to the computational nodes
  velInterpolationKernel(mesh.NnonPmlElements,
                         mesh.o_nonPmlElements,
                         o_IM,
                         o_rhsX,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_meshVelx,
                         o_meshVely);  
}
