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

  switch(testCase){

  default:
  case 1:{
    // printf("MoveMesh: BOX\n");
    const dlong nx = 1;
    const dlong ny = 1;
    const dlong nt = 1;
    const dlong Lx = 20;  // BOX DIMX
    const dlong Ly = 20;  // BOX DIMY
    const dfloat t0 = sqrt(200);
    const dfloat Ax = 0.95;
    const dfloat Ay = 0.95;

    const dfloat omega = 2 * PI * nt / t0;
    // const dfloat S     = sin(omega * T);
    const dfloat kx    = 2 * PI * nx / Lx;
    const dfloat ky    = 2 * PI * ny / Ly;

    // Explicit deformation
    explicitDeformationKernel(mesh.NnonPmlElements,
                              mesh.o_nonPmlElements,
                              T,
                              Ax,
                              Ay,
                              kx,
                              ky,
                              omega,
                              o_VX0,
                              o_rhsX,
                              o_Vx);

    break;
  } // end case 1 (BOX)

  case 2:{
    printf("MoveMesh: PLUNGING AIRFOIL\n");
    // Plunging airfoil
    const dfloat u_inf = 0.1;
    const dfloat xc = 0.5, yc = 0.0;;
    const dfloat r1 = 1.00, r2 = 4.7;

    // // Slow Plunge
    const dfloat H0 = 0.08;  // Plunging amplitude
    const dfloat Sr = 0.46;  // Strouhal number

    const dfloat omega = u_inf*Sr/H0;
    // const dfloat f  = 0.01;  // Plunging frequency
    // const dfloat omega = 2*PI*f;
    const dfloat ex = 0.0, ey = 1.0; // Plunging directions

    const dfloat inv_dr = 1.0 / (r2 - r1);
    const dfloat h = H0 * sin(omega*T);
    const dfloat dhdt = H0 * omega * cos(omega*T);

    // Explicit deformation
    explicitDeformationKernel(mesh.NnonPmlElements,
                              mesh.o_nonPmlElements,
                              T,
                              xc,
                              yc,
                              ex,
                              ey,
                              h,
                              dhdt,
                              r1,
                              r2,
                              inv_dr,
                              o_VX0,
                              o_rhsX,
                              o_Vx);

    break;

  } // end case 2 (PLUNGINGAIRFOIL)

  case 3:{
    // printf("MoveMesh: SOLVE MESH\n");
    // AA: This may only solve for nonPml elements
    dlong Ntotal = (meshN1.Nelements+meshN1.totalHaloPairs)*meshN1.Np*mdsNfields;

    // Create the solution and rhs vectors in every direction
    deviceMemory<dfloat> o_rhsVx = platform.reserve<dfloat>(Ntotal);
    deviceMemory<dfloat> o_vxL   = platform.reserve<dfloat>(Ntotal);
    deviceMemory<dfloat> o_rhsVy, o_vyL;

    // set solution vector to zero
    platform.linAlg().set(meshN1.Nelements*meshN1.Np*mdsNfields, (dfloat)0.0, o_rhsVx);
    platform.linAlg().set(meshN1.Nelements*meshN1.Np*mdsNfields, (dfloat)0.0, o_vxL);

    // Create gather arrays
    deviceMemory<dfloat> o_GrhsVx = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);
    deviceMemory<dfloat> o_Gvx    = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);
    deviceMemory<dfloat> o_GrhsVy, o_Gvy;

    if(mdsSolver.deform_laplace){
      o_rhsVy = platform.reserve<dfloat>(Ntotal);
      o_vyL   = platform.reserve<dfloat>(Ntotal);
      platform.linAlg().set(meshN1.Nelements*meshN1.Np*1, (dfloat)0.0, o_rhsVy);
      platform.linAlg().set(meshN1.Nelements*meshN1.Np*1, (dfloat)0.0, o_vyL);

      o_GrhsVy = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);
      o_Gvy    = platform.reserve<dfloat>(mdsSolver.Ndofs+mdsSolver.Nhalo);
    }

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
                 o_rhsVx,
                 o_rhsVy);

    int maxIter = 5000;
    int verbose = 0;


    // Gather - Solve - Scatter
    mdsSolver.ogsMasked.Gather(o_GrhsVx, o_rhsVx, mdsNfields, ogs::Add, ogs::Trans);
    mdsSolver.ogsMasked.Gather(o_Gvx, o_vxL, mdsNfields, ogs::Add, ogs::NoTrans);

    if(mdsSolver.deform_laplace){

      mdsSolver.ogsMasked.Gather(o_GrhsVy, o_rhsVy, 1, ogs::Add, ogs::Trans);
      mdsSolver.ogsMasked.Gather(o_Gvy, o_vyL, 1, ogs::Add, ogs::NoTrans);

      int Nitery;
      Niter  = mdsSolver.Solve(mdsLinearSolver, o_Gvx, o_GrhsVx, mdsTOL, maxIter, verbose);
      Nitery = mdsSolver.Solve(mdsLinearSolver, o_Gvy, o_GrhsVy, mdsTOL, maxIter, verbose);
      // printf("Total Number of Iterations in x: %d\n", Niter);
      // printf("Total Number of Iterations in y: %d\n", Nitery);

      mdsSolver.ogsMasked.Scatter(o_vxL, o_Gvx, 1, ogs::NoTrans);
      mdsSolver.ogsMasked.Scatter(o_vyL, o_Gvy, 1, ogs::NoTrans);
      o_GrhsVx.free(); o_Gvx.free();
      o_GrhsVy.free(); o_Gvy.free();

    } else if(mdsSolver.deform_linElastic){
      Niter = mdsSolver.Solve(mdsLinearSolver, o_Gvx, o_GrhsVx, mdsTOL, maxIter, verbose);
      // printf("Total Number of Iterations: %d\n", Niter);
      mdsSolver.ogsMasked.Scatter(o_vxL, o_Gvx, mdsNfields, ogs::NoTrans);
      o_GrhsVx.free(); o_Gvx.free();      
    }

    aleBCKernel(meshN1.Nelements,
                meshN1.o_x,
                meshN1.o_y,
                meshN1.o_z,
                T,
                mdsSolver.o_mapB,
                o_rhsX,
                o_vxL, 
                o_vyL); 

    break;

  } // end case 3 (Solve Mesh)

  case 4:{ 
    // printf("MoveMesh: 3D TGV\n");
    // Explicit deformation for 3D TGV
    const dfloat A  = 3.141592654359/6.;
    const dfloat TG = 20;
    explicitDeformationKernel(mesh.NnonPmlElements,
                              mesh.o_nonPmlElements,
                              T,
                              A,
                              TG,
                              o_VX0,
                              o_rhsX,
                              o_Vx);

    break;
  } // end case 4 (3D TGV)

  case 5: {
    // printf("MoveMesh: CARANGIFORM FISH\n");
    // Explicit deformation of the carangiform fish
    const dfloat r1 = 0.20, r2 = 0.80;
    const dfloat inv_dr = 1.0 / (r2 - r1);
    const dfloat U_inf = 0.1; // swimming speed

    const dfloat a0 = 0.02, a1 = -0.08, a2 = 0.16;
    const dfloat k = 2*PI; // Wave number
    const dfloat Af = 0.2; // Peak to peak tailbeat amplitude 0.2/L
    const dfloat St = 0.91; // Strouhal number
    const dfloat f = St * U_inf / Af; // tailbeat frequency
    const dfloat omega = 2.0*PI*f;

    // Explicit deformation
    explicitDeformationKernel(mesh.NnonPmlElements,
                              mesh.o_nonPmlElements,
                              T,
                              a0,
                              a1,
                              a2,
                              r1,
                              r2,
                              inv_dr,
                              omega,
                              k,
                              o_VX0,
                              o_rhsX,
                              o_Vx);
    break;
  } // end case 5 (carangiform fish)
}

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
                         o_meshVely,  
                         o_meshVelz);  
}
