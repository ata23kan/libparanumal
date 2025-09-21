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
	const dlong nx = 1;
  const dlong ny = 1;
  const dlong nt = 1;
  const dlong Lx = 20;  // BOX DIMX
  const dlong Ly = 20;  // BOX DIMY
  const dfloat t0 = sqrt(50);
  const dfloat Ax = 0.5;
  const dfloat Ay = 0.5;

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
