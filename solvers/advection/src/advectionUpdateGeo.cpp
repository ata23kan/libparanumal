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

void advection_t::UpdateGeo(const dfloat rk_dt){

	dlong Nverts = meshN1.Nelements*meshN1.Np;

	// Intermediate RK stage positions of vertices
	deviceMemory<dfloat> o_xrk = platform.reserve<dfloat>(Nverts);
	deviceMemory<dfloat> o_yrk = platform.reserve<dfloat>(Nverts);

	// x^rk = x^n + c_rk*dt*w^n+1
	platform.linAlg().zaxpy(Nverts, 1.0, meshN1.o_x, rk_dt, o_meshVelx, o_xrk);
	platform.linAlg().zaxpy(Nverts, 1.0, meshN1.o_y, rk_dt, o_meshVely, o_yrk);

  updateVgeoKernel(mesh.Nelements,
                   o_xrk,
                   o_yrk,
                   meshN1.o_z, // AA: change this 
                   meshN1.o_wJ,
                   meshN1.o_vgeo,
                   meshN1.o_ggeo,
                   mesh.o_wJ,
                   mesh.o_vgeo,
                   mesh.o_ggeo);

  updateSgeoKernel(mesh.Nelements,
                   o_xrk,
                   o_yrk,
                   meshN1.o_z,  // AA: change
                   meshN1.o_sgeo,
                   mesh.o_sgeo);

  o_xrk.free(); o_yrk.free();


} // end of function UpdateGeo