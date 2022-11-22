/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "esdg.hpp"

//evaluate ODE rhs = f(q,t)
void esdg_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  switch(mesh.elementType){
  case Mesh::TRIANGLES:
    rhsfTri2D(o_Q, o_RHS, T);
    break;
  case Mesh::QUADRILATERALS:
    rhsfQuad2D(o_Q, o_RHS, T);
    break;
  default:
    break;
  }
  
}

void esdg_t::rhsfTri2D(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){


    // entropy stable
  esInterpolateKernel(mesh.Nelements,
		      gamma,
		      o_EToE,
		      mesh.o_EToB,
		      o_esIqfT,
		      o_esFqT,
		      mesh.o_vgeo,
		      o_Q,
		      o_esQc,
		      o_esQe);

  // fix elements with entropy violations using Hadamard product
  esVolumeKernel(mesh.Nelements,
		 gamma,
		 mesh.o_vgeo,
		 o_esQNrT,
		 o_esQNsT,
		 o_esPqLfT,
		 o_esQe,
		 o_RHS);
  
  // surface terms
  esSurfaceKernel(mesh.Nelements,
		  T,
		  gamma,
		  mesh.o_sgeo,
		  o_esVmapM,
		  o_esVmapP,
		  o_EToB,
		  o_esLfT,
		  o_esX,
		  o_esY,
		  o_esZ,
		  o_esQc,
		  o_esQe,
		  o_esQp,
		  o_esQcrecon,
		  o_RHS);
  


}

void esdg_t::rhsfQuad2D(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

#if 0
  // fix elements with entropy violations using Hadamard product
  esVolumeKernel(mesh.Nelements,
		 mesh.o_esVgeo,
		 o_esVfgeo,
		 o_esFgeo,
		 o_esD1D,
		 o_esVf1D,
		 o_esLf1D,
		 o_esQ,
		 o_esQf,
		 o_RHS,
		 o_RHSf);

  // surface fluxes
  esSurfaceernel(mesh.Nelements,
		 mesh.o_esVgeo,
		 o_esFgeo,
		 o_esMapPq,
		 o_esBcFlag,
		 o_esLf1D,
		 o_esQf,
		 o_RHSf,
		 o_RHS);
  
#endif
}
