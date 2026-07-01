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

void bns_t::constantErrorNorm(memory<dfloat>& Q, dfloat rbar, dfloat Ubar, dfloat Vbar){

	dfloat rNorm = 0.f, uNorm = 0.f, vNorm = 0.f;

	memory<dfloat> r(mesh.Np);
	memory<dfloat> u(mesh.Np);
	memory<dfloat> v(mesh.Np);

	memory<dfloat> Mr(mesh.Np);
	memory<dfloat> Mu(mesh.Np);
	memory<dfloat> Mv(mesh.Np);

	for(dlong e=0;e<mesh.Nelements;e++){
		for(dlong n=0;n<mesh.Np;n++){

			const dlong id = n + e*mesh.Np*Nfields;

			r[n] =   Q[id + 0*mesh.Np] - rbar;
			u[n] = c*Q[id + 1*mesh.Np] - rbar*Ubar;
			v[n] = c*Q[id + 2*mesh.Np] - rbar*Vbar;
			// r[n] =   Q[id + 0*mesh.Np]    - rbar;
			// u[n] = c*Q[id + 1*mesh.Np]/rm - Ubar;
			// v[n] = c*Q[id + 2*mesh.Np]/rm - Vbar;


		}

		for(dlong n=0;n<mesh.Np;n++){
			dfloat r_Mr = 0.0, r_Mu = 0.0, r_Mv = 0.0;

			for(int k=0;k<mesh.Np;k++){
				r_Mr += mesh.MM[n+k*mesh.Np]*r[k]; 
				r_Mu += mesh.MM[n+k*mesh.Np]*u[k]; 
				r_Mv += mesh.MM[n+k*mesh.Np]*v[k]; 
			}
			Mr[n] = r_Mr;
			Mu[n] = r_Mu;
			Mv[n] = r_Mv;
		}

		for(dlong n=0; n<mesh.Np;n++){
				rNorm += r[n]*Mr[n]; 
				uNorm += u[n]*Mu[n]; 
				vNorm += v[n]*Mv[n]; 
		}
	}

	comm.Allreduce(rNorm, Comm::Sum);
	comm.Allreduce(uNorm, Comm::Sum);
	comm.Allreduce(vNorm, Comm::Sum);

	rNorm = sqrt(rNorm);
	uNorm = sqrt(uNorm);
	vNorm = sqrt(vNorm);

	if(mesh.rank==0){
		// printf("Density error norm: %e,  Ux error norm: %e, Uy error norm: %e\n", rNorm, uNorm, vNorm);
		printf("Density error norm: %e,  x-momentum error norm: %e, y-momentum error norm: %e\n", rNorm, uNorm, vNorm);
	}

}