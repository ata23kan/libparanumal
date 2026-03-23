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

void bns_t::EnergyTGV(memory<dfloat>& Q, std::string fileName, dfloat time){

	FILE *fp;
	dfloat KE = 0.0;

	if(mesh.rank==0)
		fp = fopen(fileName.c_str(), "a+");

	memory<dfloat> r(mesh.Np);
	memory<dfloat> u(mesh.Np);
	memory<dfloat> v(mesh.Np);
	memory<dfloat> w(mesh.Np);

	mesh.o_wJ.copyTo(mesh.wJ);

	for(int e=0;e<mesh.Nelements;e++){
		for(int n=0;n<mesh.Np;n++){

			const int id = n + e*mesh.Np*Nfields;

			r[n] =   Q[id + 0*mesh.Np];
			u[n] = c*Q[id + 1*mesh.Np]/Q[id + 0*mesh.Np];
			v[n] = c*Q[id + 2*mesh.Np]/Q[id + 0*mesh.Np];
			w[n] = c*Q[id + 3*mesh.Np]/Q[id + 0*mesh.Np];
		}

		for(int n=0;n<mesh.Np;n++){
			dfloat Mu = 0.0, Mv = 0.0, Mw = 0.0;

			for(int k=0;k<mesh.Np;k++){
				Mu += mesh.MM[n+k*mesh.Np]*u[k]; 
				Mv += mesh.MM[n+k*mesh.Np]*v[k]; 
				Mw += mesh.MM[n+k*mesh.Np]*w[k]; 
			}
			KE += mesh.wJ[e]*r[n]*(u[n]*Mu + v[n]*Mv + w[n]*Mw);
		}
	}

	comm.Allreduce(KE, Comm::Sum);
	const dfloat volScale = 0.5/(6.28318530718*6.28318530718*6.28318530718); 

	KE *= volScale;

	if(mesh.rank==0){
		fprintf(fp, "%5.4f %13.12f\n", time, KE);
		fclose(fp);
	}
}