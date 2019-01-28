
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

#include "elliptic.h"

template < int p_Nq >
dfloat ellipticSerialUpdate1NBPCGKernel(const hlong Nelements,
					const dfloat * __restrict__ cpu_invDegree,
					const dfloat * __restrict__ cpu_z,
					const dfloat * __restrict__ cpu_Z,
					const dfloat beta,
					dfloat * __restrict__ cpu_p,
					dfloat * __restrict__ cpu_s){
  

#define p_Np (p_Nq*p_Nq*p_Nq)

  cpu_p  = (dfloat*)__builtin_assume_aligned(cpu_p,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_s  = (dfloat*)__builtin_assume_aligned(cpu_s,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_z  = (dfloat*)__builtin_assume_aligned(cpu_z,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_Z  = (dfloat*)__builtin_assume_aligned(cpu_Z,  USE_OCCA_MEM_BYTE_ALIGN) ;

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN) ;
  
  dfloat pdots = 0;
  
  for(hlong e=0;e<Nelements;++e){
    for(int i=0;i<p_Np;++i){
      const hlong n = e*p_Np+i;

      dfloat pn = cpu_p[n];
      dfloat zn = cpu_z[n];
      dfloat sn = cpu_s[n];
      dfloat Zn = cpu_Z[n];

      pn = zn + beta*pn;
      sn = Zn + beta*sn;

      dfloat invDeg = cpu_invDegree[n];
      
      pdots += pn*sn*invDeg;

      cpu_p[n] = pn;
      cpu_s[n] = sn;
    }
  }

#undef p_Np
  
  return pdots;
}
				     
dfloat ellipticSerialUpdate1NBPCG(const int Nq, const hlong Nelements,
			       occa::memory &o_invDegree, occa::memory &o_z, occa::memory &o_Z, const dfloat beta,
			       occa::memory &o_p, occa::memory &o_s){

  const dfloat * __restrict__ cpu_z  = (dfloat*)__builtin_assume_aligned(o_z.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_Z  = (dfloat*)__builtin_assume_aligned(o_Z.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  dfloat * __restrict__ cpu_p  = (dfloat*)__builtin_assume_aligned(o_p.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  dfloat * __restrict__ cpu_s  = (dfloat*)__builtin_assume_aligned(o_s.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  const dfloat * __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(o_invDegree.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  dfloat pdots = 0;
  
  switch(Nq){
  case  2: pdots = ellipticSerialUpdate1NBPCGKernel <  2 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break; 
  case  3: pdots = ellipticSerialUpdate1NBPCGKernel <  3 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case  4: pdots = ellipticSerialUpdate1NBPCGKernel <  4 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case  5: pdots = ellipticSerialUpdate1NBPCGKernel <  5 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case  6: pdots = ellipticSerialUpdate1NBPCGKernel <  6 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case  7: pdots = ellipticSerialUpdate1NBPCGKernel <  7 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case  8: pdots = ellipticSerialUpdate1NBPCGKernel <  8 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case  9: pdots = ellipticSerialUpdate1NBPCGKernel <  9 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case 10: pdots = ellipticSerialUpdate1NBPCGKernel < 10 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case 11: pdots = ellipticSerialUpdate1NBPCGKernel < 11 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  case 12: pdots = ellipticSerialUpdate1NBPCGKernel < 12 > (Nelements, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s); break;
  }

  return pdots;
}

void ellipticNonBlockingUpdate1NBPCG(elliptic_t *elliptic,
				     occa::memory &o_z, occa::memory &o_Z, const dfloat beta,
				     occa::memory &o_p, occa::memory &o_s,
				     dfloat *localpdots, dfloat *globalpdots, MPI_Request *request){
  
  const cgOptions_t cgOptions = elliptic->cgOptions;
  
  mesh_t *mesh = elliptic->mesh;

  dfloat pdots = 0;
  
  if(cgOptions.serial==1 && cgOptions.continuous==1){
    
    pdots = ellipticSerialUpdate1NBPCG(mesh->Nq, mesh->Nelements, 
				       elliptic->o_invDegree,
				       o_z, o_Z, beta, o_p, o_s);

  }
  else{
    if(!cgOptions.continuous){ // e.g. IPDG
      printf("EXITING - NBPCG NOT IMPLEMENTED FOR IPDG\n");
      MPI_Finalize();
      exit(-1);
    }else{
    
      // p <= z + beta*p
      // s <= Z + beta*s
      // dot(p,s)
      elliptic->update1NBPCGKernel(mesh->Nelements*mesh->Np, elliptic->NblocksUpdatePCG,
				   elliptic->o_invDegree, o_z, o_Z, beta, o_p, o_s, elliptic->o_tmppdots);
      
      elliptic->o_tmpNormr.copyTo(elliptic->tmppdots);
      
      *localpdots = 0;
      for(int n=0;n<elliptic->NblocksUpdatePCG;++n){
	*localpdots += elliptic->tmppdots[n];
      }
    }
  }
  
  *globalpdots = 0;
  if(cgOptions.enableReductions)      
    MPI_Iallreduce(localpdots, globalpdots, 1, MPI_DFLOAT, MPI_SUM, mesh->comm, request);
  else
    *globalpdots = 1;

}


template < int p_Nq >
void ellipticSerialUpdate2NBPCGKernel(const hlong Nelements,
				      const dfloat * __restrict__ cpu_invDegree,
				      const dfloat * __restrict__ cpu_s,
				      const dfloat * __restrict__ cpu_S,
				      const dfloat alpha,
				      dfloat * __restrict__ cpu_r,
				      dfloat * __restrict__ cpu_z,
				      dfloat *rdotz,
				      dfloat *zdotz){
  
  
#define p_Np (p_Nq*p_Nq*p_Nq)
  
  cpu_r  = (dfloat*)__builtin_assume_aligned(cpu_r,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_z  = (dfloat*)__builtin_assume_aligned(cpu_z,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_s  = (dfloat*)__builtin_assume_aligned(cpu_s,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_S  = (dfloat*)__builtin_assume_aligned(cpu_S,  USE_OCCA_MEM_BYTE_ALIGN) ;
  
  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN) ;
  
  *rdotz = 0;
  *zdotz = 0;
  
  for(hlong e=0;e<Nelements;++e){
    for(int i=0;i<p_Np;++i){
      const hlong n = e*p_Np+i;

      dfloat rn = cpu_r[n];
      dfloat zn = cpu_z[n];
      dfloat sn = cpu_s[n];
      dfloat Sn = cpu_S[n];

      dfloat invDeg = cpu_invDegree[n];
      
      rn = rn - alpha*sn;
      zn = zn - alpha*Sn;

      *rdotz += rn*zn*invDeg;
      *zdotz += zn*zn*invDeg;

      cpu_r[n] = rn;
      cpu_z[n] = zn;
    }
  }

#undef p_Np

}
				     
void ellipticSerialUpdate2NBPCG(const int Nq, const hlong Nelements,
				occa::memory &o_invDegree, occa::memory &o_s, occa::memory &o_S, const dfloat alpha,
				occa::memory &o_r, occa::memory &o_z,
				dfloat *rdotz, dfloat *zdotz){
  
  const dfloat * __restrict__ cpu_s  = (dfloat*)__builtin_assume_aligned(o_s.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_S  = (dfloat*)__builtin_assume_aligned(o_S.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  dfloat * __restrict__ cpu_r  = (dfloat*)__builtin_assume_aligned(o_r.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  dfloat * __restrict__ cpu_z  = (dfloat*)__builtin_assume_aligned(o_z.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  const dfloat * __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(o_invDegree.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  
  switch(Nq){
  case  2: ellipticSerialUpdate2NBPCGKernel <  2 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break; 
  case  3: ellipticSerialUpdate2NBPCGKernel <  3 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case  4: ellipticSerialUpdate2NBPCGKernel <  4 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case  5: ellipticSerialUpdate2NBPCGKernel <  5 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case  6: ellipticSerialUpdate2NBPCGKernel <  6 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case  7: ellipticSerialUpdate2NBPCGKernel <  7 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case  8: ellipticSerialUpdate2NBPCGKernel <  8 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case  9: ellipticSerialUpdate2NBPCGKernel <  9 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case 10: ellipticSerialUpdate2NBPCGKernel < 10 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case 11: ellipticSerialUpdate2NBPCGKernel < 11 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  case 12: ellipticSerialUpdate2NBPCGKernel < 12 > (Nelements, cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, rdotz, zdotz); break;
  }
}

void ellipticNonBlockingUpdate2NBPCG(elliptic_t *elliptic,
				     occa::memory &o_s, occa::memory &o_S, const dfloat alpha,
				     occa::memory &o_r, occa::memory &o_z,
				     dfloat *localdots, dfloat *globaldots, MPI_Request *request){
  
  const cgOptions_t cgOptions = elliptic->cgOptions;
  
  mesh_t *mesh = elliptic->mesh;

  localdots[0] = 0; // rdotz
  localdots[1] = 0; // zdotz
  
  if(cgOptions.serial==1 && cgOptions.continuous==1){
    
    ellipticSerialUpdate2NBPCG(mesh->Nq, mesh->Nelements, 
			       elliptic->o_invDegree,
			       o_s, o_S, alpha, o_r, o_z,
			       localdots, localdots+1);
    
  }
  else{
    if(!cgOptions.continuous){ // e.g. IPDG
      printf("EXITING - NBPCG NOT IMPLEMENTED FOR IPDG\n");
      MPI_Finalize();
      exit(-1);
    }else{
      // r <= r - alpha*s
      // z <= z - alpha*S
      // dot(r,z)
      // dot(z,z)
      elliptic->update2NBPCGKernel(mesh->Nelements*mesh->Np, elliptic->NblocksUpdatePCG,
				   elliptic->o_invDegree, o_s, o_S, alpha, o_r, o_z, elliptic->o_tmprdotz, elliptic->o_tmpzdotz);
      
      elliptic->o_tmprdotz.copyTo(elliptic->tmprdotz);
      elliptic->o_tmpzdotz.copyTo(elliptic->tmpzdotz);
      
      localdots[0] = 0;
      localdots[1] = 0;
      for(int n=0;n<elliptic->NblocksUpdatePCG;++n){
	localdots[0] += elliptic->tmprdotz[n];
	localdots[1] += elliptic->tmpzdotz[n];
      }
    }
  }
  
  globaldots[0] = 1;
  globaldots[1] = 1;
  if(cgOptions.enableReductions)      
    MPI_Iallreduce(localdots, globaldots, 2, MPI_DFLOAT, MPI_SUM, mesh->comm, request);
  
}

