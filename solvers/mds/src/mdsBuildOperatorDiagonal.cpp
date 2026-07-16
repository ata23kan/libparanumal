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

#include "mds.hpp"

void mds_t::BuildOperatorDiagonal(memory<dfloat>& diagA){

  if(Comm::World().rank()==0) {printf("Building diagonal...");fflush(stdout);}

  memory<dfloat> diagAL(mesh.Np*mesh.Nelements);

  switch(mesh.elementType){
  case Mesh::TRIANGLES:
    BuildOperatorDiagonalLaplacianTri2D(diagAL); break;
  case Mesh::TETRAHEDRA:
    BuildOperatorDiagonalLaplacianTet3D(diagAL); break;
  case Mesh::QUADRILATERALS:
    LIBP_FORCE_ABORT("Quadrilaterals are not supported yet for mesh deformation!");
  case Mesh::HEXAHEDRA:
    LIBP_FORCE_ABORT("Hexahedral elements are not supported yet for mesh deformation!");
  }

  // gather the diagonal to assemble it
  ogsMasked.Gather(diagA, diagAL, 1, ogs::Add, ogs::Trans);

  if(Comm::World().rank()==0) printf("done.\n");
}

void mds_t::BuildOperatorDiagonalLaplacianTri2D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong gbase = eM*mesh.Nggeo;
    dfloat Grr = mesh.ggeo[gbase + mesh.G00ID];
    dfloat Grs = mesh.ggeo[gbase + mesh.G01ID];
    dfloat Gss = mesh.ggeo[gbase + mesh.G11ID];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
        A[eM*mesh.Np+n]  = Grr*mesh.Srr[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grs*mesh.Srs[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gss*mesh.Sss[n+n*mesh.Np];
      } else {
        A[eM*mesh.Np+n] = 1; //just put a 1 so A is invertable
      }
    }
  }
}

void mds_t::BuildOperatorDiagonalLaplacianTet3D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong gbase = eM*mesh.Nggeo;
    dfloat Grr = mesh.ggeo[gbase + mesh.G00ID];
    dfloat Grs = mesh.ggeo[gbase + mesh.G01ID];
    dfloat Grt = mesh.ggeo[gbase + mesh.G02ID];
    dfloat Gss = mesh.ggeo[gbase + mesh.G11ID];
    dfloat Gst = mesh.ggeo[gbase + mesh.G12ID];
    dfloat Gtt = mesh.ggeo[gbase + mesh.G22ID];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
        A[eM*mesh.Np+n]  = Grr*mesh.Srr[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grs*mesh.Srs[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grt*mesh.Srt[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gss*mesh.Sss[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gst*mesh.Sst[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gtt*mesh.Stt[n+n*mesh.Np];
      } else {
        A[eM*mesh.Np+n] = 1; //just put a 1 so A is invertable
      }
    }
  }
}
