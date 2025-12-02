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

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

void mds_t::BuildOperatorMatrix(parAlmond::parCOO& A) {

  switch(mesh.elementType){
  case Mesh::TRIANGLES:
    BuildOperatorMatrixLaplacianTri2D(A); break;
  case Mesh::QUADRILATERALS:
  {
    LIBP_FORCE_ABORT("Quadrilaterals are not supported yet for mesh deformation!");
    // if(mesh.dim==2)
    //   BuildOperatorMatrixContinuousQuad2D(A);
    // else
    //   BuildOperatorMatrixContinuousQuad3D(A);

    // break;
  }
  case Mesh::TETRAHEDRA:
    LIBP_FORCE_ABORT("Tetrahedrons are not supported yet for mesh deformation!");
    // BuildOperatorMatrixContinuousTet3D(A); break;
  case Mesh::HEXAHEDRA:
    LIBP_FORCE_ABORT("Hexahedral elements are not supported yet for mesh deformation!");
    // BuildOperatorMatrixContinuousHex3D(A); break;
  }
}

void mds_t::BuildOperatorMatrixLaplacianTri2D(parAlmond::parCOO& A) {

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked.Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1, 0);
  A.globalColStarts.malloc(mesh.size+1, 0);
  mesh.comm.Allgather(Ngather, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  // Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements;

  memory<parAlmond::parCOO::nonZero_t> sendNonZeros(nnzLocal);
  memory<int> AsendCounts (mesh.size, 0);
  memory<int> ArecvCounts (mesh.size);
  memory<int> AsendOffsets(mesh.size+1);
  memory<int> ArecvOffsets(mesh.size+1);

  memory<dfloat> Srr = mesh.Srr;
  memory<dfloat> Srs = mesh.Srs;
  memory<dfloat> Sss = mesh.Sss;
  memory<dfloat> MM  = mesh.MM ;

  if(Comm::World().rank()==0) {printf("Building full FEM matrix...");fflush(stdout);}

  //Build unassembed non-zeros
  dlong cnt =0;
  for (dlong e=0;e<mesh.Nelements;e++) {
    dfloat Grr = mesh.ggeo[e*mesh.Nggeo + mesh.G00ID];
    dfloat Grs = mesh.ggeo[e*mesh.Nggeo + mesh.G01ID];
    dfloat Gss = mesh.ggeo[e*mesh.Nggeo + mesh.G11ID];
    // dfloat J   = mesh.wJ[e];

    for (int n=0;n<mesh.Np;n++) {
      if (maskedGlobalNumbering[e*mesh.Np + n]<0) continue; //skip masked nodes
      for (int m=0;m<mesh.Np;m++) {
        if (maskedGlobalNumbering[e*mesh.Np + m]<0) continue; //skip masked nodes

        dfloat val = 0.;

        val += Grr*Srr[m+n*mesh.Np];
        val += Grs*Srs[m+n*mesh.Np];
        val += Gss*Sss[m+n*mesh.Np];
        // val += J*lambda*MM[m+n*mesh.Np];

        dfloat nonZeroThreshold = 1e-7;
        if (fabs(val)>nonZeroThreshold) {
          // pack non-zero
          sendNonZeros[cnt].val = val;
          sendNonZeros[cnt].row = maskedGlobalNumbering[e*mesh.Np + n];
          sendNonZeros[cnt].col = maskedGlobalNumbering[e*mesh.Np + m];
          cnt++;
        }
      }
    }
  }

  // sort by row ordering
  sort(sendNonZeros.ptr(), sendNonZeros.ptr()+cnt,
       [](const parAlmond::parCOO::nonZero_t& a,
          const parAlmond::parCOO::nonZero_t& b) {
         if (a.row < b.row) return true;
         if (a.row > b.row) return false;

         return a.col < b.col;
        });

  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<cnt;++n) {
    const hlong id = sendNonZeros[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  mesh.comm.Alltoall(AsendCounts, ArecvCounts);

  // find send and recv offsets for gather
  A.nnz = 0;
  AsendOffsets[0] = 0;
  ArecvOffsets[0] = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  A.entries.malloc(A.nnz);

  // determine number to receive
  mesh.comm.Alltoallv(sendNonZeros, AsendCounts, AsendOffsets,
                      A.entries,    ArecvCounts, ArecvOffsets);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  sort(A.entries.ptr(), A.entries.ptr()+A.nnz,
       [](const parAlmond::parCOO::nonZero_t& a,
          const parAlmond::parCOO::nonZero_t& b) {
         if (a.row < b.row) return true;
         if (a.row > b.row) return false;

         return a.col < b.col;
       });

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<A.nnz;++n){
    if(A.entries[n].row == A.entries[cnt].row &&
       A.entries[n].col == A.entries[cnt].col){
       A.entries[cnt].val += A.entries[n].val;
    }
    else{
      ++cnt;
      A.entries[cnt] = A.entries[n];
    }
  }
  if (A.nnz) cnt++;
  A.nnz = cnt;

  if(Comm::World().rank()==0) printf("done.\n");
}