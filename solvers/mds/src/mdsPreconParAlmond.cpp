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

#include "mdsPrecon.hpp"

//AMG preconditioner via parAlmond
void ParAlmondPrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  //hand off to parAlmond
  parAlmond.Operator(o_r, o_Mr);
}

ParAlmondPrecon::ParAlmondPrecon(mds_t& _mds):
  mds(_mds), settings(_mds.settings),
  parAlmond(mds.platform, settings, mds.mesh.comm) {

  //build full A matrix and pass to parAlmond
  if (Comm::World().rank()==0){
    printf("-----------------------------Multigrid AMG Setup--------------------------------------------\n");
  }
  parAlmond::parCOO A(mds.platform, mds.mesh.comm);
  // } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    mds.BuildOperatorMatrix(A);
  // }

  //populate null space unit vector
  int rank = mds.mesh.rank;
  int size = mds.mesh.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = static_cast<dlong>(A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);
  memory<pfloat> null(numLocalRows);
  for (dlong i=0;i<numLocalRows;i++) {
    null[i] = 1.0/sqrt(TotalRows);
  }

  parAlmond.AMGSetup(A, false, null, 0.);

  parAlmond.Report();
  
  //The csr matrix at the top level of parAlmond may have a larger
  // halo region than the matrix free kernel. Adjust if necessary
  dlong parAlmondNrows = parAlmond.getNumRows(0);
  dlong parAlmondNcols = parAlmond.getNumCols(0);
  dlong parAlmondNhalo = parAlmondNcols - parAlmondNrows;
  _mds.Nhalo = std::max(_mds.Nhalo, parAlmondNhalo);
}
