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

// Overlapping additive Schwarz with patch problems consisting of the
//  entire local mesh + 1 ring overlap, solved with a local multigrid
//  precon and coarse problem consisting of the global degree 1
//  problem, solved with parAlmond
void OASPrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  pfloat one = 1., zero = 0;
  
  if (mesh.N>1) {
    deviceMemory<pfloat> o_rPatch = mds.platform.reserve<pfloat>(mdsPatch.Ndofs);

    linAlg_t& linAlg = mds.platform.linAlg();

    if (mds.disc_c0) {
      dlong Ntotal = mesh.Np*(mesh.Nelements+mesh.totalRingElements);
      deviceMemory<pfloat> o_rPatchL = mds.platform.reserve<pfloat>(Ntotal);
      //Scatter to localDof ordering, exchange ring halo,
      // then compress the ring mesh to globalDofs order.
      // TODO: Technically, these steps couple be fused
      // to a single operation, but currently theres no easy way
      // as the ordering of globalDofs between the original mesh
      // partition and the ring mesh could be different
      mds.ogsMasked.Scatter(o_rPatchL, o_r, 1, ogs::NoTrans);
      mesh.ringHalo.Exchange(o_rPatchL, mesh.Np);
      mdsPatch.ogsMasked.Gather(o_rPatch, o_rPatchL, 1, ogs::Add, ogs::NoTrans);
    } else {
      o_rPatch.copyFrom(o_r, mds.Ndofs, 0, properties_t("async", true));
      mesh.ringHalo.Exchange(o_rPatch, mesh.Np);
    }

    //Apply local patch precon
    // TODO: This is blocking due to H<->D transfers.
    //       Should modify precons so size=1 is non-blocking
    deviceMemory<pfloat> o_zPatch = mds.platform.reserve<pfloat>(mdsPatch.Ndofs);
    preconPatch.Operator(o_rPatch, o_zPatch);

    dlong NcolsC = parAlmond.getNumCols(0);
    deviceMemory<pfloat> o_rC = mds.platform.reserve<pfloat>(NcolsC);
    deviceMemory<pfloat> o_zC = mds.platform.reserve<pfloat>(NcolsC);

    //Coarsen problem to N=1 and pass to parAlmond
    level.coarsen(o_r, o_rC);

    parAlmond.Operator(o_rC, o_zC);

    //Add contributions from all patches together
    if (mds.disc_c0) {
      dlong Ntotal = mesh.Np*(mesh.Nelements+mesh.totalRingElements);
      deviceMemory<pfloat> o_zPatchL = mds.platform.reserve<pfloat>(Ntotal);

      mdsPatch.ogsMasked.Scatter(o_zPatchL, o_zPatch, 1, ogs::NoTrans);
      ogsMaskedRing.GatherScatter(o_zPatchL, 1, ogs::Add, ogs::Sym);

      // Weight by overlap degree, zPatch = patchWeight*zPatch
      Ntotal=mesh.Nelements*mesh.Np;
      linAlg.amx(Ntotal, one, o_patchWeight, o_zPatchL);

      mds.ogsMasked.Gather(o_Mr, o_zPatchL, 1, ogs::Add, ogs::NoTrans);

    } else {
      mesh.ringHalo.Combine(o_zPatch, mesh.Np);

      // Weight by overlap degree, Mr = patchWeight*zPatch
      linAlg.amxpy(mds.Ndofs, one, o_patchWeight, o_zPatch, zero, o_Mr);
    }

    // Add prologatated coarse solution
    level.prolongate(o_zC, o_Mr);
  } else {
    //if N=1 just call the coarse solver
    parAlmond.Operator(o_r, o_Mr);
  }

  // zero mean of RHS
  if(mds.allNeumann) mds.ZeroMean(o_Mr);

}

OASPrecon::OASPrecon(mds_t& _mds):
  mds(_mds), mesh(_mds.mesh), settings(_mds.settings),
  parAlmond(mds.platform, settings, mesh.comm) {

  //build the one ring mesh
  if (mesh.N>1) {
    if (Comm::World().rank()==0){
      printf("-----------------------------Multigrid Degree %2d Patch--------------------------------------\n", mesh.N);
    }
    meshPatch = mesh.SetupRingPatch();
    mdsPatch = mds.SetupRingPatch(meshPatch);
    preconPatch.Setup<MultiGridPrecon>(mdsPatch);

    //compute patch overlap weighting
    patchWeight.malloc(meshPatch.Nelements*meshPatch.Np);
    for (int i=0;i<meshPatch.Nelements*meshPatch.Np;i++)
      patchWeight[i] = 1.0;

    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      //share the masked version of the global id numbering
      memory<hlong> maskedRingGlobalIds(meshPatch.Nelements*meshPatch.Np);
      maskedRingGlobalIds.copyFrom(mds.maskedGlobalIds, mesh.Nelements*mesh.Np);
      mesh.ringHalo.Exchange(maskedRingGlobalIds, mesh.Np);

      //mask ring
      for (dlong n=0;n<mdsPatch.Nmasked;n++)
        maskedRingGlobalIds[mdsPatch.maskIds[n]] = 0;

      //use the masked ids to make another gs handle
      int verbose = 0;
      bool unique = true; //flag a unique node in every gather node
      ogsMaskedRing.Setup(meshPatch.Nelements*meshPatch.Np,
                          maskedRingGlobalIds, mesh.comm,
                          ogs::Signed, ogs::Auto,
                          unique, verbose, mds.platform);

      //determine overlap of each node with masked ogs
      ogsMaskedRing.GatherScatter(patchWeight, 1, ogs::Add, ogs::Sym);

    } else {
      //determine overlap by combining halos
      mesh.ringHalo.Combine(patchWeight, mesh.Np);
    }

    //invert
    for (int i=0;i<meshPatch.Nelements*meshPatch.Np;i++)
      patchWeight[i] = (patchWeight[i] > 0.0) ? 1.0/patchWeight[i] : 0.0;

    o_patchWeight = mds.platform.malloc<pfloat>(patchWeight);
  }

  //build the coarse precon
  int Nc = 1;  //hard code
  int NpCoarse = mesh.Np;
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      NpCoarse = ((Nc+1)*(Nc+2))/2; break;
    case Mesh::QUADRILATERALS:
      NpCoarse = (Nc+1)*(Nc+1); break;
    case Mesh::TETRAHEDRA:
      NpCoarse = ((Nc+1)*(Nc+2)*(Nc+3))/6; break;
    case Mesh::HEXAHEDRA:
      NpCoarse = (Nc+1)*(Nc+1)*(Nc+1); break;
  }

  //build mesh and mds objects for this degree
  mesh_t meshC = mesh.SetupNewDegree(Nc);
  mds_t mdsC = mds.SetupNewDegree(meshC);

  //build full A matrix and pass to parAlmond
  if (Comm::World().rank()==0){
    printf("-----------------------------Multigrid AMG Setup--------------------------------------------\n");
  }
  parAlmond::parCOO A(mds.platform, meshC.comm);
  if (settings.compareSetting("DISCRETIZATION", "IPDG"))
    mdsC.BuildOperatorMatrixIpdg(A);
  else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS"))
    mdsC.BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = meshC.rank;
  int size = meshC.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = (dlong) (A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);
  memory<pfloat> null(numLocalRows);
  for (dlong i=0;i<numLocalRows;i++) {
    null[i] = 1.0/sqrt(TotalRows);
  }

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond.AMGSetup(A, mdsC.allNeumann, null,mdsC.allNeumannPenalty);

  if (mesh.N>1) {
    //make an MG level to get prologation and coarsener
    dlong Nrows, Ncols;
    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      Nrows = mds.ogsMasked.Ngather;
      Ncols = Nrows + mds.gHalo.Nhalo;
    } else {
      Nrows = mesh.Nelements*mesh.Np;
      Ncols = Nrows + mesh.totalHaloPairs*mesh.Np;
    }

    level = MGLevel(mds, Nrows, Ncols, Nc, NpCoarse);
    level.meshC = meshC;
    level.mdsC = mdsC;
  }

  //report
  parAlmond.Report();
}
