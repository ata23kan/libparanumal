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

// Inverse Mass Matrix preconditioner
MassMatrixPrecon::MassMatrixPrecon(mds_t& _mds):
  mds(_mds), mesh(_mds.mesh), settings(_mds.settings) {

  //sanity checking
  LIBP_ABORT("MASSMATRIX preconditioner is only available for triangle and tetrhedra elements. Use JACOBI instead.",
             mesh.elementType!=Mesh::TRIANGLES && mesh.elementType!=Mesh::TETRAHEDRA);

  LIBP_ABORT("MASSMATRIX preconditioner is unavailble when lambda=0.",
             mds.lambda==0);

  o_pfloat_invMM = mds.platform.malloc<pfloat>(mesh.pfloat_invMM);

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  int blockMax = 256;
  if (mds.platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;
  kernelInfo["defines/" "dfloat"]= pfloatString;

  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {
    blockJacobiKernel = mds.platform.buildKernel(DMDS "/okl/mdsPreconBlockJacobi.okl",
                                     "blockJacobi", kernelInfo);
  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    partialBlockJacobiKernel = mds.platform.buildKernel(DMDS "/okl/mdsPreconBlockJacobi.okl",
                                     "partialBlockJacobi", kernelInfo);
  }
}

void MassMatrixPrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  pfloat one = 1.0, zero = 0.0;

  pfloat invLambda = 1./mds.lambda;
  
  linAlg_t& linAlg = mds.platform.linAlg();

  if (mds.disc_c0) {//C0
    dlong Ntotal = mds.ogsMasked.Ngather + mds.gHalo.Nhalo;
    deviceMemory<pfloat> o_rtmp = mds.platform.reserve<pfloat>(Ntotal);
    deviceMemory<pfloat> o_MrL  = mds.platform.reserve<pfloat>(mesh.Np*mesh.Nelements);

    // rtmp = invDegree.*r
    linAlg.amxpy(mds.Ndofs, one, mds.o_weightG, o_r, zero, o_rtmp);

    mds.gHalo.ExchangeStart(o_rtmp, 1);

    if(mesh.NlocalGatherElements/2)
      partialBlockJacobiKernel(mesh.NlocalGatherElements/2,
                               mesh.o_localGatherElementList,
                               mds.o_GlobalToLocal,
                               invLambda, mesh.o_pfloat_vgeo, o_pfloat_invMM,
                               o_rtmp, o_MrL);

    // finalize halo exchange
    mds.gHalo.ExchangeFinish(o_rtmp, 1);

    if(mesh.NglobalGatherElements)
      partialBlockJacobiKernel(mesh.NglobalGatherElements,
                               mesh.o_globalGatherElementList,
                               mds.o_GlobalToLocal,
                               invLambda, mesh.o_pfloat_vgeo, o_pfloat_invMM,
                               o_rtmp, o_MrL);

    //gather result to Aq
    mds.ogsMasked.GatherStart(o_Mr, o_MrL, 1, ogs::Add, ogs::Trans);

    if((mesh.NlocalGatherElements+1)/2){
      partialBlockJacobiKernel((mesh.NlocalGatherElements+1)/2,
                               mesh.o_localGatherElementList+mesh.NlocalGatherElements/2,
                               mds.o_GlobalToLocal,
                               invLambda, mesh.o_pfloat_vgeo, o_pfloat_invMM,
                               o_rtmp, o_MrL);
    }

    mds.ogsMasked.GatherFinish(o_Mr, o_MrL, 1, ogs::Add, ogs::Trans);

    // Mr = invDegree.*Mr
    linAlg.amx(mds.Ndofs, one, mds.o_weightG, o_Mr);

  } else {
    //IPDG
    blockJacobiKernel(mesh.Nelements, invLambda, mesh.o_pfloat_vgeo, o_pfloat_invMM, o_r, o_Mr);
  }

  // zero mean of RHS
  if(mds.allNeumann) mds.ZeroMean(o_Mr);
}
