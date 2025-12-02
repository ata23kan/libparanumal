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
#include "mdsPrecon.hpp"

void mds_t::Setup(platform_t& _platform, mesh_t& _mesh,
                       settings_t& _settings, dfloat _lambda, dfloat _mu,
                       const int _NBCTypes, const memory<int> _BCType){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;
  lambda = _lambda;
  mu = _mu;

  if (settings.compareSetting("DEFORMATION METHOD", "LAPLACIAN")){
    Nfields = 1;
    deform_laplace = 1;
    deform_linElastic = 0;
  } else if (settings.compareSetting("DEFORMATION METHOD", "LINEARELASTIC")){
    Nfields = (mesh.dim==3) ? 3:2;
    deform_laplace = 0;
    deform_linElastic = 1;
  }

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);
  ogs::InitializeKernels(platform, ogs::Pfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"add", "sum", "scale",
        "axpy", "zaxpy",
        "amx", "amxpy", "zamxpy",
        "adx", "adxpy", "zadxpy",
        "innerProd", "norm2", "d2p", "p2d"});


  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  // Boundary Type translation. Just defaults.
  NBCTypes = _NBCTypes;
  BCType.malloc(NBCTypes);
  BCType.copyFrom(_BCType);

  //setup boundary flags and make mask and masked ogs
  BoundarySetup();

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DMDS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DMDS "/data/mdsBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DMDS "/data/mdsBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  properties_t kernelInfoDouble = kernelInfo;
  kernelInfoDouble["defines/dfloat"] = "double";
  kernelInfoDouble["defines/dfloat4"] = "double4";

  properties_t kernelInfoFloat = kernelInfo;
  kernelInfoFloat["defines/dfloat"] = "float";
  kernelInfoFloat["defines/dfloat4"] = "float4";

  // Ax kernel
  if (settings.compareSetting("DEFORMATION METHOD","LAPLACIAN")) {
    fileName   = oklFilePrefix + "mdsAxLaplacian" + suffix + oklFileSuffix;
    // if(mesh.elementType==Mesh::HEXAHEDRA){
    //   if(mesh.settings.compareSetting("ELEMENT MAP", "TRILINEAR"))
    //   kernelName = "mdsPartialAxTrilinear" + suffix;
    //   else
    //     kernelName = "mdsPartialAxLaplacian" + suffix;
    // } else{
      kernelName = "mdsPartialAxLaplacian" + suffix;
    // }
  } else if (settings.compareSetting("DEFORMATION METHOD", "LINEARELASTIC")){
    fileName = oklFilePrefix + "mdsAxLinElastic" + suffix + oklFileSuffix;
    kernelName = "mdsPartialAxLinElastic" + suffix;
  }

    partialAxKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfoDouble);


    floatPartialAxKernel = platform.buildKernel(fileName, kernelName,
                                                kernelInfoFloat);


  /* Preconditioner Setup */
  Ndofs = ogsMasked.Ngather*Nfields;
  Nhalo = gHalo.Nhalo*Nfields;

  if(deform_laplace){
    if (settings.compareSetting("PRECONDITIONER", "NONE"))
      precon.Setup<IdentityPrecon>(Ndofs);
    else if(settings.compareSetting("PRECONDITIONER", "PARALMOND"))
      precon.Setup<ParAlmondPrecon>(*this);
  } else if (deform_linElastic){
    precon.Setup<IdentityPrecon>(Ndofs);
  }

  // if (settings.compareSetting("PRECONDITIONER", "JACOBI"))
  //   precon.Setup<JacobiPrecon>(*this);
  // else if(settings.compareSetting("PRECONDITIONER", "MASSMATRIX"))
  //   precon.Setup<MassMatrixPrecon>(*this);
  // else if(settings.compareSetting("PRECONDITIONER", "MULTIGRID"))
  //   precon.Setup<MultiGridPrecon>(*this);
  // else if(settings.compareSetting("PRECONDITIONER", "SEMFEM"))
  //   precon.Setup<SEMFEMPrecon>(*this);
  // else if(settings.compareSetting("PRECONDITIONER", "OAS"))
  //   precon.Setup<OASPrecon>(*this);
  // else if(settings.compareSetting("PRECONDITIONER", "NONE"))
  //   precon.Setup<IdentityPrecon>(Ndofs);
}
