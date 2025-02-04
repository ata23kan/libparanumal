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

mds_t mds_t::SetupNewDegree(mesh_t& meshC){

  //if asking for the same degree, return the original solver
  if (meshC.N == mesh.N) return *this;

  //shallow copy
  mds_t mds = *this;

  mds.mesh = meshC;

  /*setup trace halo exchange */
  mds.traceHalo = meshC.HaloTraceSetup(Nfields);

  //setup boundary flags and make mask and masked ogs
  mds.BoundarySetup();

  // OCCA build stuff
  properties_t kernelInfo = meshC.props; //copy base occa properties

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  std::string oklFilePrefix = DMDS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (meshC.dim==2)
    boundaryHeaderFileName = std::string(DMDS "/data/mdsBoundary2D.h");
  else if (meshC.dim==3)
    boundaryHeaderFileName = std::string(DMDS "/data/mdsBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1,blockMax/meshC.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  properties_t kernelInfoDouble = kernelInfo;
  kernelInfoDouble["defines/dfloat"]= "double";
  kernelInfoDouble["defines/dfloat4"]= "double4";

  properties_t kernelInfoFloat = kernelInfo;
  kernelInfoFloat["defines/dfloat"]= "float";
  kernelInfoFloat["defines/dfloat4"]= "float4";

  // Ax kernel
  fileName   = oklFilePrefix + "mdsAx" + suffix + oklFileSuffix;
  if(meshC.elementType==Mesh::HEXAHEDRA){
    if(mesh.settings.compareSetting("ELEMENT MAP", "TRILINEAR"))
      kernelName = "mdsPartialAxTrilinear" + suffix;
    else
      kernelName = "mdsPartialAx" + suffix;
  } else{
    kernelName = "mdsPartialAx" + suffix;
  }

  mds.partialAxKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfoDouble);

  mds.floatPartialAxKernel = platform.buildKernel(fileName, kernelName,
                                                         kernelInfoFloat);

  mds.Ndofs = mds.ogsMasked.Ngather*Nfields;
  mds.Nhalo = mds.gHalo.Nhalo*Nfields;

  mds.precon = precon_t();

  return mds;
}
