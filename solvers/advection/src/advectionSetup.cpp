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

#include "advection.hpp"

void advection_t::Setup(platform_t& _platform, mesh_t& _mesh,
                         advectionSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = mesh.comm;
  settings = _settings;

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup linear algebra module
  // platform.linAlg().InitKernels({"innerProd", "max"});
  platform.linAlg().InitKernels({"innerProd", "axpy", "max", "set"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(1); //one field

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, 1, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, 1, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, 1, platform, comm);
  }

  // Setup mesh deformation solver
  // bc = 1 -> walls
  // bc = 2 -> outflow
  // bc = 3 -> wallm
  int NBCTypes = 4;
  memory<int> mdsBCType(NBCTypes);
  mdsBCType[0] = 0;
  mdsBCType[1] = 1;
  mdsBCType[2] = 1;
  mdsBCType[3] = 2;

  // Build low order mesh for deformation
  meshN1 = mesh.SetupNewDegree(1);

  // Build interpolation matrix to high order mesh
  mesh.DegreeRaiseMatrixTri2D(meshN1.N, mesh.N, IM);
  // memory<dfloat> IMT(meshN1.Np*mesh.Np);
  // linAlg_t::matrixTranspose(mesh.Np, meshN1.Np, IM, meshN1.Np, IMT, mesh.Np);
  o_IM = platform.malloc<dfloat>(IM);

  // for(int m=0; m<mesh.Np;m++){
  //   for(int n=0;n<meshN1.Np;n++){
  //     int id = m*meshN1.Np + n;
  //     printf("%f ", IM[id]);
  //   }
  //   printf("\n");
  // }

  // printf("\n");
  // for(int m=0; m<meshN1.Np;m++){
  //   for(int n=0;n<mesh.Np;n++){
  //     int id = m*mesh.Np + n;
  //     printf("%f ", IMT[id]);
  //   }
  //   printf("\n");
  // }
  // std::exit(EXIT_SUCCESS);

  mdsSettings = _settings.extractMdsSettings();

  lambda = 1.0;
  mu = 0.35;

  mdsSolver.Setup(platform, meshN1, mdsSettings,
                  lambda, mu, NBCTypes, mdsBCType);

  mdsNfields = mdsSolver.Nfields;

  dlong mdsNLocal = mdsSolver.Ndofs;
  dlong mdsNhalo = mdsSolver.Nhalo;

  if (mdsSettings.compareSetting("LINEAR SOLVER","NBPCG")){
    mdsLinearSolver.Setup<LinearSolver::nbpcg<dfloat> >(mdsNLocal, mdsNhalo, platform, mdsSettings, comm);
  } else if (mdsSettings.compareSetting("LINEAR SOLVER","NBFPCG")){
    mdsLinearSolver.Setup<LinearSolver::nbfpcg<dfloat> >(mdsNLocal, mdsNhalo, platform, mdsSettings, comm);
  } else if (mdsSettings.compareSetting("LINEAR SOLVER","PCG")){
    mdsLinearSolver.Setup<LinearSolver::pcg<dfloat> >(mdsNLocal, mdsNhalo, platform, mdsSettings, comm);
  } else if (mdsSettings.compareSetting("LINEAR SOLVER","PGMRES")){
    mdsLinearSolver.Setup<LinearSolver::pgmres<dfloat> >(mdsNLocal, mdsNhalo, platform, mdsSettings, comm);
  } else if (mdsSettings.compareSetting("LINEAR SOLVER","PMINRES")){
    mdsLinearSolver.Setup<LinearSolver::pminres<dfloat> >(mdsNLocal, mdsNhalo, platform, mdsSettings, comm);
  }

  // solver tolerances

  //Solver tolerances
  if (sizeof(dfloat)==sizeof(double)) {
    mdsTOL = 1.0E-3;
  } else {
    mdsTOL = 1.0E-3;
  }

  // mesh velocity at the interpolation nodes
  meshVelx.malloc(Nlocal+Nhalo);
  meshVely.malloc(Nlocal+Nhalo);
  o_meshVelx = platform.malloc<dfloat>(Nlocal+Nhalo);
  o_meshVely = platform.malloc<dfloat>(Nlocal+Nhalo);

  o_dx = platform.reserve<dfloat>(meshN1.Np*meshN1.Nelements);
  o_dy = platform.reserve<dfloat>(meshN1.Np*meshN1.Nelements);
  if (mesh.dim==3){
    o_dz = platform.reserve<dfloat>(meshN1.Np*meshN1.Nelements);
  }

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo);
  o_q = platform.malloc<dfloat>(Nlocal+Nhalo);

  mesh.MassMatrixKernelSetup(1); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties
  properties_t kernelInfoN1 = meshN1.props;

  kernelInfo["defines/" "p_NpN1"] = meshN1.Np;

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;
  kernelInfoN1["includes"] += dataFileName;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int maxNodesN1 = std::max(meshN1.Np, (meshN1.Nfp*meshN1.Nfaces));
  kernelInfoN1["defines/" "p_maxNodes"]= maxNodesN1;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockVN1 = std::max(1, blockMax/meshN1.Np);
  kernelInfoN1["defines/" "p_NblockV"]= NblockV;
  kernelInfoN1["defines/" "p_NblockVN1"]= NblockVN1;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int NblockSN1 = std::max(1, blockMax/maxNodesN1);
  kernelInfoN1["defines/" "p_NblockS"]= NblockSN1;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  std::string oklFilePrefix = DADVECTION "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  kernelInfo["defines/ p_Nfields"] = mdsNfields;
  kernelInfoN1["defines/ p_Nfields"] = mdsNfields;

  int Nmax = std::max(meshN1.Np, meshN1.Nfaces*meshN1.Nfp);
  kernelInfoN1["defines/" "p_Nmax"]= Nmax;

  // Mesh Deformation kernels
  fileName   = oklFilePrefix + "advectionAleRhs" + suffix + oklFileSuffix;
  if (mdsSettings.compareSetting("DEFORMATION METHOD", "LINEARELASTIC")){
    kernelName = "aleRhsLinElastic" + suffix;
    aleRhsKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
  }else if(mdsSettings.compareSetting("DEFORMATION METHOD", "LAPLACIAN")){

    kernelName = "aleRhsLaplace" + suffix;
    aleRhsKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
  }


  kernelName  = "aleBC" + suffix;
  aleBCKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);

  fileName  = oklFilePrefix + "advectionUpdateGeometricFactors" + suffix + oklFileSuffix;
  kernelName = "updateVgeo" + suffix;
  updateVgeoKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "updateSgeo" + suffix;
  updateSgeoKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName = oklFilePrefix + "advectionInterpolateDeformation" + suffix + oklFileSuffix;
  kernelName = "interpolateVelocity" + suffix;
  velInterpolationKernel = platform.buildKernel(fileName, kernelName, kernelInfo); // kernelInfo of high order

  kernelName = "interpolatePosition" + suffix;
  posInterpolationKernel = platform.buildKernel(fileName, kernelName, kernelInfo); // kernelInfo of high order

  // kernels from volume file
  fileName   = oklFilePrefix + "advectionVolume" + suffix + oklFileSuffix;
  kernelName = "advectionVolume" + suffix;

  volumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "advectionAleVolume" + suffix;
  aleVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  // kernels from surface file
  fileName   = oklFilePrefix + "advectionSurface" + suffix + oklFileSuffix;
  kernelName = "advectionSurface" + suffix;

  surfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "advectionAleSurface" + suffix;
  aleSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);


  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "advectionInitialCondition2D" + oklFileSuffix;
    kernelName = "advectionInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "advectionInitialCondition3D" + oklFileSuffix;
    kernelName = "advectionInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "advectionMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "advectionMaxWaveSpeed" + suffix;

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
}
