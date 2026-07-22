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

void bns_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  bnsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //get physical paramters
  settings.getSetting("SPEED OF SOUND", c);
  settings.getSetting("VISCOSITY", nu);
  RT     = c*c;
  tauInv = RT/nu;

  Nfields    = (mesh.dim==3) ? 10:6;
  Npmlfields = mesh.dim*Nfields;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup cubature
  mesh.CubatureSetup();

  //Setup PML
  PmlSetup();

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;

  semiAnalytic = 0;
  if (settings.compareSetting("TIME INTEGRATOR","SARK4")
    ||settings.compareSetting("TIME INTEGRATOR","SARK5")
    ||settings.compareSetting("TIME INTEGRATOR","SAAB3")
    ||settings.compareSetting("TIME INTEGRATOR","MRSAAB3"))
    semiAnalytic = 1;

  //semi-analytic exponential coefficients
  memory<dfloat> lambda(Nfields);
  for (int i=0;i<mesh.dim+1;i++) lambda[i] = 0.0;
  for (int i=mesh.dim+1;i<Nfields;i++) lambda[i] = -tauInv;

  //make array of time step estimates for each element
  memory<dfloat> EtoDT(mesh.Nelements);
  dfloat vmax = MaxWaveSpeed();
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat h = mesh.ElementCharacteristicLength(e);
    dfloat dtAdv  = h/(vmax*(mesh.N+1.)*(mesh.N+1.));
    dfloat dtVisc = 1.0/tauInv;

    if (semiAnalytic)
      EtoDT[e] = dtAdv;
    else
      EtoDT[e] = std::min(dtAdv, dtVisc);

    /*
    Artificial warping of time step size for multirate testing
    */
#if 0
    dfloat x=0., y=0, z=0;
    for (int v=0;v<mesh.Nverts;v++){
      x += mesh.EX[e*mesh.Nverts+v];
      y += mesh.EY[e*mesh.Nverts+v];
      if (mesh.dim==3)
        z += mesh.EZ[e*mesh.Nverts+v];
    }
    x /=mesh.Nverts;
    y /=mesh.Nverts;
    z /=mesh.Nverts;

    dfloat f = std::min(fabs(x),fabs(y));
    if (mesh.dim==3)
      f = std::min(f,fabs(z));

    f = std::max(0.5, f);
    EtoDT[e] *= f;
#endif
  }

  mesh.mrNlevels=0;
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3")) {
    mesh.MultiRateSetup(EtoDT);
    mesh.MultiRatePmlSetup();
    multirateTraceHalo = mesh.MultiRateHaloTraceSetup(Nfields);
  }

  if (settings.compareSetting("TIME INTEGRATOR","MRAB3")){
    timeStepper.Setup<TimeStepper::mrab3>(mesh.Nelements, mesh.NpmlElements,
                                          mesh.totalHaloPairs,
                                          mesh.Np, Nfields, Npmlfields,
                                          platform, mesh);
  } else if (settings.compareSetting("TIME INTEGRATOR","MRSAAB3")){
    timeStepper.Setup<TimeStepper::mrsaab3>(mesh.Nelements, mesh.NpmlElements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, Nfields, Npmlfields,
                                            lambda, platform, mesh);
  } else if (settings.compareSetting("TIME INTEGRATOR","SAAB3")) {
    timeStepper.Setup<TimeStepper::saab3>(mesh.Nelements, mesh.NpmlElements,
                                          mesh.totalHaloPairs,
                                          mesh.Np, Nfields, Npmlfields,
                                          lambda, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements, mesh.NpmlElements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, Nfields, Npmlfields,
                                        platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements, mesh.NpmlElements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, Npmlfields,
                                           platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements, mesh.NpmlElements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, Npmlfields,
                                           platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK4")) {
    timeStepper.Setup<TimeStepper::sark4>(mesh.Nelements, mesh.NpmlElements,
                                          mesh.totalHaloPairs,
                                          mesh.Np, Nfields, Npmlfields,
                                          lambda, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK5")) {
    timeStepper.Setup<TimeStepper::sark5>(mesh.Nelements, mesh.NpmlElements,
                                          mesh.totalHaloPairs,
                                          mesh.Np, Nfields, Npmlfields,
                                          lambda, platform, comm);
  } else {
    LIBP_FORCE_ABORT("Requested TIME INTEGRATOR not found.");
  }

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  // Setup mesh deformation solver
  // bc = 1 -> moving boundary, bc = 2 -> stationary boundary 
  // wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, physical bounding box 6
  int NBCTypes = 11;
  memory<int> mdsBCType(NBCTypes);
  mdsBCType[0] = 0;
  mdsBCType[1] = 1;
  mdsBCType[2] = 2;
  mdsBCType[3] = 3;
  mdsBCType[4] = 2;
  mdsBCType[5] = 2;
  mdsBCType[6] = 2;
  mdsBCType[7] = 2;
  mdsBCType[8] = 2;
  mdsBCType[9] = 2;
  mdsBCType[10] = 2; // Stationary wall boundary for flexible rod

  // Build low order mesh for deformation
  meshN1 = mesh.SetupNewDegree(1);
  properties_t kernelInfoN1 = meshN1.props;

  // Build interpolation matrix to high order mesh
  if(mesh.dim==2){
    mesh.DegreeRaiseMatrixTri2D(meshN1.N, mesh.N, IM);
  } else if(mesh.dim==3){
    mesh.DegreeRaiseMatrixTet3D(meshN1.N, mesh.N, IM);
  }
  o_IM = platform.malloc<dfloat>(IM);

  mdsSettings = _settings.extractMdsSettings();

  mdsLambda = 1.0; // TODO: Why are these not coming from the settings -AA
  mdsMu = 0.35;    // TODO: Why are these not coming from the settings -AA

  mdsSolver.Setup(platform, meshN1, mdsSettings,
                  mdsLambda, mdsMu, NBCTypes, mdsBCType);

  mdsNfields = mdsSolver.Nfields;
  kernelInfoN1["defines/" "p_Nfields"] = mdsNfields;

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

  //Solver tolerances
  if (sizeof(dfloat)==sizeof(double)) {
    mdsTOL = 1.0E-3;
  } else {
    mdsTOL = 1.0E-3;
  }

  // Setup ALE velocity
  dlong NlocalAle = mesh.Nelements*mesh.Np;
  dlong NhaloAle = mesh.totalHaloPairs*mesh.Np;

  // printf("Nelements: %d\n", mesh.Nelements);
  // printf("Nelements N1: %d\n", meshN1.Nelements);
  // printf("nnonpmlNelements: %d\n", mesh.NnonPmlElements);
  // printf("nnonpmlNelements N1: %d\n", meshN1.NnonPmlElements);
  // printf("npmlNelements: %d\n", mesh.NpmlElements);
  // std::exit(EXIT_SUCCESS);

  // mesh velocity at the interpolation nodes
  meshVelx.calloc(NlocalAle+NhaloAle);
  meshVely.calloc(NlocalAle+NhaloAle);
  o_meshVelx = platform.malloc<dfloat>(meshVelx);
  o_meshVely = platform.malloc<dfloat>(meshVely);

  if(mesh.dim==3){
    meshVelz.calloc(NlocalAle+NhaloAle);
    o_meshVelz = platform.malloc<dfloat>(meshVelz);
  }  

  o_VX  = platform.reserve<dfloat>(meshN1.Np*meshN1.Nelements*mesh.dim);
  o_VX0 = platform.reserve<dfloat>(meshN1.Np*meshN1.Nelements*mesh.dim);

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo);
  o_q = platform.malloc<dfloat>(Nlocal+Nhalo);

  pmlq.malloc(mesh.NpmlElements*mesh.Np*Npmlfields);
  o_pmlq = platform.malloc<dfloat>(mesh.NpmlElements*mesh.Np*Npmlfields);

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;
  kernelInfoN1["includes"] += dataFileName;

  // ALE First order mesh properties
  kernelInfo["defines/" "p_NpN1"] = meshN1.Np;
  kernelInfo["defines/ p_NfieldsN1"] = mdsNfields;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_Npmlfields"]= Npmlfields;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int Nmax = std::max(meshN1.Np, meshN1.Nfaces*meshN1.Nfp);
  kernelInfoN1["defines/" "p_Nmax"]= Nmax;

  int blockMax = 256;
  if (platform.device.mode()=="CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockVN1 = std::max(1, blockMax/meshN1.Np);
  kernelInfoN1["defines/" "p_NblockV"]= NblockVN1;
  kernelInfo["defines/" "p_NblockVN1"]= NblockVN1;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int NblockCub = std::max(1, blockMax/mesh.cubNp);
  kernelInfo["defines/" "p_NblockCub"]= NblockCub;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DBNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "bnsVolume" + suffix + oklFileSuffix;
  kernelName = "bnsVolume" + suffix;
  volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  if (pmlcubature) {
    kernelName = "bnsPmlVolumeCub" + suffix;
    pmlVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  } else {
    kernelName = "bnsPmlVolume" + suffix;
    pmlVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  }

  // kernels from relaxation file
  fileName   = oklFilePrefix + "bnsRelaxation" + suffix + oklFileSuffix;
  kernelName = "bnsRelaxation" + suffix;
  relaxationKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  if (pmlcubature) {
    kernelName = "bnsPmlRelaxationCub" + suffix;
    pmlRelaxationKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    pmlRelaxationKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }


  // kernels from surface file
  fileName   = oklFilePrefix + "bnsSurface" + suffix + oklFileSuffix;
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3")) {
    kernelName = "bnsMRSurface" + suffix;
    surfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    kernelName = "bnsMRPmlSurface" + suffix;
    pmlSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    kernelName = "bnsSurface" + suffix;
    surfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    kernelName = "bnsPmlSurface" + suffix;
    pmlSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  // vorticity calculation
  fileName   = oklFilePrefix + "bnsVorticity" + suffix + oklFileSuffix;
  kernelName = "bnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);

  // // Q-Criterion calculation
  // fileName   = oklFilePrefix + "bnsQCriterion" + suffix + oklFileSuffix;
  // kernelName = "bnsQCriterion" + suffix;

  // qcriterionKernel = platform.buildKernel(fileName, kernelName,
  //                                    kernelInfo);                                     

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "bnsInitialCondition2D" + oklFileSuffix;
    initialConditionKernel = platform.buildKernel(fileName,
                                                  "bnsInitialCondition2D",
                                                  kernelInfo);
    pmlInitialConditionKernel = platform.buildKernel(fileName,
                                                  "bnsPmlInitialCondition2D",
                                                  kernelInfo);

    // ALE Initial Position
    kernelName = "bnsInitialPosition2D";
    initialPositionKernel  = platform.buildKernel(fileName, kernelName, kernelInfo);

  } else {
    fileName   = oklFilePrefix + "bnsInitialCondition3D" + oklFileSuffix;
    initialConditionKernel = platform.buildKernel(fileName,
                                                  "bnsInitialCondition3D",
                                                  kernelInfo);
    pmlInitialConditionKernel = platform.buildKernel(fileName,
                                                  "bnsPmlInitialCondition3D",
                                                  kernelInfo);

    // ALE Initial Position
    kernelName = "bnsInitialPosition3D";
    initialPositionKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
  }

  // ALE Kernels
  fileName  = oklFilePrefix + "bnsExplicitDeformation" + suffix + oklFileSuffix;
  if(settings.compareSetting("ALE TEST", "BOX")){
    testCase = 1;
    kernelName = "explicitDeformationBox" + suffix;
    explicitDeformationKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
  } else if(settings.compareSetting("ALE TEST", "PLUNGINGAIRFOIL")){
    testCase = 2;
    kernelName = "explicitDeformationAirfoil" + suffix;
    explicitDeformationKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);    
  } else if(settings.compareSetting("ALE TEST", "SOLVEMESH")){
    testCase = 3;
  } else if(settings.compareSetting("ALE TEST", "TGV")){
    testCase = 4;
    kernelName = "explicitDeformationTGV" + suffix;
    explicitDeformationKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
  } else if(settings.compareSetting("ALE TEST", "CARANGIFORMFISH")){
    testCase = 5;
    kernelName = "explicitDeformationFish" + suffix;
    explicitDeformationKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
  } else {
    LIBP_FORCE_ABORT("Requested ALE TEST not found.");
  }


  fileName = oklFilePrefix + "bnsInterpolateDeformation" + suffix + oklFileSuffix;
  kernelName = "interpolateVelocity" + suffix;
  velInterpolationKernel = platform.buildKernel(fileName, kernelName, kernelInfo); // kernelInfo of high order

  kernelName = "interpolatePosition" + suffix;
  posInterpolationKernel = platform.buildKernel(fileName, kernelName, kernelInfo); // kernelInfo of high order

  fileName  = oklFilePrefix + "bnsUpdateGeometricFactors" + suffix + oklFileSuffix;
  kernelName = "updateVgeo" + suffix;
  updateVgeoKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "updateSgeo" + suffix;
  updateSgeoKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "bnsAleSurface" + suffix + oklFileSuffix;
  kernelName = "bnsAleSurface" + suffix;
  aleSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "bnsAleVolume" + suffix + oklFileSuffix;
  kernelName = "bnsAleVolume" + suffix;
  aleVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "bnsAleRhs" + suffix + oklFileSuffix;
  if (mdsSettings.compareSetting("DEFORMATION METHOD", "LINEARELASTIC")){
    kernelName = "aleRhsLinElastic" + suffix;
    aleRhsKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
    kernelName = "aleBCLinElastic" + suffix;
    aleBCKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
  }else if(mdsSettings.compareSetting("DEFORMATION METHOD", "LAPLACIAN")){
    kernelName = "aleRhsLaplace" + suffix;
    aleRhsKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
    kernelName = "aleBCLaplace" + suffix;
    aleBCKernel = platform.buildKernel(fileName, kernelName, kernelInfoN1);
  }
}
