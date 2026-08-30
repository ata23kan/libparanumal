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
#include "timer.hpp"

void mds_t::Run(){

  //setup linear algebra module
  platform.linAlg().InitKernels({"set"});

  //setup linear solver
  hlong NglobalDofs;
  NglobalDofs = ogsMasked.NgatherGlobal*Nfields;


  linearSolver_t<dfloat> linearSolver;
  if (settings.compareSetting("LINEAR SOLVER","NBPCG")){
    linearSolver.Setup<LinearSolver::nbpcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","NBFPCG")){
    linearSolver.Setup<LinearSolver::nbfpcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PCG")){
    linearSolver.Setup<LinearSolver::pcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PGMRES")){
    linearSolver.Setup<LinearSolver::pgmres<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PMINRES")){
    linearSolver.Setup<LinearSolver::pminres<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  }

  properties_t kernelInfo = mesh.props; //copy base occa properties

  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DMDS "/data/mdsBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DMDS "/data/mdsBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int Nmax = std::max(mesh.Np, mesh.Nfaces*mesh.Nfp);
  kernelInfo["defines/" "p_Nmax"]= Nmax;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  
  std::string oklFilePrefix = DMDS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // fileName   = oklFilePrefix + "mdsRhs" + suffix + oklFileSuffix;
  // kernelName = "mdsRhs" + suffix;
  // kernel_t forcingKernel = platform.buildKernel(fileName, kernelName,
  //                                                   kernelInfo);

  kernel_t rhsBCKernel, addBCKernel;

  if(settings.compareSetting("DEFORMATION METHOD", "LAPLACIAN")){
    fileName   = oklFilePrefix + "mdsRhsBCLaplacian" + suffix + oklFileSuffix;
    kernelName = "mdsRhsBCLaplacian" + suffix;
  
    rhsBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

    fileName   = oklFilePrefix + "mdsAddBCLaplacian" + suffix + oklFileSuffix;
    kernelName = "mdsAddBCLaplacian" + suffix;

    addBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
  } else if(settings.compareSetting("DEFORMATION METHOD", "LINEARELASTIC")){
    fileName   = oklFilePrefix + "mdsRhsBCLinElastic" + suffix + oklFileSuffix;
    kernelName = "mdsRhsBCLinElastic" + suffix;
  
    rhsBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

    fileName   = oklFilePrefix + "mdsAddBCLinElastic" + suffix + oklFileSuffix;
    kernelName = "mdsAddBCLinElastic" + suffix;

    addBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  }


  //create occa buffers
  dlong Nall = Nfields*mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);

  memory<dfloat> rxL(Nall);
  memory<dfloat> xL(Nall);
  memory<dfloat> ryL(Nall);
  memory<dfloat> yL(Nall);
  memory<dfloat> rzL(Nall);
  memory<dfloat> zL(Nall);

  deviceMemory<dfloat> o_rxL = platform.reserve<dfloat>(Nall);
  deviceMemory<dfloat> o_xL = platform.reserve<dfloat>(Nall);
  deviceMemory<dfloat> o_ryL;
  deviceMemory<dfloat> o_yL;
  deviceMemory<dfloat> o_rzL;
  deviceMemory<dfloat> o_zL;

  deviceMemory<dfloat> o_rx, o_x, o_ry, o_y, o_rz, o_z;
  dlong Ng = ogsMasked.Ngather;
  dlong Nghalo = gHalo.Nhalo;
  dlong Ngall  = Nfields*(Ng+Nghalo);
  o_rx = platform.reserve<dfloat>(Ngall);
  o_x  = platform.reserve<dfloat>(Ngall);

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  //Set x to zero
  platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_xL);

  if(deform_laplace){
    o_ryL = platform.reserve<dfloat>(Nall);
    o_yL  = platform.reserve<dfloat>(Nall);
    platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_yL);
    o_ry = platform.reserve<dfloat>(Ngall);
    o_y  = platform.reserve<dfloat>(Ngall);
    if(mesh.dim==3){
      o_rzL = platform.reserve<dfloat>(Nall);
      o_zL  = platform.reserve<dfloat>(Nall);
      platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_zL);
      o_rz = platform.reserve<dfloat>(Ngall);
      o_z  = platform.reserve<dfloat>(Ngall);
    }
  }

  rhsBCKernel(mesh.Nelements,
              mesh.o_wJ,
              mesh.o_ggeo,
              mesh.o_sgeo,
              mesh.o_vgeo,
              mesh.o_D,
              mesh.o_S,
              mesh.o_Se,
              mesh.o_MM,
              mesh.o_vmapM,
              mesh.o_sM,
              lambda,
              mu,
              o_gamma,
              mesh.o_x,
              mesh.o_y,
              mesh.o_z,
              o_mapB,
              o_rxL,
              o_ryL,
              o_rzL);

  // gather rhs to globalDofs if c0
  ogsMasked.Gather(o_rx, o_rxL, Nfields, ogs::Add, ogs::Trans);
  ogsMasked.Gather(o_x, o_xL, Nfields, ogs::Add, ogs::NoTrans);

  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-8 : 1.0e-5;

  // timePoint_t start = GlobalPlatformTime(platform);
  timePoint_t start;
  int iter_u;
  if(deform_laplace){
    ogsMasked.Gather(o_ry, o_ryL, 1, ogs::Add, ogs::Trans);
    ogsMasked.Gather(o_y, o_yL, 1, ogs::Add, ogs::NoTrans);

    start = GlobalPlatformTime(platform);
    iter_u = Solve(linearSolver, o_x, o_rx, tol, maxIter, verbose);
    int iter_v = Solve(linearSolver, o_y, o_ry, tol, maxIter, verbose);

    ogsMasked.Scatter(o_xL, o_x, 1, ogs::NoTrans);
    ogsMasked.Scatter(o_yL, o_y, 1, ogs::NoTrans);

    if(mesh.dim==3){
      ogsMasked.Gather(o_rz, o_rzL, 1, ogs::Add, ogs::Trans);
      ogsMasked.Gather(o_z, o_zL, 1, ogs::Add, ogs::NoTrans);
      int iter_w = Solve(linearSolver, o_z, o_rz, tol, maxIter, verbose);
      ogsMasked.Scatter(o_zL, o_z, 1, ogs::NoTrans);
    }
  } else if(deform_linElastic){

    start = GlobalPlatformTime(platform);

    iter_u = Solve(linearSolver, o_x, o_rx, tol, maxIter, verbose);

    ogsMasked.Scatter(o_xL, o_x, Nfields, ogs::NoTrans);
  }

  memory<dfloat> Q(mesh.dim*mesh.Np*mesh.Nelements);
  deviceMemory<dfloat> o_Q;
  o_Q = platform.reserve<dfloat>(mesh.dim*mesh.Np*mesh.Nelements);

  //fill masked nodes with BC data
  addBCKernel(mesh.Nelements,
              mesh.o_x,
              mesh.o_y,
              mesh.o_z,
              o_mapB,
              o_Q,
              o_xL,
              o_yL,
              o_zL);

  timePoint_t end = GlobalPlatformTime(platform);
  double elapsedTime = ElapsedTime(start, end);

  if ((mesh.rank==0) && verbose){
    printf("%d, " hlongFormat ", %g, %d, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
           mesh.N,
           NglobalDofs,
           elapsedTime,
           iter_u,
           elapsedTime/(NglobalDofs),
           NglobalDofs*((dfloat)iter_u/elapsedTime),
           (char*) settings.getSetting("PRECONDITIONER").c_str());
  }

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_u_%04d.vtu", name.c_str(), mesh.rank);

    o_Q.copyTo(Q);
    PlotNewMesh2(Q, fname);
  }

  // output norm of final solution
  {
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_MxL = platform.reserve<dfloat>(Nentries);
    mesh.MassMatrixApply(o_xL, o_MxL);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_xL, o_MxL, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}
