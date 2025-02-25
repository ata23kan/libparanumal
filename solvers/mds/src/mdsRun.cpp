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
  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    NglobalDofs = ogsMasked.NgatherGlobal*Nfields;
  } else {
    NglobalDofs = mesh.NelementsGlobal*mesh.Np*Nfields;
  }

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

  fileName   = oklFilePrefix + "mdsRhs" + suffix + oklFileSuffix;
  kernelName = "mdsRhs" + suffix;
  kernel_t forcingKernel = platform.buildKernel(fileName, kernelName,
                                                    kernelInfo);

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

  memory<dfloat> ruL(Nall);
  memory<dfloat> rvL(Nall);
  memory<dfloat> xuL(Nall);
  memory<dfloat> xvL(Nall);
  deviceMemory<dfloat> o_ruL = platform.reserve<dfloat>(Nall);
  deviceMemory<dfloat> o_xuL = platform.reserve<dfloat>(Nall);
  
  // deviceMemory<dfloat> o_rvL;
  deviceMemory<dfloat> o_xvL;

  // o_rvL = platform.reserve<dfloat>(Nall);
  // o_xvL = platform.reserve<dfloat>(Nall);

  deviceMemory<dfloat> o_ru, o_rv, o_xu, o_xv;
  dlong Ng = ogsMasked.Ngather;
  dlong Nghalo = gHalo.Nhalo;
  dlong Ngall  = Nfields*(Ng+Nghalo);
  o_ru = platform.reserve<dfloat>(Ngall);
  o_xu = platform.reserve<dfloat>(Ngall);

  // if(settings.compareSetting("DEFORMATION METHOD", "LAPLACIAN")){
  //   o_rv = platform.reserve<dfloat>(Ngall);
  //   o_xv = platform.reserve<dfloat>(Ngall);
  // }

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  //Set x to zero
  platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_xuL);

  // if(settings.compareSetting("DEFORMATION METHOD", "LAPLACIAN")){
  //   platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_xvL);
  // }

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
              mesh.o_x,
              mesh.o_y,
              mesh.o_z,
              o_mapB,
              o_ruL);
              // o_rvL);

  // gather rhs to globalDofs if c0
  ogsMasked.Gather(o_ru, o_ruL, Nfields, ogs::Add, ogs::Trans);
  ogsMasked.Gather(o_xu, o_xuL, Nfields, ogs::Add, ogs::NoTrans);

  // if(settings.compareSetting("DEFORMATION METHOD","LAPLACIAN")){
  //   ogsMasked.Gather(o_rv, o_rvL, Nfields, ogs::Add, ogs::Trans);
  //   ogsMasked.Gather(o_xv, o_xvL, Nfields, ogs::Add, ogs::NoTrans);
  // }

  int maxIter = 50;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;

  timePoint_t start = GlobalPlatformTime(platform);

  //call the solver
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-3 : 1.0e-5;
  int iter_u = Solve(linearSolver, o_xu, o_ru, tol, maxIter, verbose);

  // if(settings.compareSetting("DEFORMATION METHOD","LAPLACIAN")){
  //   int iter_v = Solve(linearSolver, o_xv, o_rv, tol, maxIter, verbose);
  // }

  //add the boundary data to the masked nodes
  // scatter x to LocalDofs if c0
  ogsMasked.Scatter(o_xuL, o_xu, Nfields, ogs::NoTrans);

  // if(settings.compareSetting("DEFORMATION METHOD","LAPLACIAN")){
  //   ogsMasked.Scatter(o_xvL, o_xv, Nfields, ogs::NoTrans);
  // }

  deviceMemory<dfloat> o_Q;
  o_Q = platform.reserve<dfloat>(Nfields*mesh.Np*mesh.Nelements);
  // if(settings.compareSetting("DEFORMATION METHOD", "LINEARELASTIC")){
  // }

  //fill masked nodes with BC data
  addBCKernel(mesh.Nelements,
              mesh.o_x,
              mesh.o_y,
              mesh.o_z,
              o_mapB,
              o_Q,
              o_xuL,
              o_xvL);

  
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

    o_Q.copyTo(xuL);

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_u_%04d.vtu", name.c_str(), mesh.rank);
    // PlotNewMesh(xuL, xvL, fname);
    PlotNewMesh2(xuL, fname);

    // sprintf(fname, "%s_v_%04d.vtu", name.c_str(), mesh.rank);
    // PlotFields(xvL, fname);
  }

  // output norm of final solution
  {
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_MxL = platform.reserve<dfloat>(Nentries);
    mesh.MassMatrixApply(o_xuL, o_MxL);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_xuL, o_MxL, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}
