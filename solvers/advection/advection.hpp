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

#ifndef ADVECTION_HPP
#define ADVECTION_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "mds.hpp"

#define DADVECTION LIBP_DIR"/solvers/advection/"

using namespace libp;

class advectionSettings_t: public settings_t {
public:
  advectionSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);

  mdsSettings_t extractMdsSettings();
};

class advection_t: public solver_t {
public:
  mesh_t mesh;
  timeStepper_t timeStepper;

  int plotcount = 0;

  ogs::halo_t traceHalo;

  // Mesh deformation for ALE
  mdsSettings_t mdsSettings;
  mds_t mdsSolver;
  linearSolver_t<dfloat> mdsLinearSolver;
  dfloat mdsTOL, lambda, mu;
  int mdsNfields;
  memory<dfloat> mapB;  // boundary flag of face nodes
  deviceMemory<dfloat> o_mapB;
  int Niter;

  mesh_t meshN1;
  memory<dfloat> IM;          // interpolation to higher order
  deviceMemory<dfloat> o_IM;

  memory<dfloat> q, meshVelx, meshVely;
  deviceMemory<dfloat> o_q, o_meshVelx, o_meshVely;
  deviceMemory<dfloat> o_VX, o_VX0;  // Vertex positions to be updated in ALE
  deviceMemory<dfloat> o_uxV, o_uyV, o_uzV; // Mesh velocity at the verteces obtained by elliptic solve

  kernel_t volumeKernel;
  kernel_t surfaceKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  kernel_t aleVolumeKernel, aleSurfaceKernel;
  kernel_t aleRhsKernel, aleBCKernel;
  kernel_t updateVgeoKernel, updateSgeoKernel;
  kernel_t velInterpolationKernel, posInterpolationKernel;
  kernel_t initialPositionKernel;
  kernel_t explicitDeformationKernel;

  kernel_t convertPrimitiveKernel, convertConservativeKernel;

  advection_t() = default;
  advection_t(platform_t &_platform, mesh_t &_mesh,
              advectionSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& platform, mesh_t& mesh,
             advectionSettings_t& settings);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Q, const std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  void MeshSolve(deviceMemory<dfloat>& o_VX, deviceMemory<dfloat>& o_rhsX, const dfloat T);

  void UpdateGeo(deviceMemory<dfloat>& o_VX);

  void UpdateX(deviceMemory<dfloat>& o_VX);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);
};

#endif
