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

#ifndef MDS_HPP
#define MDS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "linAlg.hpp"
#include "precon.hpp"
#include "linearSolver.hpp"
#include "parAlmond.hpp"

#define DMDS LIBP_DIR"/solvers/mds/"

using namespace libp;

class mdsSettings_t: public settings_t {
public:
  mdsSettings_t() = default;
  mdsSettings_t(const comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};
void mdsAddRunSettings(settings_t& settings);
void mdsAddSettings(settings_t& settings,
                         const std::string prefix="");

class mds_t: public solver_t {
public:
  mesh_t mesh;

  dlong Ndofs, Nhalo;
  int Nfields;

  dfloat lambda;
  dfloat mu;
  dfloat tau;

  memory<dfloat> gamma;
  deviceMemory<dfloat> o_gamma;

  int deform_laplace, deform_linElastic;

  ogs::halo_t traceHalo;

  precon_t precon;

  // NOTE pfloat
  memory<pfloat> weight, weightG;
  deviceMemory<pfloat> o_weight, o_weightG;

  //C0-FEM mask data
  ogs::ogs_t ogsMasked;
  ogs::halo_t gHalo;
  memory<int> mapB;      // boundary flag of face nodes
  deviceMemory<int> o_mapB;

  dlong Nmasked;
  memory<dlong> maskIds;
  memory<hlong> maskedGlobalIds;
  memory<hlong> maskedGlobalNumbering;
  memory<dlong> GlobalToLocal;

  deviceMemory<dlong> o_maskIds;
  deviceMemory<dlong> o_GlobalToLocal;

  int NBCTypes;
  memory<int> BCType;
  memory<int> EToB;
  deviceMemory<int> o_EToB;

  int allNeumann;
  dfloat allNeumannPenalty;
  dfloat allNeumannScale;

  kernel_t maskKernel;
  kernel_t partialAxKernel;
  kernel_t partialGradientKernel;
  kernel_t partialIpdgKernel;
  // kernel_t AxKernel;

  kernel_t floatPartialAxKernel;
  kernel_t floatPartialGradientKernel;
  kernel_t floatPartialIpdgKernel;


  mds_t() = default;
  mds_t(platform_t &_platform, mesh_t &_mesh,
              settings_t& _settings, dfloat _lambda, dfloat _mu,
              const int _NBCTypes, const memory<int> _BCType) {
    Setup(_platform, _mesh, _settings, _lambda, _mu, _NBCTypes, _BCType);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             settings_t& _settings, dfloat _lambda, dfloat _mu,
             const int _NBCTypes, const memory<int> _BCType);

  void BoundarySetup();

  void SetupDiffusionCoefficient();
  void SetupDiffusionCoefficientTri2D();

  void Run();

  int Solve(linearSolver_t<dfloat>& linearSolver, deviceMemory<dfloat> &o_x, deviceMemory<dfloat> &o_r,
            const dfloat tol, const int MAXIT, const int verbose);

  void PlotFields(memory<dfloat>& Q, std::string fileName);
  void PlotNewMesh(memory<dfloat>& Qx, memory<dfloat>& Qy, std::string fileName);
  void PlotNewMesh2(memory<dfloat>& Q, std::string fileName);

  void Operator(deviceMemory<double>& o_q, deviceMemory<double>& o_Aq);
  void Operator(deviceMemory<float>& o_q, deviceMemory<float>& o_Aq);

  void BuildOperatorMatrix(parAlmond::parCOO& A);
  void BuildOperatorMatrixLaplacianTri2D(parAlmond::parCOO& A);
  void BuildOperatorMatrixLaplacianTet3D(parAlmond::parCOO& A);

  void BuildOperatorDiagonal(memory<dfloat>& diagA);
  void BuildOperatorDiagonalLaplacianTri2D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalLaplacianTet3D(memory<dfloat>& diagA);

  mds_t SetupNewDegree(mesh_t& meshF);

  mds_t SetupRingPatch(mesh_t& meshPatch);

  void ZeroMean(deviceMemory<double> &o_q);
  void ZeroMean(deviceMemory<float> &o_q);
};


#endif

