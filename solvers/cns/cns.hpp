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

#ifndef CNS_HPP
#define CNS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "stab.hpp"

#define DCNS LIBP_DIR"/solvers/cns/"

using namespace libp;

class cnsSettings_t: public settings_t {
public:
  cnsSettings_t(comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);

  stabSettings_t extractStabSettings();
};



class cns_t: public solver_t {
public:
  mesh_t mesh;
  stab_t stab {};
  stabSettings_t stabSettings; 

  int Nfields;
  int Ngrads;

  dfloat Re, Ma, Pr;
  dfloat mu, cp, cv, R;
  dfloat gamma, igamma, gammaM1, igammaM1, gammaP1, igammaP1;

  int viscType;
  int Nph; // number of physical parameters 
  int MUID, GMID, RRID, PRID, CPID, CVID; 
  int EXID, TRID, TSID, CSID; 

  // Physical coefficients
  memory<dfloat> refState; 

  // Physical coefficients
  memory<dfloat> pCoeff; 
  deviceMemory<dfloat> o_pCoeff; 

  int cubature;
  int isothermal;

  timeStepper_t timeStepper;

  ogs::halo_t fieldTraceHalo;
  ogs::halo_t gradTraceHalo;



  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;
  kernel_t cubatureVolumeKernel;
  kernel_t cubatureSurfaceKernel;

  kernel_t gradVolumeKernel;
  kernel_t gradSurfaceKernel;

  kernel_t vorticityKernel;

  kernel_t constrainKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  cns_t() = default;
  cns_t(platform_t &_platform, mesh_t &_mesh,
              cnsSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             cnsSettings_t& _settings);

  void Run();

  // Set reference values, thermodynamics, nondimensional values
  void setPhysics(properties_t & kernelInfo);
  void setArtificialDiffusion(properties_t & kernelInfo);

  // Set detectors
  void SetDetector(properties_t & kernelInfo);


  void Report(dfloat time, int tstep) override;

  void PlotFields(memory<dfloat> Q, memory<dfloat> V, std::string fileName);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void rhsNoStab(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void rhsArtDiff(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
};

#endif
