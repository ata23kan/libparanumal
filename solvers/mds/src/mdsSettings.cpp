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

//settings for mds solver
mdsSettings_t::mdsSettings_t(const comm_t& _comm):
  settings_t(_comm) {

  //common settings used when the mds solver
  // is used inside another solver
  mdsAddSettings(*this);
  parAlmond::AddSettings(*this);
}

void mdsAddRunSettings(settings_t& settings) {
  settings.newSetting("DATA FILE",
                      "data/mdsTest2D.h",
                      "Boundary and Initial conditions header");

  settings.newSetting("LAMBDA",
                      "1.0",
                      "First Lame Parameter");

  settings.newSetting("MU",
                      "0.35",
                      "Second Lame Parameter");

  settings.newSetting("OUTPUT TO FILE",
                      "FALSE",
                      "Flag for writing fields to VTU files",
                      {"TRUE", "FALSE"});

  settings.newSetting("OUTPUT FILE NAME",
                      "mds");
}

void mdsAddSettings(settings_t& settings,
                         const std::string prefix) {
  settings.newSetting(prefix+"DISCRETIZATION",
                      "CONTINUOUS",
                      "Type of Finite Element Discretization",
                      {"CONTINUOUS"});

  settings.newSetting(prefix+"DEFORMATION METHOD",
                      "LAPLACIAN",
                      "Type of Finite Element Discretization",
                      {"LAPLACIAN", "LINEARELASTIC"});

  settings.newSetting(prefix+"LINEAR SOLVER",
                      "PCG",
                      "Iterative Linear Solver to use for solve",
                      {"PCG", "FPCG", "NBPCG", "NBFPCG", "PGMRES", "PMINRES"});

  settings.newSetting(prefix+"LINEAR SOLVER STOPPING CRITERION",
                      "ABS/REL-INITRESID",
                      "Stopping criterion for the linear solver",
                      {"ABS/REL-INITRESID", "ABS/REL-RHS-2NORM"});

  settings.newSetting(prefix+"PRECONDITIONER",
                      "NONE",
                      "Preconditioning Strategy",
                      {"NONE", "JACOBI", "PARALMOND"});

  settings.newSetting(prefix+"VERBOSE",
                      "FALSE",
                      "Enable verbose output",
                      {"TRUE", "FALSE"});
}

void mdsSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "mds Settings:\n\n";
    reportSetting("DATA FILE");

    reportSetting("LAMBDA");
    reportSetting("DISCRETIZATION");
    reportSetting("LINEAR SOLVER");
    reportSetting("PRECONDITIONER");

    if (compareSetting("PRECONDITIONER","PARALMOND"))
      parAlmond::ReportSettings(*this);

    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");
  }
}

void mdsSettings_t::parseFromFile(platformSettings_t& platformSettings,
                                       meshSettings_t& meshSettings,
                                       const std::string filename) {
  //read all settings from file
  settings_t s(comm);
  s.readSettingsFromFile(filename);

  for(auto it = s.settings.begin(); it != s.settings.end(); ++it) {
    setting_t& set = it->second;
    const std::string name = set.getName();
    const std::string val = set.getVal<std::string>();
    if (platformSettings.hasSetting(name))
      platformSettings.changeSetting(name, val);
    else if (meshSettings.hasSetting(name))
      meshSettings.changeSetting(name, val);
    else if (hasSetting(name)) //self
      changeSetting(name, val);
    else  {
      LIBP_FORCE_ABORT("Unknown setting: [" << name << "] requested");
    }
  }
}
