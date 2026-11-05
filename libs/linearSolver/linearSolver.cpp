/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Anthony Austin

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

#include "linearSolver.hpp"

namespace libp {

template <typename T>
int linearSolver_t<T>::Solve(operator_t& linearOperator,
			     operator_t& precon,
			     deviceMemory<T>& o_x,
			     deviceMemory<T>& o_rhs,
			     const T tol,
			     const int MAXIT,
			     const int verbose) {
  assertInitialized();

  const int timing = ls->timing;
  timePoint_t lap;
  if (timing) lap = PlatformTime(ls->platform);

  ig->FormInitialGuess(o_x, o_rhs);
  LinearSolver::TimeLap(ls->platform, timing, LinearSolver::TIME_IG_FORM, lap);

  int iters = ls->Solve(linearOperator, precon, o_x, o_rhs, tol, MAXIT, verbose);

  // the solver has charged its own phases; restart the lap after it
  if (timing) lap = PlatformTime(ls->platform);

  ig->Update(linearOperator, o_x, o_rhs);
  LinearSolver::TimeLap(ls->platform, timing, LinearSolver::TIME_IG_UPDATE, lap);

  if (timing) {
    LinearSolver::solveNsolves++;
    LinearSolver::solveNiterations += iters;
  }

  return iters;
}

template <typename T>
void linearSolver_t<T>::MakeDefaultInitialGuessStrategy() {
  ig = std::make_shared<InitialGuess::Zero<T> >(ls->N, ls->platform, ls->settings, ls->comm);
}

template <typename T>
bool linearSolver_t<T>::isInitialized() {
  return (ls!=nullptr && ig!=nullptr);
}

template <typename T>
void linearSolver_t<T>::assertInitialized() {
  LIBP_ABORT("LinearSolver not initialized",
             ls==nullptr);
  LIBP_ABORT("InitialGuess not initialized",
             ig==nullptr);
}

template class linearSolver_t<double>;
template class linearSolver_t<float>;

namespace LinearSolver {

double    solveTime[NTIMES]  = {0.0};
long long solveCalls[NTIMES] = {0};
long long solveNsolves = 0, solveNiterations = 0;

void TimeLap(platform_t& platform, const int timing, const int id, timePoint_t& lap) {
  if (!timing) return;
  timePoint_t now = PlatformTime(platform);
  solveTime[id] += ElapsedTime(lap, now);
  solveCalls[id]++;
  lap = now;
}

void PrintSolveTimes(comm_t comm) {

  static const char* names[NTIMES] = {
    "initial guess (form)",
    "pool reserve",
    "A*x (initial residual)",
    "axpy (r=b-Ax, p=z+beta p)",
    "preconditioner",
    "norms / inner products",
    "A*p",
    "UpdatePCG (x, r, r.r)",
    "initial guess (update)"};

  // per-rank totals differ by load imbalance; report the slowest rank
  memory<double> times(NTIMES);
  for (int n=0;n<NTIMES;++n) times[n] = solveTime[n];
  comm.Allreduce(times, Comm::Max);

  double total = 0.0;
  for (int n=0;n<NTIMES;++n) total += times[n];

  if (comm.rank()==0) {
    printf("Linear solve phase times (max over ranks): %lld solves, %lld iterations (%.2f per solve)\n",
           solveNsolves, solveNiterations,
           (solveNsolves>0) ? static_cast<double>(solveNiterations)/solveNsolves : 0.0);
    printf("  %-28s %12s %7s %12s %10s\n", "phase", "time [s]", "share", "calls", "us/call");
    for (int n=0;n<NTIMES;++n) {
      printf("  %-28s %12.4e %6.1f%% %12lld %10.2f\n", names[n], times[n],
             (total>0.0) ? 100.0*times[n]/total : 0.0, solveCalls[n],
             (solveCalls[n]>0) ? 1.0e6*times[n]/solveCalls[n] : 0.0);
    }
    printf("  %-28s %12.4e\n", "sum of phases", total);
  }
}

} //namespace LinearSolver

} //namespace libp
