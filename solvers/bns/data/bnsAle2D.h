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

//mean flow
#define RBAR 1.0
#define UBAR 0.2
#define VBAR 0.0

//Heaving Airfoil
#define H 0.08
#define FREQ 0.01
#define PI 3.14159265

// Initial conditions
#define bnsInitialConditions2D(c, nu, t, x, y, r, u, v, s11, s12, s22) \
{                                         \
  *(r) = RBAR;                            \
  *(u) = UBAR;                            \
  *(v) = VBAR;                            \
  *(s11) = 0.0;                           \
  *(s12) = 0.0;                           \
  *(s22) = 0.0;                           \
}

// // Initial conditions
// #define bnsInitialConditions2D(c, nu, t, x, y, r, u, v, s11, s12, s22) \
// {                                         \
//   *(r) = 1 + exp(-3*(x*x+y*y));           \
//   *(u) = exp(-3*(x*x+y*y));               \
//   *(v) = exp(-3*(x*x+y*y));               \
//   *(s11) = 0.0;                           \
//   *(s12) = 0.0;                           \
//   *(s22) = 0.0;                           \
// }

// Body force
#define bnsBodyForce2D(c, nu, t, x, y, r, u, v, fx, fy) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define bnsBoundaryConditions2D(bc, c, nu, \
                                t, x, y, nx, ny, \
                                rM, uM, vM, s11M, s12M, s22M, \
                                rB, uB, vB, s11B, s12B, s22B) \
{                                      \
  if(bc==1){                           \
    *(rB) = rM;                        \
    *(uB) = 0.0;                       \
    *(vB) = 0.0;                       \
    *(s11B) = 0.0;                     \
    *(s12B) = 0.0;                     \
    *(s22B) = 0.0;                     \
  } else if(bc==2){                    \
    *(rB) = RBAR;                      \
    *(uB) = UBAR;                      \
    *(vB) = VBAR;                      \
    *(s11B) = 0.0;                     \
    *(s12B) = 0.0;                     \
    *(s22B) = 0.0;                     \
  } else if(bc==3){                    \
    *(rB) = RBAR;                      \
    *(uB) = uM;                        \
    *(vB) = vM;                        \
    *(s11B) = s11M;                    \
    *(s12B) = s12M;                    \
    *(s22B) = s22M;                    \
  } else if(bc==4||bc==5){             \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM)*nx;     \
    *(vB) = vM - (nx*uM+ny*vM)*ny;     \
    *(s11B) = s11M;                    \
    *(s12B) = s12M;                    \
    *(s22B) = s22M;                    \
  }                                    \
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define bnsAleBoundaryConditions2D(bc, c, nu, \
                                t, x, y, nx, ny, \
                                mVx, mVy, \
                                rM, uM, vM, s11M, s12M, s22M, \
                                rB, uB, vB, s11B, s12B, s22B) \
{                                      \
  if(bc==1){                           \
    *(rB) = rM;                        \
    *(uB) = mVx;                       \
    *(vB) = mVy;                       \
    *(s11B) = 0.0;                     \
    *(s12B) = 0.0;                     \
    *(s22B) = 0.0;                     \
  } else if(bc==2){                    \
    *(rB) = RBAR;                      \
    *(uB) = UBAR;                      \
    *(vB) = VBAR;                      \
    *(s11B) = 0.0;                     \
    *(s12B) = 0.0;                     \
    *(s22B) = 0.0;                     \
  } else if(bc==3){                    \
    *(rB) = RBAR;                      \
    *(uB) = uM;                        \
    *(vB) = vM;                        \
    *(s11B) = s11M;                    \
    *(s12B) = s12M;                    \
    *(s22B) = s22M;                    \
  } else if(bc==4||bc==5){             \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM)*nx;     \
    *(vB) = vM - (nx*uM+ny*vM)*ny;     \
    *(s11B) = s11M;                    \
    *(s12B) = s12M;                    \
    *(s22B) = s22M;                    \
  }                                    \
}

// Mesh Deformation boundary
/* moving wall 1, stationary 2 */
#define bnsMeshBoundary2D(bc, t, x, y, mx, my) \
{ \
  if(bc==1){                                 \
    *(mx)=0.0;                               \
    *(my)=H*2*PI*FREQ*cos(2*PI*FREQ*t);      \
  } else if(bc==2){                          \
    *(mx)=0.0;                               \
    *(my)=0.0;                               \
  } else if(bc==3){                          \
    *(mx)=0.0;                               \
    *(my)=0.0;                               \
  } else if(bc==6){                          \
    *(mx)=0.0;                               \
    *(my)=0.0;                               \
  }                                          \
}
