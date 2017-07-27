  
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void meshSurfaceGeometricFactorsTetP3D(mesh3D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 5;
  mesh->sgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nsgeo*mesh->Nfaces, sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    iint id = e*mesh->Nverts;
    dfloat xe1 = mesh->EX[id+0], ye1 = mesh->EY[id+0], ze1 = mesh->EZ[id+0];
    dfloat xe2 = mesh->EX[id+1], ye2 = mesh->EY[id+1], ze2 = mesh->EZ[id+1];
    dfloat xe3 = mesh->EX[id+2], ye3 = mesh->EY[id+2], ze3 = mesh->EZ[id+2];
    dfloat xe4 = mesh->EX[id+3], ye4 = mesh->EY[id+3], ze4 = mesh->EZ[id+3];

    /* Jacobian matrix */
    dfloat xr = 0.5*(xe2-xe1), xs = 0.5*(xe3-xe1), xt = 0.5*(xe4-xe1);
    dfloat yr = 0.5*(ye2-ye1), ys = 0.5*(ye3-ye1), yt = 0.5*(ye4-ye1);
    dfloat zr = 0.5*(ze2-ze1), zs = 0.5*(ze3-ze1), zt = 0.5*(ze4-ze1);

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
    dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
    dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
    dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

    if(J<0) printf("bugger: got negative geofac\n");
    
    /* face 1 */
    iint base = mesh->Nsgeo*mesh->Nfaces*e;
    dfloat nx1 = -tx;
    dfloat ny1 = -ty;
    dfloat nz1 = -tz;
    dfloat sJ1 = norm(nx1,ny1,nz1);

    mesh->sgeo[base+NXID] = nx1/sJ1;
    mesh->sgeo[base+NYID] = ny1/sJ1;
    mesh->sgeo[base+NZID] = nz1/sJ1;
    mesh->sgeo[base+SJID] = sJ1*J;
    mesh->sgeo[base+IJID] = 1./J;

    /* face 2 */
    base += mesh->Nsgeo;
    dfloat nx2 = -sx;
    dfloat ny2 = -sy;
    dfloat nz2 = -sz;
    dfloat sJ2 = norm(nx2,ny2,nz2);

    mesh->sgeo[base+NXID] = nx2/sJ2;
    mesh->sgeo[base+NYID] = ny2/sJ2;
    mesh->sgeo[base+NZID] = nz2/sJ2;
    mesh->sgeo[base+SJID] = sJ2*J;
    mesh->sgeo[base+IJID] = 1./J;

    /* face 3 */
    base += mesh->Nsgeo;
    dfloat nx3 = rx+sx+tx;
    dfloat ny3 = ry+sy+ty;
    dfloat nz3 = rz+sz+tz;
    dfloat sJ3 = norm(nx3,ny3,nz3);

    mesh->sgeo[base+NXID] = nx3/sJ3;
    mesh->sgeo[base+NYID] = ny3/sJ3;
    mesh->sgeo[base+NZID] = nz3/sJ3;
    mesh->sgeo[base+SJID] = sJ3*J;
    mesh->sgeo[base+IJID] = 1./J;

    /* face 4 */
    base += mesh->Nsgeo;
    dfloat nx4 = -rx;
    dfloat ny4 = -ry;
    dfloat nz4 = -rz;
    dfloat sJ4 = norm(nx4,ny4,nz4);

    mesh->sgeo[base+NXID] = nx4/sJ4;
    mesh->sgeo[base+NYID] = ny4/sJ4;
    mesh->sgeo[base+NZID] = nz4/sJ4;
    mesh->sgeo[base+SJID] = sJ4*J;
    mesh->sgeo[base+IJID] = 1./J;
  }
}
