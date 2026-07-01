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

void bns_t::Report(dfloat time, int tstep){

  static int frame=0;
  static int error_frame = 0;
  static int force_frame = 0;

  //compute q.M*q
  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
  deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
  mesh.MassMatrixApply(o_q, o_Mq);

  dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));
  o_Mq.free();

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.4f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    //compute vorticity
    deviceMemory<dfloat> o_Vort = platform.reserve<dfloat>(mesh.dim*mesh.Nelements*mesh.Np);
    vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_Vort);

    memory<dfloat> Vort(mesh.dim*mesh.Nelements*mesh.Np);

    //compute Q-criterion
    deviceMemory<dfloat> o_QCrit = platform.reserve<dfloat>(mesh.Nelements*mesh.Np);
    qcriterionKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_QCrit);
    memory<dfloat> QCrit(mesh.Nelements*mesh.Np);

    // copy data back to host
    o_q.copyTo(q);
    o_Vort.copyTo(Vort);
    o_QCrit.copyTo(QCrit);

    // Copy the new mesh positions from ALE
    mesh.o_x.copyTo(mesh.x);
    mesh.o_y.copyTo(mesh.y);
    if(mesh.dim==3)
      mesh.o_z.copyTo(mesh.z);

    if(testCase==4){
      char fnameTgv[BUFSIZ];
      char nameTgv[] = "KE";
      sprintf(fnameTgv, "%s_%04d_%04d.txt", nameTgv, mesh.rank, mesh.Np);
      EnergyTGV(q, std::string(fnameTgv), time);
    }

    std::string name;
    // output field files
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];

    // //Re36000
    // // Output interval is 0.071429. Tolerance is 0.07142857 / 2.
    // const dfloat output_tolerance = 0.035714;

    // // Snapshot 1: Start of 5th cycle (t/T = 0)
    // if(time < (28.571429 + output_tolerance) && time > (28.571429 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    // // Snapshot 2: Quarter of 5th cycle (t/T = 1/4)
    // if(time < (30.357143 + output_tolerance) && time > (30.357143 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    // // Snapshot 3: Half of 5th cycle (t/T = 2/4)
    // if(time < (32.142857 + output_tolerance) && time > (32.142857 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    // // Snapshot 4: Three-quarters of 5th cycle (t/T = 3/4)
    // if(time < (33.928571 + output_tolerance) && time > (33.928571 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    //Re17000
    //// Output interval is 0.068966. Tolerance is 0.068966 / 2.
    //const dfloat output_tolerance = 0.034483;
    //
    //// Snapshot 1: Start of 5th cycle (t/T = 0)
    //if(time < (27.586207 + output_tolerance) && time > (27.586207 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    //// Snapshot 2: Quarter of 5th cycle (t/T = 1/4)
    //if(time < (29.310345 + output_tolerance) && time > (29.310345 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    //// Snapshot 3: Half of 5th cycle (t/T = 2/4)
    //if(time < (31.034483 + output_tolerance) && time > (31.034483 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    //// Snapshot 4: Three-quarters of 5th cycle (t/T = 3/4)
    //if(time < (32.758621 + output_tolerance) && time > (32.758621 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    // Re5800
    // Snapshot 1: Start of 5th cycle (t/T = 0)
    //if(time < (23.529412 + 1e-05) && time > (23.529412 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}
    //
    //// Snapshot 2: Quarter of 5th cycle (t/T = 1/4)
    //if(time < (25.000000 + 1e-05) && time > (25.000000 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}
    //
    //// Snapshot 3: Half of 5th cycle (t/T = 2/4)
    //if(time < (26.470588 + 1e-05) && time > (26.470588 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}
    //
    //// Snapshot 4: Three-quarters of 5th cycle (t/T = 3/4)
    //if(time < (27.941176 + 1e-05) && time > (27.941176 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    // Re5800 (u_sim = 0.1, Coarse Output: 20 points/cycle)
    // Target: 9th Period Snapshots
    const dfloat output_tolerance = 0.147059; 
    
    // Snapshot 1: Start of 9th cycle (t = 8.00 * T)
    if(time < (47.058824 + output_tolerance) && time > (47.058824 - output_tolerance)){
      sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
      PlotFields(q, QCrit, std::string(fname));
    }

    // Snapshot 2: Quarter of 9th cycle (t = 8.25 * T)
    if(time < (48.529412 + output_tolerance) && time > (48.529412 - output_tolerance)){
      sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
      PlotFields(q, QCrit, std::string(fname));
    }

    // Snapshot 3: Half of 9th cycle (t = 8.50 * T)
    if(time < (50.000000 + output_tolerance) && time > (50.000000 - output_tolerance)){
      sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
      PlotFields(q, QCrit, std::string(fname));
    }

    // Snapshot 4: Three-quarters of 9th cycle (t = 8.75 * T)
    if(time < (51.470588 + output_tolerance) && time > (51.470588 - output_tolerance)){
      sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
      PlotFields(q, QCrit, std::string(fname));
    }

    // Re 2400
    // Snapshot 1: Start of 5th cycle (t/T = 0)
    //if(time < (19.047619 + 1e-05) && time > (19.047619 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    //// Snapshot 2: Quarter of 5th cycle (t/T = 1/4)
    //if(time < (20.238095 + 1e-05) && time > (20.238095 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    //// Snapshot 3: Half of 5th cycle (t/T = 2/4)
    //if(time < (21.428571 + 1e-05) && time > (21.428571 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    //// Snapshot 4: Three-quarters of 5th cycle (t/T = 3/4)
    //if(time < (22.619048 + 1e-05) && time > (22.619048 - 1e-05)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}

    // // Re2400 (u_sim = 0.1, Coarse Output: 20 points/cycle)
    // // Target: 9th Period Snapshots
    // const dfloat output_tolerance = 0.119048; 
    
    // // Snapshot 1: Start of 9th cycle (t/T = 0)
    // if(time < (38.095238 + output_tolerance) && time > (38.095238 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    // // Snapshot 2: Quarter of 9th cycle (t/T = 1/4)
    // if(time < (39.285714 + output_tolerance) && time > (39.285714 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    // // Snapshot 3: Half of 9th cycle (t/T = 2/4)
    // if(time < (40.476190 + output_tolerance) && time > (40.476190 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    // // Snapshot 4: Three-quarters of 9th cycle (t/T = 3/4)
    // if(time < (41.666667 + output_tolerance) && time > (41.666667 - output_tolerance)){
    //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //   PlotFields(q, QCrit, std::string(fname));
    // }

    //Re720
    //const dfloat output_tolerance = 0.017857;
    //// Snapshot 1: Start of 5th cycle (t/T = 0)
    //if(time < (14.285714 + output_tolerance) && time > (14.285714 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}
    //
    //// Snapshot 2: Quarter of 5th cycle (t/T = 1/4)
    //if(time < (15.178571 + output_tolerance) && time > (15.178571 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}
    //
    //// Snapshot 3: Half of 5th cycle (t/T = 2/4)
    //if(time < (16.071429 + output_tolerance) && time > (16.071429 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}
    //
    //// Snapshot 4: Three-quarters of 5th cycle (t/T = 3/4)
    //if(time < (16.964286 + output_tolerance) && time > (16.964286 - output_tolerance)){
    //  sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
    //  PlotFields(q, QCrit, std::string(fname));
    //}
      
    if(testCase==2){
      ComputeForces(time);
    }
    if(testCase==5){
      ComputeForces(time);
    }
    
    if(testCase==1){
      char fname2[BUFSIZ];
      sprintf(fname2, "%s_error_%04d_%04d.vtu", name.c_str(), mesh.rank, error_frame++);

      PlotConstantError(q, 1.f, 1.f, 1.f, std::string(fname2));
      constantErrorNorm(q, 1.f, 1.f, 1.f);
    }

  }

  /*
  if(bns->dim==3){
    if(options.compareArgs("OUTPUT FILE FORMAT","ISO")){

      for (int gr=0; gr<bns->isoGNgroups; gr++){

        bns->isoNtris[0] = 0;
        bns->o_isoNtris.copyFrom(bns->isoNtris);
        if(mesh->nonPmlNelements){
        bns->isoSurfaceKernel(mesh->nonPmlNelements,    // Numner of elements
                              mesh->o_nonPmlElementIds,    // Element Ids
                              bns->isoField,               // which field to use for isosurfacing
                              bns->isoColorField,          // which field to use for isosurfacing
                              bns->isoGNlevels[gr],        // number of isosurface levels
                              bns->o_isoGLvalues[gr],      // array of isosurface levels
                              bns->isoMaxNtris,            // maximum number of generated triangles
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              bns->o_q,
                              bns->o_Vort,
                              bns->o_VortMag,
                              bns->o_plotInterp,
                              bns->o_plotEToV,
                              bns->o_isoNtris,             // output: number of generated triangles
                              bns->o_isoq);                // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)

      }
        // find number of generated triangles
        bns->o_isoNtris.copyTo(bns->isoNtris);
        bns->isoNtris[0] = mymin(bns->isoNtris[0], bns->isoMaxNtris);

        //
        printf("Rank:%2d Group:%2d Triangles:%8d\n", mesh->rank, bns->isoNtris[0], gr);
        //
        int offset = 0;
        bns->o_isoq.copyTo(bns->isoq, bns->isoNtris[0]*(mesh->dim+bns->isoNfields)*3*sizeof(dfloat), offset);

        char fname[BUFSIZ];
        string outName;
        options.getArgs("OUTPUT FILE NAME", outName);


        if(options.compareArgs("OUTPUT FILE FORMAT", "WELD"))
        {
          int Ntris1 = bns->isoNtris[0];
          int Ntris2 = bnsWeldTriVerts(bns, Ntris1, bns->isoq);

          printf("Welding triangles:%8d to:%8d\n", Ntris1, Ntris2);
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, gr, mesh->rank, bns->frame);
          bnsIsoWeldPlotVTU(bns,  fname);
        }
        else
        {
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, gr, mesh->rank, bns->frame);
          bnsIsoPlotVTU(bns, bns->isoNtris[0], bns->isoq, fname);
        }
      }
      bns->frame++;
    }
  }
  */
}


void bns_t::ComputeForces(const dfloat T){

  // Hard coded weights i.e. sum(MM1D,1), MM1D = inv(V1)'*V1
  
  memory<dfloat> W; W.malloc(mesh.Nfp, 0);

  memory<dfloat> Mf; Mf.malloc(mesh.Nfp*mesh.Nfp,0); 

  if(mesh.dim==2){
    const int Np = (mesh.N+1); 
    memory<dfloat> r, V;
    mesh.Nodes1D(mesh.N, r); //Gauss-Legendre-Lobatto nodes
    mesh.Vandermonde1D(mesh.N, r, V); 
    mesh.MassMatrix1D(Np, V, Mf);   
  }else if(mesh.dim==3){
    const int Np = (mesh.N+1)*(mesh.N+2)/2; 
    memory<dfloat> r, s, V;
    mesh.NodesTri2D(mesh.N, r, s); //Gauss-Legendre-Lobatto nodes
    mesh.VandermondeTri2D(mesh.N, r, s, V); 
    mesh.MassMatrixTri2D(Np, V, Mf);   
  }

  for(int i=0; i< mesh.Nfp; i++){
    dfloat sum = 0; 
      for(int j=0; j< mesh.Nfp; j++){
        sum += Mf[i + j*mesh.Nfp]; 
      }
      W[i] = sum; 
      // printf(" %.8f \n", W[i]);
  }
  //

  // ALE updates geometric factors on device, so refresh host-side copies before writing.
  mesh.o_sgeo.copyTo(mesh.sgeo);

  const dfloat rhoRef = 1.0;
  const dfloat pRef = RT*rhoRef;

  dfloat FVx = 0.0, FVy = 0.0, FVz = 0.0; 
  dfloat FPx = 0.0, FPy = 0.0, FPz = 0.0;
  dfloat FVx_fin = 0.0, FPx_fin = 0.0, FVx_body = 0.0, FPx_body = 0.0;
  const dfloat x_cutoff = 0.85; // Hard coded cutoff to separate the main body and fin surfaces
  for(int e=0;e<mesh.Nelements;++e){
  
    for(int f=0;f<mesh.Nfaces;++f){
  
      int bc = mesh.EToB[e*mesh.Nfaces+f];
  
      if(bc==1){ // Hard coded here i.e. wall
  
        for(int n=0;n<mesh.Nfp; n++){
          dfloat nx=0.0, ny=0.0, nz=0.0, sJ = 0.0; 
          const dlong sid = e*mesh.Nfaces*mesh.Nsgeo + f*mesh.Nsgeo;
          nx = mesh.sgeo[sid+mesh.NXID];
          ny = mesh.sgeo[sid+mesh.NYID];
          nz = mesh.dim==3 ? mesh.sgeo[sid+mesh.NZID]: 0.0;

          sJ = mesh.sgeo[sid+mesh.SJID];

          const dlong vid  = e*mesh.Nfp*mesh.Nfaces + f*mesh.Nfp + n;
          const dlong idM  = mesh.vmapM[vid];
          // load traces
          const int vidM = idM%mesh.Np;

          // Use mapped boundary node, not face-local index, for body/tail split.
          const dfloat x_node = mesh.x[idM];

          const dlong qidM = e*mesh.Np*Nfields + vidM;

          // Read trace values
          dfloat q1  = 0.0, q2  = 0.0, q3  = 0.0, q4  = 0.0, q5  = 0.0; 
          dfloat q6  = 0.0, q7  = 0.0, q8  = 0.0, q9  = 0.0, q10 = 0.0; 

          // Read trace values
          q1  = q[qidM + 0*mesh.Np];
          q2  = q[qidM + 1*mesh.Np];
          q3  = q[qidM + 2*mesh.Np];
          q4  = q[qidM + 3*mesh.Np];
          q5  = q[qidM + 4*mesh.Np];
          q6  = q[qidM + 5*mesh.Np];

          if(mesh.dim==3){
            q7  = q[qidM + 6*mesh.Np];
            q8  = q[qidM + 7*mesh.Np];
            q9  = q[qidM + 8*mesh.Np];
            q10 = q[qidM + 9*mesh.Np];
          }


          dfloat P, s11, s12, s13, s22, s23, s33; 

          if(mesh.dim==2){
            P  =   RT*q1; 
            s11 = -RT*(sqrt(2) * q5 - q2*q2/ q1);
            s22 = -RT*(sqrt(2) * q6 - q3*q3/ q1);
            s12 = -RT*(          q4 - q2*q3/ q1);
            
            FVx += -W[n]*sJ*(s11*nx + s12*ny);   
            FVy += -W[n]*sJ*(s12*nx + s22*ny);   
            FPx +=  W[n]*sJ*(P*nx);   
            FPy +=  W[n]*sJ*(P*ny);  

          }else if(mesh.dim==3){

            P  =   RT*(q1-rhoRef); 

            s11 = -RT*(sqrt(2) * q8 - q2*q2/ q1);
            s22 = -RT*(sqrt(2) * q9 - q3*q3/ q1);
            s33 = -RT*(sqrt(2) * q10 - q4*q4/ q1);
            
            s12 = -RT*(          q5 - q2*q3/ q1);
            s13 = -RT*(          q6 - q2*q4/ q1);
            s23 = -RT*(          q7 - q3*q4/ q1);

            const dfloat dFVx = -W[n]*sJ*(s11*nx + s12*ny + s13*nz);   
            const dfloat dFPx =  W[n]*sJ*(P*nx);

            if (x_node > x_cutoff) {
              // This integration point sits on the fin
              FVx_fin += dFVx;
              FPx_fin += dFPx;
            } else {
              // This integration point sits on the main body
              FVx_body += dFVx;
              FPx_body += dFPx;
            }
            // FVx += -W[n]*sJ*(s11*nx + s12*ny + s13*nz);   
            // FPx +=  W[n]*sJ*(P*nx);
            // FVy += -W[n]*sJ*(s12*nx + s22*ny + s23*nz);   
            // FVz += -W[n]*sJ*(s13*nx + s23*ny + s33*nz);   
            // FPy +=  W[n]*sJ*(P*ny);  
            // FPz +=  W[n]*sJ*(P*nz);  
          }           
        }
      }
    }
  }

  dfloat gFVx_fin = FVx_fin, gFPx_fin = FPx_fin, gFVx_body = FVx_body, gFPx_body = FPx_body;
  // dfloat gFVx = FVx; 
  // dfloat gFVy = FVy; 
  // dfloat gFPx = FPx; 
  // dfloat gFPy = FPy; 

  comm.Allreduce(gFVx_fin, Comm::Sum);
  comm.Allreduce(gFPx_fin, Comm::Sum);
  comm.Allreduce(gFVx_body, Comm::Sum);
  comm.Allreduce(gFPx_body, Comm::Sum);

  // comm.Allreduce(gFVx, Comm::Sum);
  // comm.Allreduce(gFVy, Comm::Sum);
  // comm.Allreduce(gFPx, Comm::Sum);
  // comm.Allreduce(gFPy, Comm::Sum);

  // dfloat gFVz = 0.0, gFPz = 0.0;  
  // if(mesh.dim==3){
  //     gFVz = FVz;  
  //     gFPz = FPz; 
  //     comm.Allreduce(gFVz, Comm::Sum);
  //     comm.Allreduce(gFPz, Comm::Sum);
  // }

  if(mesh.rank==0){
    std::string force_name;
    settings.getSetting("OUTPUT FILE NAME", force_name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_BNSForceData_N%d.dat", force_name.c_str(), mesh.N);
    FILE *fp; fp = fopen(fname, "a");  
    // if(mesh.dim==2)
    //   fprintf(fp, "%.4e %.8e %.8e %.8e %.8e\n", T, gFVx,  gFVy,  gFPx,  gFPy); 
    
    if(mesh.dim==3)
      fprintf(fp, "%.4e %.8e %.8e %.8e %.8e\n", T, gFVx_fin, gFPx_fin, gFVx_body, gFPx_body); 
      // fprintf(fp, "%.4e %.8e %.8e %.8e %.8e %.8e %.8e\n", T, gFVx,  gFVy, gFVz, gFPx, gFPy, gFPz); 

    fclose(fp);
  }
}