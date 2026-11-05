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
#include <limits>
#include <algorithm>

// a face keyed by its (sorted) global vertices, used to match the boundary faces
// listed in the mesh file against the faces of the local elements
typedef struct {

  dlong element; // -1 for a face read from the mesh file
  int face;      // -1 for a face read from the mesh file

  hlong v[4];    // max number of face vertices

  int bctype;    // -1 for an element face

}mdsFace_t;

void mds_t::BoundarySetup(){

  /* mesh_t::ConnectBoundary only assigns a boundary flag to element faces which
     have no neighbouring element (EToE==-1). Any face which is tagged in the
     mesh file but sits *inside* the mesh therefore loses its tag: it never finds
     a partner in the boundaryFaces list and mesh.EToB stays -1 for it. The
     physical/PML interface is exactly such a face - it is tagged in gmsh (e.g.
     "BoundingBox"), but it connects a physical element to a PML element.
     Consequently mesh.mapB carries no flag there either, which is why the
     interface was still free to move.

     Recover those tags here by matching the tagged faces of the mesh file
     against the (global) vertices of every local element face, then rebuild the
     node-wise flags from scratch, the same way mesh_t::ConnectNodes does. Note
     that mesh_t's memory<> handles are reference counted, so we work on private
     copies: writing mesh.EToB/mesh.mapB would silently change the host arrays of
     the solver which owns this mesh (bns shares them with its N=1 copy). */
  const int NfaceVertices = mesh.NfaceVertices;
  const bool verboseSetup = settings.compareSetting("VERBOSE", "TRUE");

  memory<int> gmshEToB(mesh.Nelements*mesh.Nfaces);
  gmshEToB.copyFrom(mesh.EToB, mesh.Nelements*mesh.Nfaces);

  if (mesh.NboundaryFaces>0 && mesh.boundaryInfo.length()>0) {

    /* gather the tagged faces of the mesh file and all local element faces in one
       list, keyed by their sorted global vertices, and sort it. Faces which share
       a key then sit next to each other, so a tag can be handed to both elements
       which meet at a tagged interface face. */
    const dlong Nfaces = mesh.Nelements*mesh.Nfaces;
    const dlong Ntagged = static_cast<dlong>(mesh.NboundaryFaces);

    memory<mdsFace_t> faces(Nfaces+Ntagged);

    dlong cnt=0;
    for (dlong e=0;e<mesh.Nelements;++e) {
      for (int f=0;f<mesh.Nfaces;++f) {
        for (int n=0;n<NfaceVertices;++n)
          faces[cnt].v[n] = mesh.EToV[e*mesh.Nverts + mesh.faceVertices[f*NfaceVertices+n]];
        for (int n=NfaceVertices;n<4;++n) faces[cnt].v[n] = -1;

        std::sort(faces[cnt].v, faces[cnt].v+NfaceVertices, std::less<hlong>());

        faces[cnt].element = e;
        faces[cnt].face    = f;
        faces[cnt].bctype  = -1;
        ++cnt;
      }
    }

    for (dlong b=0;b<Ntagged;++b) {
      for (int n=0;n<NfaceVertices;++n)
        faces[cnt].v[n] = mesh.boundaryInfo[b*(NfaceVertices+1)+n+1];
      for (int n=NfaceVertices;n<4;++n) faces[cnt].v[n] = -1;

      std::sort(faces[cnt].v, faces[cnt].v+NfaceVertices, std::less<hlong>());

      faces[cnt].element = -1; //flags an entry read from the mesh file
      faces[cnt].face    = -1;
      faces[cnt].bctype  = static_cast<int>(mesh.boundaryInfo[b*(NfaceVertices+1)]);
      ++cnt;
    }

    std::sort(faces.ptr(), faces.ptr()+cnt,
              [&](const mdsFace_t& a, const mdsFace_t& b) {
                return std::lexicographical_compare(a.v, a.v+NfaceVertices,
                                                    b.v, b.v+NfaceVertices);
              });

    /* walk the runs of faces sharing a key: if one of them carries a tag, hand it
       to every element face in the run */
    dlong Nrecovered = 0;
    dlong start = 0;
    while (start<cnt) {
      dlong end = start+1;
      while (end<cnt && std::equal(faces[start].v, faces[start].v+NfaceVertices,
                                   faces[end].v)) ++end;

      int bc = -1;
      for (dlong m=start;m<end;++m)
        if (faces[m].element==-1) bc = std::max(bc, faces[m].bctype);

      if (bc>0) {
        for (dlong m=start;m<end;++m) {
          const dlong e = faces[m].element;
          if (e<0) continue;
          const int f = faces[m].face;
          if (gmshEToB[e*mesh.Nfaces+f]>0) continue; //already flagged
          gmshEToB[e*mesh.Nfaces+f] = bc;
          ++Nrecovered;
        }
      }

      start = end;
    }

    if (verboseSetup) {
      dlong Ntot = Nrecovered;
      mesh.comm.Allreduce(Ntot);
      if (mesh.rank==0)
        printf("mds: recovered mesh-file boundary flags on " dlongFormat
               " interior faces (e.g. PML interface)\n", Ntot);
    }
  }

  /* node-wise boundary flags of the patched EToB (mesh_t::ConnectNodes) */
  memory<int> gmshMapB((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np, -1);

  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int f=0;f<mesh.Nfaces;f++) {
      const int bc = gmshEToB[f+e*mesh.Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh.Nfp;n++) {
          const int fid = mesh.faceNodes[n+f*mesh.Nfp];
          const int bcn = gmshMapB[fid+e*mesh.Np];
          if (bcn == -1) gmshMapB[fid+e*mesh.Np] = bc;         //no flag here yet
          else           gmshMapB[fid+e*mesh.Np] = std::min(bc,bcn);
        }
      }
    }
  }

  /* propagate the flags to every copy of a flagged node (across elements and
     ranks) until nothing changes. Without this the interface nodes are masked in
     one element and free in its neighbour, so the interface still moves. */
  hlong gatherChange = 1;
  while (gatherChange>0) {
    gatherChange = 0;

    mesh.halo.Exchange(gmshMapB, mesh.Np);

    for (dlong e=0;e<mesh.Nelements;++e) {
      for (int n=0;n<mesh.Nfp*mesh.Nfaces;++n) {
        const dlong id  = e*mesh.Nfp*mesh.Nfaces + n;
        const dlong idM = mesh.vmapM[id];
        const dlong idP = mesh.vmapP[id];
        const int bcM = gmshMapB[idM];
        const int bcP = gmshMapB[idP];

        if (bcP > 0) {
          if (bcM == -1 || bcP < bcM) {
            gmshMapB[idM] = bcP;
            ++gatherChange;
          }
        }
      }
    }

    mesh.comm.Allreduce(gatherChange);
  }

  //translate the mesh's element-to-boundaryflag mapping
  EToB.malloc(mesh.Nelements*mesh.Nfaces, 0);
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int f=0;f<mesh.Nfaces;f++) {
      int bc = gmshEToB[e*mesh.Nfaces+f];
      if (bc>0) {
        int BC = BCType[bc];         //translate mesh's boundary flag
        EToB[e*mesh.Nfaces+f] = BC;  //record it
      }
    }
  }
  o_EToB = platform.malloc<int>(EToB);

  //translate the mesh's node-wise bc flag
  Nmasked = 0;
  mapB.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np, 0);

  for (int n=0;n<mesh.Nelements*mesh.Np;n++) {
    int bc = gmshMapB[n];
    if (bc>0) {
      int BC = BCType[bc];     //translate mesh's boundary flag
      mapB[n] = BC;  //record it
      if (mapB[n] > 0) Nmasked++;   //Dirichlet boundary
    }
  }
  o_mapB = platform.malloc<int>(mapB);

  maskIds.malloc(Nmasked);
  Nmasked =0; //reset
  for (dlong n=0;n<mesh.Nelements*mesh.Np;n++) {
    if (mapB[n] > 0) maskIds[Nmasked++] = n;
  }
  o_maskIds = platform.malloc<int>(maskIds);

  //make a masked version of the global id numbering
  maskedGlobalIds.malloc(mesh.Nelements*mesh.Np);
  maskedGlobalIds.copyFrom(mesh.globalIds);
  for (dlong n=0;n<Nmasked;n++) {
    maskedGlobalIds[maskIds[n]] = 0;
  }

  //use the masked ids to make another gs handle (signed so the gather is defined)
  bool verbose = settings.compareSetting("VERBOSE", "TRUE") ? true : false;
  bool unique = true; //flag a unique node in every gather node
  ogsMasked.Setup(mesh.Nelements*mesh.Np, maskedGlobalIds,
                  mesh.comm, ogs::Signed, ogs::Auto,
                  unique, verbose, platform);

  /* use the masked gs handle to define a global ordering */
  dlong Ntotal  = mesh.Np*mesh.Nelements; // number of degrees of freedom on this rank (before gathering)
  hlong Ngather = ogsMasked.Ngather;     // number of degrees of freedom on this rank (after gathering)

  // build inverse degree vectors
  // used for the weight in linear solvers (used in C0)
  weight.malloc(Ntotal, (pfloat)1.0);

  weightG.malloc(Ngather);
  ogsMasked.Gather(weightG, weight, 1, ogs::Add, ogs::Trans);

  for(dlong n=0;n<ogsMasked.Ngather;++n) {
    if (weightG[n]>0.0) weightG[n] = 1./weightG[n];
  }

  ogsMasked.Scatter(weight, weightG, 1, ogs::NoTrans);

  // TW
  o_weight  = platform.malloc<pfloat>(weight);
  o_weightG = platform.malloc<pfloat>(weightG);

  // create a global numbering system
  memory<hlong> globalIds(Ngather);

  // every gathered degree of freedom has its own global id
  hlong globalOffset=static_cast<hlong>(Ngather);
  comm.Scan(Ngather, globalOffset);
  globalOffset = globalOffset-Ngather;

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<ogsMasked.Ngather;n++) {
    globalIds[n] = n + globalOffset;
  }

  //scatter this numbering to the original nodes
  maskedGlobalNumbering.malloc(Ntotal, -1);
  ogsMasked.Scatter(maskedGlobalNumbering, globalIds, 1, ogs::NoTrans);

  /* Build halo exchange for gathered ordering */
  gHalo.SetupFromGather(ogsMasked);

  GlobalToLocal.malloc(mesh.Nelements*mesh.Np);
  ogsMasked.SetupGlobalToLocalMapping(GlobalToLocal);

  o_GlobalToLocal = platform.malloc<dlong>(GlobalToLocal);
}
