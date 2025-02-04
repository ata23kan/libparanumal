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

void mds_t::Operator(deviceMemory<double> &o_q, deviceMemory<double> &o_Aq){

  deviceMemory<double> o_MM, o_D, o_S, o_LIFT;
  deviceMemory<double> o_wJ, o_ggeo, o_sgeo, o_vgeo;

  if constexpr (std::is_same_v<dfloat,double>) {
    o_MM   = mesh.o_MM;
    o_D    = mesh.o_D;
    o_S    = mesh.o_S;
    o_LIFT = mesh.o_LIFT;

    o_wJ   = mesh.o_wJ;
    o_ggeo = mesh.o_ggeo;
    o_sgeo = mesh.o_sgeo;
    o_vgeo = mesh.o_vgeo;
  } else if (std::is_same_v<pfloat,double>) {
    o_MM   = mesh.o_pfloat_MM;
    o_D    = mesh.o_pfloat_D;
    o_S    = mesh.o_pfloat_S;
    o_LIFT = mesh.o_pfloat_LIFT;

    o_wJ   = mesh.o_pfloat_wJ;
    o_ggeo = mesh.o_pfloat_ggeo;
    o_sgeo = mesh.o_pfloat_sgeo;
    o_vgeo = mesh.o_pfloat_vgeo;
  } else {
    LIBP_FORCE_ABORT("mds_t::Operator called on type double, but double not set in types.h");
  }

  //buffer for local Ax
  deviceMemory<double> o_AqL = platform.reserve<double>(mesh.Np*mesh.Nelements);

  // int mapType = (mesh.elementType==Mesh::HEXAHEDRA &&
  //                mesh.settings.compareSetting("ELEMENT MAP", "TRILINEAR")) ? 1:0;

  // int integrationType = (mesh.elementType==Mesh::HEXAHEDRA &&
  //                        settings.compareSetting("mds INTEGRATION", "CUBATURE")) ? 1:0;

  gHalo.ExchangeStart(o_q, 1);

  if(mesh.NlocalGatherElements/2){
    // if(integrationType==0) { // GLL or non-hex
      // if(mapType==0)
        partialAxKernel(mesh.NlocalGatherElements/2,
                        mesh.o_localGatherElementList,
                        o_GlobalToLocal,
                        o_wJ,
                        o_ggeo,
                        o_D,
                        o_S,
                        o_MM,
                        static_cast<double>(lambda),
                        o_q,
                        o_AqL);
      /* NC: disabling until we re-add treatment of affine elements
      else
        partialAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList,
                        mesh.o_EXYZ, mesh.o_gllzw, mesh.o_D, mesh.o_S, mesh.o_MM, lambda, o_q, o_Aq);
      */
    // } else {
    //   partialCubatureAxKernel(mesh.NlocalGatherElements,
    //                           mesh.o_localGatherElementList,
    //                           mesh.o_cubggeo,
    //                           mesh.o_cubD,
    //                           mesh.o_cubInterpT,
    //                           lambda,
    //                           o_q,
    //                           o_Aq);
    // }
  }

  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, 1);

  if(mesh.NglobalGatherElements) {

    // if(integrationType==0) { // GLL or non-hex
      // if(mapType==0)
        partialAxKernel(mesh.NglobalGatherElements,
                        mesh.o_globalGatherElementList,
                        o_GlobalToLocal,
                        o_wJ,
                        o_ggeo,
                        o_D,
                        o_S,
                        o_MM,
                        static_cast<double>(lambda),
                        o_q,
                        o_AqL);
      /* NC: disabling until we re-add treatment of affine elements
      else
        partialAxKernel(mesh.NglobalGatherElements, mesh.o_globalGatherElementList,
                        mesh.o_EXYZ, mesh.o_gllzw, mesh.o_D, mesh.o_S, mesh.o_MM, lambda, o_q, o_Aq);
      */
    // } else {
    //   partialCubatureAxKernel(mesh.NglobalGatherElements,
    //                           mesh.o_globalGatherElementList,
    //                           mesh.o_cubggeo,
    //                           mesh.o_cubD,
    //                           mesh.o_cubInterpT,
    //                           lambda, o_q, o_Aq);
    // }
  }

  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, 1, ogs::Add, ogs::Trans);

  if((mesh.NlocalGatherElements+1)/2){
    partialAxKernel((mesh.NlocalGatherElements+1)/2,
                    mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
                    o_GlobalToLocal,
                    o_wJ,
                    o_ggeo,
                    o_D,
                    o_S,
                    o_MM,
                    static_cast<double>(lambda),
                    o_q,
                    o_AqL);
  }

  ogsMasked.GatherFinish(o_Aq, o_AqL, 1, ogs::Add, ogs::Trans);
}


void mds_t::Operator(deviceMemory<float> &o_q, deviceMemory<float> &o_Aq){

  deviceMemory<float> o_MM, o_D, o_S, o_LIFT;
  deviceMemory<float> o_wJ, o_ggeo, o_sgeo, o_vgeo;

  if constexpr (std::is_same_v<dfloat,float>) {
    o_MM   = mesh.o_MM;
    o_D    = mesh.o_D;
    o_S    = mesh.o_S;
    o_LIFT = mesh.o_LIFT;

    o_wJ   = mesh.o_wJ;
    o_ggeo = mesh.o_ggeo;
    o_sgeo = mesh.o_sgeo;
    o_vgeo = mesh.o_vgeo;
  } else if (std::is_same_v<pfloat,float>) {
    o_MM   = mesh.o_pfloat_MM;
    o_D    = mesh.o_pfloat_D;
    o_S    = mesh.o_pfloat_S;
    o_LIFT = mesh.o_pfloat_LIFT;

    o_wJ   = mesh.o_pfloat_wJ;
    o_ggeo = mesh.o_pfloat_ggeo;
    o_sgeo = mesh.o_pfloat_sgeo;
    o_vgeo = mesh.o_pfloat_vgeo;
  } else {
    LIBP_FORCE_ABORT("mds_t::Operator called on type float, but float not set in types.h");
  }

  //buffer for local Ax
  deviceMemory<float> o_AqL = platform.reserve<float>(mesh.Np*mesh.Nelements);

  gHalo.ExchangeStart(o_q, 1);

  if(mesh.NlocalGatherElements/2){
    floatPartialAxKernel(mesh.NlocalGatherElements/2,
                         mesh.o_localGatherElementList,
                         o_GlobalToLocal,
                         o_wJ,
                         o_ggeo,
                         o_D,
                         o_S,
                         o_MM,
                         static_cast<float>(lambda),
                         o_q,
                         o_AqL);
  }

  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, 1);

  if(mesh.NglobalGatherElements) {
    floatPartialAxKernel(mesh.NglobalGatherElements,
                         mesh.o_globalGatherElementList,
                         o_GlobalToLocal,
                         o_wJ,
                         o_ggeo,
                         o_D,
                         o_S,
                         o_MM,
                         static_cast<float>(lambda),
                         o_q,
                         o_AqL);
  }

  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, 1, ogs::Add, ogs::Trans);

  if((mesh.NlocalGatherElements+1)/2){
    floatPartialAxKernel((mesh.NlocalGatherElements+1)/2,
                         mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
                         o_GlobalToLocal,
                         o_wJ,
                         o_ggeo,
                         o_D,
                         o_S,
                         o_MM,
                         static_cast<float>(lambda),
                         o_q,
                         o_AqL);
  }

  ogsMasked.GatherFinish(o_Aq, o_AqL, 1, ogs::Add, ogs::Trans);
}