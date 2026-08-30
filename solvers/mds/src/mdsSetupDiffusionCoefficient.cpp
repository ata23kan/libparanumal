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

#include <algorithm>
#include <cmath>
#include <limits>


void mds_t::SetupDiffusionCoefficient() {

  switch(mesh.elementType){
  case Mesh::TRIANGLES:
    SetupDiffusionCoefficientTri2D(); break;
  case Mesh::QUADRILATERALS:
  {
    LIBP_FORCE_ABORT("Quadrilaterals are not supported yet for mesh deformation!");
    // if(mesh.dim==2)
    //   BuildOperatorMatrixContinuousQuad2D(A);
    // else
    //   BuildOperatorMatrixContinuousQuad3D(A);

    // break;
  }
  case Mesh::TETRAHEDRA:
    LIBP_FORCE_ABORT("Tetrahedrons are not supported yet for mesh deformation!");
  case Mesh::HEXAHEDRA:
    LIBP_FORCE_ABORT("Hexahedral elements are not supported yet for mesh deformation!");
    // BuildOperatorMatrixContinuousHex3D(A); break;
  }
}

namespace {

dfloat pointSegmentDistanceSquared(const dfloat x,
																	 const dfloat y,
																	 const dfloat x0,
																	 const dfloat y0,
																	 const dfloat x1,
																	 const dfloat y1) {
	const dfloat dx = x1-x0;
	const dfloat dy = y1-y0;
	const dfloat lengthSquared = dx*dx + dy*dy;

	if (lengthSquared == 0.0) {
		const dfloat px = x-x0;
		const dfloat py = y-y0;
		return px*px + py*py;
	}

	const dfloat projection = ((x-x0)*dx + (y-y0)*dy)/lengthSquared;
	const dfloat t = std::max((dfloat) 0.0, std::min((dfloat) 1.0, projection));
	const dfloat px = x - (x0 + t*dx);
	const dfloat py = y - (y0 + t*dy);
	return px*px + py*py;
}

} // namespace

void mds_t::SetupDiffusionCoefficientTri2D(){
	LIBP_ABORT("Variable diffusion coefficients currently require a 2D triangular mesh",
						 mesh.dim != 2 || mesh.elementType != Mesh::TRIANGLES);

	dlong localSegmentCount = 0;
	for (dlong e=0; e<mesh.Nelements; ++e) {
		for (int f=0; f<mesh.Nfaces; ++f) {
			if (EToB[e*mesh.Nfaces + f] > 0) ++localSegmentCount;
		}
	}

	memory<dfloat> localSegments(4*localSegmentCount);
	dlong segment = 0;
	for (dlong e=0; e<mesh.Nelements; ++e) {
		for (int f=0; f<mesh.Nfaces; ++f) {
			if (EToB[e*mesh.Nfaces + f] <= 0) continue;

			const int v0 = mesh.faceVertices[f*mesh.NfaceVertices + 0];
			const int v1 = mesh.faceVertices[f*mesh.NfaceVertices + 1];
			const dlong eOffset = e*mesh.Nverts;
			const dlong offset = 4*segment++;

			localSegments[offset + 0] = mesh.EX[eOffset + v0];
			localSegments[offset + 1] = mesh.EY[eOffset + v0];
			localSegments[offset + 2] = mesh.EX[eOffset + v1];
			localSegments[offset + 3] = mesh.EY[eOffset + v1];
		}
	}

	memory<int> segmentCounts(comm.size());
	comm.Allgather(static_cast<int>(localSegmentCount), segmentCounts);

	memory<int> segmentOffsets(comm.size());
	int globalSegmentCount = 0;
	for (int rank=0; rank<comm.size(); ++rank) {
		segmentOffsets[rank] = 4*globalSegmentCount;
		globalSegmentCount += segmentCounts[rank];
		segmentCounts[rank] *= 4;
	}
	LIBP_ABORT("Cannot construct a diffusion coefficient without boundary faces",
						 globalSegmentCount == 0);

	memory<dfloat> boundarySegments(4*globalSegmentCount);
	comm.Allgatherv(localSegments, static_cast<int>(4*localSegmentCount),
									boundarySegments, segmentCounts, segmentOffsets);

	dfloat meanSegmentLength = 0.0;
	for (int b=0; b<globalSegmentCount; ++b) {
		const dfloat dx = boundarySegments[4*b + 2] - boundarySegments[4*b + 0];
		const dfloat dy = boundarySegments[4*b + 3] - boundarySegments[4*b + 1];
		meanSegmentLength += std::sqrt(dx*dx + dy*dy);
	}
	meanSegmentLength /= globalSegmentCount;
	const dfloat gammaMin = 1.0E-4;
	const dfloat gammaMax = 1.0;
	const dfloat epsilon = meanSegmentLength;

	gamma.malloc(mesh.Nelements);
	for (dlong e=0; e<mesh.Nelements; ++e) {
		dfloat gammaElement = 0.0;
		for (int n=0; n<mesh.Np; ++n) {
			const dlong id = n + e*mesh.Np;
			dfloat distanceSquared = std::numeric_limits<dfloat>::max();
			for (int b=0; b<globalSegmentCount; ++b) {
				distanceSquared = std::min(distanceSquared,
						pointSegmentDistanceSquared(mesh.x[id], mesh.y[id],
																				boundarySegments[4*b + 0], boundarySegments[4*b + 1],
																				boundarySegments[4*b + 2], boundarySegments[4*b + 3]));
			}
			const dfloat distance = std::sqrt(distanceSquared);
			// gammaElement += 1.0/((distance + epsilon)*(distance + epsilon));
			const dfloat normalizedDistance = distance/meanSegmentLength;
			const dfloat proximity = 1.0/((normalizedDistance + 1.0)*(normalizedDistance + 1.0));

			gammaElement += proximity;		
		}
		// gamma[e] = gammaElement/mesh.Np;
		gamma[e] = gammaMin + (gammaMax - gammaMin)*gammaElement/mesh.Np;
	}

	dfloat gammaMinimum = *std::min_element(gamma.ptr(),
																			gamma.ptr() + gamma.length());
	dfloat gammaMaximum = *std::max_element(gamma.ptr(),
																			gamma.ptr() + gamma.length());
	// printf("Min Gamma: %e\n", gammaMinimum);
	// printf("Max Gamma: %e\n", gammaMaximum);
	// std::exit(EXIT_SUCCESS);

	o_gamma = platform.malloc<dfloat>(gamma);

}