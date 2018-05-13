//
// This file is part of the libWetCloth open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// We would like to hear from you if you appreciate this work.
//
// The Clear BSD License
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef SOLID_MESH_H__
#define SOLID_MESH_H__

/* Capsule generation algorithm.
 * Adapted from Paul Bourke's C implementation found here:
 * http://paulbourke.net/geometry/capsule/
 */

#include "MathDefs.h"
#include <vector>
#include <unordered_map>
#include <fstream>

class SolidMesh
{
protected:
	std::vector<Vector3s> vertices;
	std::vector<Vector3i> indices;
public:
	const std::vector<Vector3s>& getVertices() const { return vertices; }
	const std::vector<Vector3i>& getIndices() const { return indices; }
	
	void WriteOBJ(const std::string& szfn) const
	{
		std::ofstream fout(szfn);
		WriteOBJ(fout);
		fout.close();
	}
	void WriteOBJ(std::ostream& str) const
	{
		// vertices
		for(const Vector3s& v : vertices)
		{
			str << "v " << v.transpose() << std::endl;
		}
		
		for(const Vector3i& f : indices)
		{
			str << "f " << f(0)+1 << "//" << f(0)+1 << " " << f(1)+1 << "//" << f(1)+1 << " " << f(2)+1 << "//" << f(2)+1 << std::endl;
		}
	}
};

#endif
