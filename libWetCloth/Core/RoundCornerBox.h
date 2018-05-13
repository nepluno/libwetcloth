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

#ifndef RoundCornerBox_hpp
#define RoundCornerBox_hpp

#include "MathDefs.h"
#include "SolidMesh.h"
#include <vector>
#include <string>
#include <fstream>

class RoundCornerBox : public SolidMesh
{
protected:
	std::vector<int> m_index_to_verts;
	double m_radius;
	int m_N_edge;
	inline void AddVertex(int i, int j, int k, const Vector3s& pos, const Vector3s& base_pos)
	{
		int pidx = k * m_N_edge * m_N_edge + j * m_N_edge + i;
		if(m_index_to_verts[pidx] < 0) {
			int next_idx = (int) vertices.size();
			m_index_to_verts[pidx] = next_idx;
			
			Vector3s dir = pos - base_pos;
			if(dir.norm() > 0.0) {
				dir.normalize();
				vertices.push_back(base_pos + dir * m_radius);
			} else {
				vertices.push_back(pos);
			}
		}
	}
	inline int translateIndices(int i, int j, int k)
	{
		int pidx = k * m_N_edge * m_N_edge + j * m_N_edge + i;
		return m_index_to_verts[pidx];
	}
	inline void AddFace(int i, int j, int k, bool inversed)
	{
		if(inversed)
		{
			indices.push_back(Vector3i(i, k, j));
		} else {
			indices.push_back(Vector3i(i, j, k));
		}
	}
	
public:
	RoundCornerBox(int N, const Vector3s& b, const double& radius);
};


#endif /* CornerBox_hpp */
