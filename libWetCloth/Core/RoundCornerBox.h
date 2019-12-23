//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
