//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
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

#include "Icosphere.h"

Icosphere::Icosphere(int recursionLevel, const scalar& radius)
: index(0)
{
	middlePointIndexCache.clear();
	vertices.clear();
	indices.clear();
	index = 0;
	
	const scalar t = (1.0 + sqrt(5.0)) / 2.0;
	
	addVertex(Vector3s(-1,  t,  0));
	addVertex(Vector3s( 1,  t,  0));
	addVertex(Vector3s(-1, -t,  0));
	addVertex(Vector3s( 1, -t,  0));
	
	addVertex(Vector3s( 0, -1,  t));
	addVertex(Vector3s( 0,  1,  t));
	addVertex(Vector3s( 0, -1, -t));
	addVertex(Vector3s( 0,  1, -t));
	
	addVertex(Vector3s( t,  0, -1));
	addVertex(Vector3s( t,  0,  1));
	addVertex(Vector3s(-t,  0, -1));
	addVertex(Vector3s(-t,  0,  1));
	
	// 5 faces around point 0
	indices.push_back(Vector3i(0, 11, 5));
	indices.push_back(Vector3i(0, 5, 1));
	indices.push_back(Vector3i(0, 1, 7));
	indices.push_back(Vector3i(0, 7, 10));
	indices.push_back(Vector3i(0, 10, 11));
	
	indices.push_back(Vector3i(1, 5, 9));
	indices.push_back(Vector3i(5, 11, 4));
	indices.push_back(Vector3i(11, 10, 2));
	indices.push_back(Vector3i(10, 7, 6));
	indices.push_back(Vector3i(7, 1, 8));
	
	indices.push_back(Vector3i(3, 9, 4));
	indices.push_back(Vector3i(3, 4, 2));
	indices.push_back(Vector3i(3, 2, 6));
	indices.push_back(Vector3i(3, 6, 8));
	indices.push_back(Vector3i(3, 8, 9));
	
	indices.push_back(Vector3i(4, 9, 5));
	indices.push_back(Vector3i(2, 4, 11));
	indices.push_back(Vector3i(6, 2, 10));
	indices.push_back(Vector3i(8, 6, 7));
	indices.push_back(Vector3i(9, 8, 1));
	
	// refine triangles
	for (int i = 0; i < recursionLevel; i++)
	{
		std::vector<Vector3i> indices2;
		for(const Vector3i& tri : indices)
		{
			// replace triangle by 4 triangles
			int a = getMiddlePoint(tri(0), tri(1));
			int b = getMiddlePoint(tri(1), tri(2));
			int c = getMiddlePoint(tri(2), tri(0));
			
			indices2.push_back(Vector3i(tri(0), a, c));
			indices2.push_back(Vector3i(tri(1), b, a));
			indices2.push_back(Vector3i(tri(2), c, b));
			indices2.push_back(Vector3i(a, b, c));
		}
		indices = indices2;
	}
	
	for(auto& v : vertices) {
		v *= radius;
	}
}

int Icosphere::addVertex(const Vector3s& p)
{
	vertices.push_back(p.normalized());
	return index++;
}

int Icosphere::getMiddlePoint(int p1, int p2)
{
	int smallerIndex = std::min(p1, p2);
	int greaterIndex = std::max(p1, p2);
	
	uint64 key = ((uint64)(smallerIndex) << 32UL) | (uint64) greaterIndex;
	auto itr = middlePointIndexCache.find(key);
	if(itr != middlePointIndexCache.end())
	{
		return itr->second;
	}
	
	Vector3s middle = (vertices[p1] + vertices[p2]) * 0.5;
	int i = addVertex(middle);
	
	middlePointIndexCache[key] = i;
	return i;
}

