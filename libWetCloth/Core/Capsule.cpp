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

#include "Capsule.h"
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

Capsule::Capsule(int N, const scalar& radius, const scalar& halfheight)
{
	vertices.clear();
	indices.clear();
	
	vertices.resize((N + 1) * (N / 2 + 2));
	
	int N_4 = N / 4;
	int N_2 = N / 2;
	
	int index = 0;
	for(int j = 0; j <= N_4; ++j)
	{
		for(int i = 0; i <= N; ++i)
		{
			scalar theta = (scalar) i * M_PI * 2.0 / (scalar) N;
			scalar phi = -M_PI / 2.0 + M_PI * (scalar) j / (scalar) N_2;
			vertices[index](0) = radius * sin(phi) - halfheight;
			vertices[index](1) = radius * cos(phi) * sin(theta);
			vertices[index](2) = radius * cos(phi) * cos(theta);
			++index;
		}
	}
	
	for(int j = N_4; j <= N_2; ++j)
	{
		for(int i = 0; i <= N; ++i) {
			scalar theta = (scalar) i * M_PI * 2.0 / (scalar) N;
			scalar phi = -M_PI / 2.0 + M_PI * (scalar) j / (scalar) N_2;
			vertices[index](0) = radius * sin(phi) + halfheight;
			vertices[index](1) = radius * cos(phi) * sin(theta);
			vertices[index](2) = radius * cos(phi) * cos(theta);
			++index;
		}
	}
	
	assert(index == (int) vertices.size());
	
	indices.resize((N_2 + 1) * N * 2);
	index = 0;
	for(int j = 0; j <= N_2; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			int i1 = j * (N+1) + i;
			int i2 = j * (N+1) + (i + 1);
			int i3 = (j + 1) * (N + 1) + (i + 1);
			int i4 = (j + 1) * (N + 1) + i;
			indices[index++] = Vector3i(i1, i2, i3);
			indices[index++] = Vector3i(i1, i3, i4);
		}
	}
}
