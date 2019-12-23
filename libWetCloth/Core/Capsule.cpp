//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
