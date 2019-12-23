//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "RoundCylinder.h"
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>


RoundCylinder::RoundCylinder(int N, int M, const double& ra, const double& rb, const double& h)
{
	vertices.clear();
	indices.clear();
	
    const int num_verts = 2 + N * (M + 1) * 2;
	vertices.resize(num_verts);

    // top and bottom centers
    vertices[0] = Vector3s(0, h + rb, 0);
    vertices[1] = Vector3s(0, -h - rb, 0);
    
    const scalar dtheta = M_PI * 2.0 / (scalar) N;
    const scalar dphi = M_PI_2 / (scalar) M;
    
    for(int i = 0; i < N; ++i)
    {
        double theta = dtheta * (scalar) i;
        double ct = cos(theta);
        double st = sin(theta);
        
        for(int j = 0; j <= M; ++j)
        {
            double phi = dphi * (scalar) j;
            double cp = cos(phi);
            double sp = sin(phi);
            
            const int index_u = 2 + i * (M + 1) + j;
            const int index_d = index_u + N * (M + 1);
            
            vertices[index_u](0) = vertices[index_d](0) = (ra + rb * cp) * ct;
            vertices[index_u](1) = h + rb * sp;
            vertices[index_d](1) = -(h + rb * sp);
            vertices[index_u](2) = vertices[index_d](2) = (ra + rb * cp) * st;
        }
    }
    
    indices.resize(N * (M + 1) * 4);
    
    int base_idx = 0;
    // top and bottom cap
    for(int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        indices[i] = Vector3i(0, 2 + j * (M + 1) + M, 2 + i * (M + 1) + M);
        indices[i + N] = Vector3i(1, 2 + (i + N) * (M + 1) + M, 2 + (j + N) * (M + 1) + M);
    }
    
    base_idx += 2 * N;
    
    // side pads
    for(int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        indices[base_idx + i * 2 + 0] = Vector3i(2 + j * (M + 1), 2 + (i + N) * (M + 1), 2 + i * (M + 1));
        indices[base_idx + i * 2 + 1] = Vector3i(2 + j * (M + 1), 2 + (j + N) * (M + 1), 2 + (i + N) * (M + 1));
    }
    
    base_idx += 2 * N;
    
    // upper bevel
    for(int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        for(int k = 0; k < M; ++k) {
            indices[base_idx + (i * M + k) * 2 + 0] = Vector3i(2 + j * (M + 1) + (k + 1), 2 + i * (M + 1) + k, 2 + i * (M + 1) + (k + 1));
            indices[base_idx + (i * M + k) * 2 + 1] = Vector3i(2 + j * (M + 1) + (k + 1), 2 + j * (M + 1) + k, 2 + i * (M + 1) + k);
        }
    }
    
    base_idx += M * N * 2;
    
    // lower bevel
    for(int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        for(int k = 0; k < M; ++k) {
            indices[base_idx + (i * M + k) * 2 + 0] = Vector3i(2 + (j + N) * (M + 1) + k, 2 + (i + N) * (M + 1) + (k + 1), 2 + (i + N) * (M + 1) + k);
            indices[base_idx + (i * M + k) * 2 + 1] = Vector3i(2 + (j + N) * (M + 1) + k, 2 + (j + N) * (M + 1) + (k + 1), 2 + (i + N) * (M + 1) + (k + 1));
        }
    }
}
