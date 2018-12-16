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
