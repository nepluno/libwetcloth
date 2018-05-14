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

#ifndef __TWO_D_SCENE_RENDERER_H__
#define __TWO_D_SCENE_RENDERER_H__

#include <Eigen/StdVector>

#include <iostream>

#include "TwoDScene.h"
#include "MathUtilities.h"
#include "RenderingUtilities.h"

// TODO: Get display controller out of here
// TODO: Make particle system and rigid body renderers that inherit from this

class TwoDimensionalDisplayController;

struct RenderInfo
{
    bool render_particles;
    bool render_vertices;
    bool render_gauss;
    
    bool render_particle_velocity;
    bool render_vertice_velocity;
    bool render_gauss_velocity;
    scalar render_velocity_length;
    scalar render_deformation_gradient_length;
    
    bool render_cloth;
    bool render_yarn;
    bool render_levelset;
    bool render_spring;
    bool render_cohesion;

    
    bool render_buckets;

    enum NODE_VIS
    {
      RI_NV_NONE, RI_NV_CONSTANT, RI_NV_SOLID_PHI
    };
    NODE_VIS render_nodes;
    
    enum FACE_CENTER_VIS
    {
        RI_FCV_NONE, RI_FCV_CONSTANT, RI_FCV_SOLID_VOL, RI_FCV_LIQUID_VOL
    };
    FACE_CENTER_VIS render_face_centers;
    
    enum EDGE_CENTER_VIS
    {
        RI_ECV_NONE, RI_ECV_CONSTANT
    };
    EDGE_CENTER_VIS render_edge_centers;
    
    enum CELL_CENTER_VIS
    {
        RI_CCV_NONE, RI_CCV_CONSTANT, RI_CCV_LIQUID_PHI
    };
    CELL_CENTER_VIS render_cell_centers;
};

class TwoDSceneRenderer
{
public:
	
	// TODO: Gut this method
	TwoDSceneRenderer( const TwoDScene& scene );
	
	void updateParticleSimulationState( const TwoDScene& scene );
	void renderParticleSimulation( const TwoDScene& scene, const scalar& dt );

	// Returns a reference to the vector containing particle colors
	std::vector<renderingutils::Color>& getParticleColors();
	const std::vector<renderingutils::Color>& getParticleColors() const;
    
    RenderInfo& getRenderInfo();
    const RenderInfo& getRenderInfo() const;
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
private:
    RenderInfo m_info;
	
	std::vector< Vector3s > m_color_buffer;
	std::vector< Vector3s > m_group_colors;

};

#endif
