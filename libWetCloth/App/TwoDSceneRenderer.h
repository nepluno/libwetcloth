//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TWO_D_SCENE_RENDERER_H
#define TWO_D_SCENE_RENDERER_H

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
