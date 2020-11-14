//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef VISCOSITY_H
#define VISCOSITY_H

#include <Eigen/Core>

#include "MathDefs.h"
#include "Array3.h"
#include "PCGSolver/SparseMatrix.h"

class TwoDScene;

using namespace robertbridson;

namespace viscosity {

void constructViscosityMatrixRHS(
    const TwoDScene& scene, std::vector<VectorXi>& node_global_indices_x,
    std::vector<VectorXi>& node_global_indices_y,
    std::vector<VectorXi>& node_global_indices_z,
    std::vector<Vector2i>& effective_node_indices_x,
    std::vector<Vector2i>& effective_node_indices_y,
    std::vector<Vector2i>& effective_node_indices_z,
    const std::vector<VectorXs>& node_vel_src_x,
    const std::vector<VectorXs>& node_vel_src_y,
    const std::vector<VectorXs>& node_vel_src_z, SparseMatrix<scalar>& matrix,
    std::vector<scalar>& rhs, int& offset_nodes_x, int& offset_nodes_y,
    int& offset_nodes_z, const scalar& dt);

void updateViscosityRHS(const TwoDScene& scene,
                        const std::vector<VectorXi>& node_global_indices_x,
                        const std::vector<VectorXi>& node_global_indices_y,
                        const std::vector<VectorXi>& node_global_indices_z,
                        const std::vector<Vector2i>& effective_node_indices_x,
                        const std::vector<Vector2i>& effective_node_indices_y,
                        const std::vector<Vector2i>& effective_node_indices_z,
                        const std::vector<VectorXs>& node_vel_src_x,
                        const std::vector<VectorXs>& node_vel_src_y,
                        const std::vector<VectorXs>& node_vel_src_z,
                        std::vector<scalar>& rhs, int offset_nodes_x,
                        int offset_nodes_y, int offset_nodes_z,
                        const scalar& dt);

void applyNodeViscosityImplicit(
    const TwoDScene& scene, const std::vector<VectorXi>& node_global_indices_x,
    const std::vector<VectorXi>& node_global_indices_y,
    const std::vector<VectorXi>& node_global_indices_z, int offset_nodes_x,
    int offset_nodes_y, int offset_nodes_z, const SparseMatrix<scalar>& matrix,
    const std::vector<scalar>& rhs, std::vector<scalar>& soln,
    std::vector<VectorXs>& node_vel_x, std::vector<VectorXs>& node_vel_y,
    std::vector<VectorXs>& node_vel_z, scalar& residual, int& iter_out,
    const scalar& criterion, int maxiters);

void applyNodeViscosityExplicit(const TwoDScene& scene,
                                const std::vector<VectorXs>& node_vel_src_x,
                                const std::vector<VectorXs>& node_vel_src_y,
                                const std::vector<VectorXs>& node_vel_src_z,
                                std::vector<VectorXs>& node_vel_x,
                                std::vector<VectorXs>& node_vel_y,
                                std::vector<VectorXs>& node_vel_z,
                                const scalar& dt);
};  // namespace viscosity

#endif
