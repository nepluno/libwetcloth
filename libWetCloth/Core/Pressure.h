//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef PRESSURE_H
#define PRESSURE_H

#include <Eigen/Core>

#include "MathDefs.h"
#include "array3.h"
#include "pcgsolver/sparse_matrix.h"

class TwoDScene;

namespace pressure {
void constructNodeIncompressibleCondition(const TwoDScene& scene,
    std::vector< VectorXs >& node_ic,
    const std::vector< VectorXs >& node_fluid_vel_x,
    const std::vector< VectorXs >& node_fluid_vel_y,
    const std::vector< VectorXs >& node_fluid_vel_z,
    const std::vector< VectorXs >& node_elasto_vel_x,
    const std::vector< VectorXs >& node_elasto_vel_y,
    const std::vector< VectorXs >& node_elasto_vel_z);

void multiplyPressureMatrix( const TwoDScene& scene, const std::vector< VectorXs >& node_vec, std::vector< VectorXs >& out_node_vec, const std::vector< VectorXs >& node_inv_mdv_x, const std::vector< VectorXs >& node_inv_mdv_y, const std::vector< VectorXs >& node_inv_mdv_z, const std::vector< VectorXs >& node_inv_mdvs_x, const std::vector< VectorXs >& node_inv_mdvs_y, const std::vector< VectorXs >& node_inv_mdvs_z, const scalar& dt );
void solveNodePressure( const TwoDScene& scene,
                        std::vector< VectorXs >& pressure,
                        std::vector<double>& rhs,
                        robertbridson::SparseMatrix<scalar>& matrix,
                        std::vector< VectorXi >& node_global_indices,
                        const std::vector< VectorXs >& node_psi_fs_x,
                        const std::vector< VectorXs >& node_psi_fs_y,
                        const std::vector< VectorXs >& node_psi_fs_z,
                        const std::vector< VectorXs >& node_psi_sf_x,
                        const std::vector< VectorXs >& node_psi_sf_y,
                        const std::vector< VectorXs >& node_psi_sf_z,
                        const std::vector< VectorXs >& node_fluid_vel_x,
                        const std::vector< VectorXs >& node_fluid_vel_y,
                        const std::vector< VectorXs >& node_fluid_vel_z,
                        const std::vector< VectorXs >& node_elasto_vel_x,
                        const std::vector< VectorXs >& node_elasto_vel_y,
                        const std::vector< VectorXs >& node_elasto_vel_z,
                        const std::vector< VectorXs >& node_inv_C_x,
                        const std::vector< VectorXs >& node_inv_C_y,
                        const std::vector< VectorXs >& node_inv_C_z,
                        const std::vector< VectorXs >& node_inv_Cs_x,
                        const std::vector< VectorXs >& node_inv_Cs_y,
                        const std::vector< VectorXs >& node_inv_Cs_z,
                        const std::vector< VectorXs >& node_mfhdvm_hdvm_x, // (M_f+hDVm)^{-1}hDVm
                        const std::vector< VectorXs >& node_mfhdvm_hdvm_y,
                        const std::vector< VectorXs >& node_mfhdvm_hdvm_z,
                        const std::vector< VectorXs >& node_mshdvm_hdvm_x, // (M_s+hDVm)^{-1}hDVm
                        const std::vector< VectorXs >& node_mshdvm_hdvm_y,
                        const std::vector< VectorXs >& node_mshdvm_hdvm_z,
                        const scalar& dt,
                        const scalar& criterion,
                        int maxiters );


void constructJacobiPreconditioner( const TwoDScene& scene, std::vector< VectorXs >& out_node_vec, const std::vector< VectorXs >& node_inv_mdv_x, const std::vector< VectorXs >& node_inv_mdv_y, const std::vector< VectorXs >& node_inv_mdv_z, const std::vector< VectorXs >& node_inv_mdvs_x, const std::vector< VectorXs >& node_inv_mdvs_y, const std::vector< VectorXs >& node_inv_mdvs_z, const scalar& dt );
void allocateNodes( const TwoDScene& scene, std::vector< VectorXs >& out_node_vec );
void localSolveJacobi( const std::vector< VectorXs >& node_vec, std::vector< VectorXs >& out_node_vec, const std::vector< VectorXs >& jacobi_precond);

void applySurfTensionFluid( TwoDScene& scene, const std::vector< VectorXs >& surf_tension_vec, std::vector< VectorXs >& rhs_vec_x, std::vector< VectorXs >& rhs_vec_y, std::vector< VectorXs >& rhs_vec_z, const std::vector< VectorXs >& node_vol_x, const std::vector< VectorXs >& node_vol_y, const std::vector< VectorXs >& node_vol_z );

void applyPressureGradsFluid( TwoDScene& scene,
                              const std::vector< VectorXs >& pressure_vec,
                              std::vector< VectorXs >& rhs_vec_x,
                              std::vector< VectorXs >& rhs_vec_y,
                              std::vector< VectorXs >& rhs_vec_z,
                              const std::vector< VectorXs >& node_inv_mdv_x,
                              const std::vector< VectorXs >& node_inv_mdv_y,
                              const std::vector< VectorXs >& node_inv_mdv_z,
                              const scalar& dt );

void applyPressureGradsElastoRHS(
  TwoDScene& scene,
  const std::vector< VectorXs >& pressure_vec,
  std::vector< VectorXs >& rhs_vec_x,
  std::vector< VectorXs >& rhs_vec_y,
  std::vector< VectorXs >& rhs_vec_z,
  const std::vector< VectorXs >& node_mfhdvm_hdvm_x, // (M_f+hDVm)^{-1}hDVm
  const std::vector< VectorXs >& node_mfhdvm_hdvm_y,
  const std::vector< VectorXs >& node_mfhdvm_hdvm_z,
  const scalar& dt );

void applyPressureGradsFluidRHS(
  TwoDScene& scene,
  const std::vector< VectorXs >& pressure_vec,
  std::vector< VectorXs >& rhs_vec_x,
  std::vector< VectorXs >& rhs_vec_y,
  std::vector< VectorXs >& rhs_vec_z,
  const std::vector< VectorXs >& node_mshdvm_hdvm_x, // (M_s+hDVm)^{-1}hDVm
  const std::vector< VectorXs >& node_mshdvm_hdvm_y,
  const std::vector< VectorXs >& node_mshdvm_hdvm_z,
  const scalar& dt );


void extrapolate( const TwoDScene& scene, std::vector< VectorXuc >& node_valid, std::vector< VectorXs >& node_vel);

void computePorePressureGrads(const TwoDScene& scene, std::vector< VectorXs >& rhs_vec_x, std::vector< VectorXs >& rhs_vec_y, std::vector< VectorXs >& rhs_vec_z, const std::vector< VectorXs >& node_vol_x, const std::vector< VectorXs >& node_vol_y, const std::vector< VectorXs >& node_vol_z, const scalar& dt);
};

#endif
