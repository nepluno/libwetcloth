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

#ifndef __PRESSURE_H__
#define __PRESSURE_H__

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
