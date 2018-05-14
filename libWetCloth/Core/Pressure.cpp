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


#include "Pressure.h"
#include "TwoDScene.h"
#include "ThreadUtils.h"
#include "MathUtilities.h"
#include "AlgebraicMultigrid.h"

#include <numeric>

const scalar theta_criterion = 0.01;

//#define CHECK_AMGPCG_RESULT

namespace pressure
{
	
	void computePorePressureGrads(const TwoDScene& scene, std::vector< VectorXs >& rhs_vec_x, std::vector< VectorXs >& rhs_vec_y, std::vector< VectorXs >& rhs_vec_z, const std::vector< VectorXs >& node_vol_x, const std::vector< VectorXs >& node_vol_y, const std::vector< VectorXs >& node_vol_z, const scalar& dt)
	{
		const Sorter& buckets = scene.getParticleBuckets();
		
		const scalar dx = scene.getCellSize();
		const scalar dV = dx * dx * dx;
		
		const std::vector< VectorXs >& node_pore_pressure = scene.getNodePorePressureP();
		const std::vector< VectorXs >& node_sat_x = scene.getNodeSaturationX();
		const std::vector< VectorXs >& node_sat_y = scene.getNodeSaturationY();
		const std::vector< VectorXs >& node_sat_z = scene.getNodeSaturationZ();
		
        const std::vector< VectorXs >& node_weight_x = scene.getNodeSolidWeightX();
        const std::vector< VectorXs >& node_weight_y = scene.getNodeSolidWeightY();
        const std::vector< VectorXs >& node_weight_z = scene.getNodeSolidWeightZ();
        
		const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
		const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
		const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
        
        const std::vector< VectorXs >& orient_x = scene.getNodeOrientationX();
        const std::vector< VectorXs >& orient_y = scene.getNodeOrientationY();
        const std::vector< VectorXs >& orient_z = scene.getNodeOrientationZ();
        
        const std::vector< VectorXs >& shape_factor_x = scene.getNodeShapeFactorX();
        const std::vector< VectorXs >& shape_factor_y = scene.getNodeShapeFactorY();
        const std::vector< VectorXs >& shape_factor_z = scene.getNodeShapeFactorZ();
		
		buckets.for_each_bucket([&] (int bucket_idx) {
			VectorXs& bucket_rhs_x = rhs_vec_x[bucket_idx];
			VectorXs& bucket_rhs_y = rhs_vec_y[bucket_idx];
			VectorXs& bucket_rhs_z = rhs_vec_z[bucket_idx];
			
			const int num_rhs_x = bucket_rhs_x.size();
			const int num_rhs_y = bucket_rhs_y.size();
			const int num_rhs_z = bucket_rhs_z.size();
			
			for(int i = 0; i < num_rhs_x; ++i) {
                if(node_vol_x[bucket_idx](i) > 0.0 && node_weight_x[bucket_idx](i) > 0.0) {
                    const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);
                    const scalar pp0 = (indices[0] < 0 || indices[1] < 0) ? 0.0 : node_pore_pressure[ indices[0] ][ indices[1] ];
                    const scalar pp1 = (indices[2] < 0 || indices[3] < 0) ? 0.0 : node_pore_pressure[ indices[2] ][ indices[3] ];
                    
                    const scalar pressure_grad = (pp1 - pp0) / dx;
                    scalar term = pressure_grad * node_vol_x[bucket_idx](i) * dt * node_sat_x[bucket_idx](i) * node_weight_x[bucket_idx](i) * mathutils::get_rotated_pore_coeff_x(orient_x[bucket_idx].segment<3>(i * 3), shape_factor_x[bucket_idx](i));
                    
                    bucket_rhs_x(i) += term;
                }
			}
			
			for(int i = 0; i < num_rhs_y; ++i) {
                if(node_vol_y[bucket_idx](i) > 0.0 && node_weight_y[bucket_idx](i) > 0.0) {
                    const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);
                    const scalar pp0 = (indices[0] < 0 || indices[1] < 0) ? 0.0 : node_pore_pressure[ indices[0] ][ indices[1] ];
                    const scalar pp1 = (indices[2] < 0 || indices[3] < 0) ? 0.0 : node_pore_pressure[ indices[2] ][ indices[3] ];
                    
                    const scalar pressure_grad = (pp1 - pp0) / dx;
                    scalar term = pressure_grad * node_vol_y[bucket_idx](i) * dt * node_sat_y[bucket_idx](i) * node_weight_y[bucket_idx](i) * mathutils::get_rotated_pore_coeff_y(orient_y[bucket_idx].segment<3>(i * 3), shape_factor_y[bucket_idx](i));
                    
                    bucket_rhs_y(i) += term;
                }
			}
			
			for(int i = 0; i < num_rhs_z; ++i) {
                if(node_vol_z[bucket_idx](i) > 0.0 && node_weight_z[bucket_idx](i) > 0.0) {
                    const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);
                    const scalar pp0 = (indices[0] < 0 || indices[1] < 0) ? 0.0 : node_pore_pressure[ indices[0] ][ indices[1] ];
                    const scalar pp1 = (indices[2] < 0 || indices[3] < 0) ? 0.0 : node_pore_pressure[ indices[2] ][ indices[3] ];
                    
                    const scalar pressure_grad = (pp1 - pp0) / dx;
                    scalar term = pressure_grad * node_vol_z[bucket_idx](i) * dt * node_sat_z[bucket_idx](i) * node_weight_z[bucket_idx](i) * mathutils::get_rotated_pore_coeff_z(orient_z[bucket_idx].segment<3>(i * 3), shape_factor_z[bucket_idx](i));;
                    
                    bucket_rhs_z(i) += term;
                }
			}
		});
	}
	
    void constructDivEquationBiCGSTAB(const TwoDScene& scene,
                                      std::vector< VectorXs >& node_rhs,
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
                                      const scalar& dt)
    {
        const Sorter& buckets = scene.getParticleBuckets();
        
        const int num_buckets = buckets.size();
        
        if((int) node_rhs.size() != num_buckets) node_rhs.resize(num_buckets);
        
        const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
        const std::vector< VectorXi >& pressure_neighbors = scene.getPressureNeighbors();
        const std::vector< VectorXs >& node_fluid_vol_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_fluid_vol_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_fluid_vol_z = scene.getNodeFluidVolZ();
        const std::vector< VectorXs >& node_elasto_vol_x = scene.getNodeVolX();
        const std::vector< VectorXs >& node_elasto_vol_y = scene.getNodeVolY();
        const std::vector< VectorXs >& node_elasto_vol_z = scene.getNodeVolZ();
        const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
        const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
        const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
        const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
        const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
        const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
        const std::vector< VectorXs >& node_solid_vel_x = scene.getNodeSolidVelX();
        const std::vector< VectorXs >& node_solid_vel_y = scene.getNodeSolidVelY();
        const std::vector< VectorXs >& node_solid_vel_z = scene.getNodeSolidVelZ();
        const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
        
        const scalar dx = scene.getCellSize();
        const scalar coeff = dt / (dx * dx);
        
        buckets.for_each_bucket([&] (int bucket_idx) {
            VectorXs& bucket_node_rhs = node_rhs[bucket_idx];
            const VectorXs& bucket_liquid_phi = node_liquid_phi[bucket_idx];
            const VectorXi& bucket_pn = pressure_neighbors[bucket_idx];
            const VectorXs& bucket_pos_p = scene.getNodePosP(bucket_idx);
            
            const int num_rhs = bucket_liquid_phi.size();
            
            bucket_node_rhs.resize(num_rhs);
            
            for(int i = 0; i < num_rhs; ++i)
            {
                bucket_node_rhs(i) = 0.0;
                
                const scalar center_phi = bucket_liquid_phi(i);
                if(center_phi > 0.0) continue;
                const Vector3s center_pos = bucket_pos_p.segment<3>(i * 3);
                
                // left neighbor
                const int bucket_idx_left = bucket_pn[i * 12 + 0];
                const int node_idx_left = bucket_pn[i * 12 + 1];
                if(bucket_idx_left >= 0 && node_idx_left >= 0) {
                    const scalar& w = node_solid_weight_x[bucket_idx_left][node_idx_left];
                    bucket_node_rhs(i) += ((node_fluid_vel_x[bucket_idx_left][node_idx_left] *
                                            node_psi_fs_x[bucket_idx_left][node_idx_left] +
                                            node_elasto_vel_x[bucket_idx_left][node_idx_left] *
                                            node_psi_sf_x[bucket_idx_left][node_idx_left]) * w +
                                           node_solid_vel_x[bucket_idx_left][node_idx_left] * (1.0 - w)) / dx;
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.5 * dx, 0.0, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                    bucket_node_rhs(i) += (scene.interpolateBucketFluidVelX(nnp) * w + scene.interpolateBucketSolidVelX(nnp) * (1.0 - w)) / dx;
                    
                    const scalar left_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( dx, 0., 0. ));
                    
                    if(left_phi < 0.0) {
                        scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                        bucket_node_rhs(i) += term * scene.interpolateBucketPressure(center_pos - Vector3s( dx, 0., 0. ));
                    }
                }
                
                // right neighbor
                const int bucket_idx_right = bucket_pn[i * 12 + 2];
                const int node_idx_right = bucket_pn[i * 12 + 3];
                if(bucket_idx_right >= 0 && node_idx_right >= 0) {
                    const scalar& w = node_solid_weight_x[bucket_idx_right][node_idx_right];
                    bucket_node_rhs(i) -= ((node_fluid_vel_x[bucket_idx_right][node_idx_right] *
                                            node_psi_fs_x[bucket_idx_right][node_idx_right] +
                                            node_elasto_vel_x[bucket_idx_right][node_idx_right] *
                                            node_psi_sf_x[bucket_idx_right][node_idx_right]) * w +
                                           node_solid_vel_x[bucket_idx_right][node_idx_right] * (1.0 - w)) / dx;
                    
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.5 * dx, 0.0, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                    bucket_node_rhs(i) -= (scene.interpolateBucketFluidVelX(nnp) * w + scene.interpolateBucketSolidVelX(nnp) * (1.0 - w)) / dx;
                    
                    const scalar right_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( dx, 0., 0. ));
                    if(right_phi < 0.0) {
                        scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                        bucket_node_rhs(i) += term * scene.interpolateBucketPressure(center_pos + Vector3s( dx, 0., 0. ));
                    }
                }
                
                // bottom neighbor
                const int bucket_idx_bottom = bucket_pn[i * 12 + 4];
                const int node_idx_bottom = bucket_pn[i * 12 + 5];
                if(bucket_idx_bottom >= 0 && node_idx_bottom >= 0) {
                    const scalar& w = node_solid_weight_y[bucket_idx_bottom][node_idx_bottom];
                    bucket_node_rhs(i) += ((node_fluid_vel_y[bucket_idx_bottom][node_idx_bottom] *
                                            node_psi_fs_y[bucket_idx_bottom][node_idx_bottom] +
                                            node_elasto_vel_y[bucket_idx_bottom][node_idx_bottom] *
                                            node_psi_sf_y[bucket_idx_bottom][node_idx_bottom]) * w +
                                           node_solid_vel_y[bucket_idx_bottom][node_idx_bottom] * (1.0 - w)) / dx;
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.0, 0.5 * dx, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                    bucket_node_rhs(i) += (scene.interpolateBucketFluidVelY(nnp) * w + scene.interpolateBucketSolidVelY(nnp) * (1.0 - w)) / dx;
                    
                    const scalar bottom_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., dx, 0. ));
                    
                    if(bottom_phi < 0.0) {
                        scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                        bucket_node_rhs(i) += term * scene.interpolateBucketPressure(center_pos - Vector3s(0, dx, 0));
                    }
                }
                
                // top neighbor
                const int bucket_idx_top = bucket_pn[i * 12 + 6];
                const int node_idx_top = bucket_pn[i * 12 + 7];
                if(bucket_idx_top >= 0 && node_idx_top >= 0) {
                    const scalar& w = node_solid_weight_y[bucket_idx_top][node_idx_top];
                    bucket_node_rhs(i) -= ((node_fluid_vel_y[bucket_idx_top][node_idx_top] *
                                            node_psi_fs_y[bucket_idx_top][node_idx_top] +
                                            node_elasto_vel_y[bucket_idx_top][node_idx_top] *
                                            node_psi_sf_y[bucket_idx_top][node_idx_top]) * w +
                                           node_solid_vel_y[bucket_idx_top][node_idx_top] * (1.0 - w)) / dx;
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.0, 0.5 * dx, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                    bucket_node_rhs(i) -= (scene.interpolateBucketFluidVelY(nnp) * w + scene.interpolateBucketSolidVelY(nnp) * (1.0 - w)) / dx;
                    
                    const scalar top_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., dx, 0. ));
                    
                    if(top_phi < 0.0) {
                        scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                        bucket_node_rhs(i) += term * scene.interpolateBucketPressure(center_pos + Vector3s(0, dx, 0));
                    }
                }
                
                // near neighbor
                const int bucket_idx_near = bucket_pn[i * 12 + 8];
                const int node_idx_near = bucket_pn[i * 12 + 9];
                if(bucket_idx_near >= 0 && node_idx_near >= 0) {
                    const scalar& w = node_solid_weight_z[bucket_idx_near][node_idx_near];
                    bucket_node_rhs(i) += ((node_fluid_vel_z[bucket_idx_near][node_idx_near] *
                                            node_psi_fs_z[bucket_idx_near][node_idx_near] +
                                            node_elasto_vel_z[bucket_idx_near][node_idx_near] *
                                            node_psi_sf_z[bucket_idx_near][node_idx_near]) * w +
                                           node_solid_vel_z[bucket_idx_near][node_idx_near] * (1.0 - w)) / dx;
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.0, 0.0, 0.5 * dx);
                    const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                    bucket_node_rhs(i) += (scene.interpolateBucketFluidVelZ(nnp) * w + scene.interpolateBucketSolidVelZ(nnp) * (1.0 - w)) / dx;
                    
                    const scalar near_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., 0., dx ));
                    
                    if(near_phi < 0.0) {
                        scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                        bucket_node_rhs(i) += term * scene.interpolateBucketPressure(center_pos - Vector3s(0, 0, dx));
                    }
                }
                
                // far neighbor
                const int bucket_idx_far = bucket_pn[i * 12 + 10];
                const int node_idx_far = bucket_pn[i * 12 + 11];
                if(bucket_idx_far >= 0 && node_idx_far >= 0) {
                    const scalar& w = node_solid_weight_z[bucket_idx_far][node_idx_far];
                    bucket_node_rhs(i) -= ((node_fluid_vel_z[bucket_idx_far][node_idx_far] *
                                            node_psi_fs_z[bucket_idx_far][node_idx_far] +
                                            node_elasto_vel_z[bucket_idx_far][node_idx_far] *
                                            node_psi_sf_z[bucket_idx_far][node_idx_far]) * w +
                                           node_solid_vel_z[bucket_idx_far][node_idx_far] * (1.0 - w)) / dx;
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.0, 0.0, 0.5 * dx);
                    const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                    bucket_node_rhs(i) -= (scene.interpolateBucketFluidVelZ(nnp) * w + scene.interpolateBucketSolidVelZ(nnp) * (1.0 - w)) / dx;
                    
                    const scalar far_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., 0., dx ));
                    
                    if(far_phi < 0.0) {
                        scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                        bucket_node_rhs(i) += term * scene.interpolateBucketPressure(center_pos + Vector3s(0, 0, dx));
                    }
                }
            }
        });
    }
    
    void constructNodeIncompressibleCondition(const TwoDScene& scene,
                                              std::vector< VectorXs >& node_ic,
                                              const std::vector< VectorXs >& node_fluid_vel_x,
                                              const std::vector< VectorXs >& node_fluid_vel_y,
                                              const std::vector< VectorXs >& node_fluid_vel_z,
                                              const std::vector< VectorXs >& node_elasto_vel_x,
                                              const std::vector< VectorXs >& node_elasto_vel_y,
                                              const std::vector< VectorXs >& node_elasto_vel_z)
	{
		const Sorter& buckets = scene.getParticleBuckets();
		
		const int num_buckets = buckets.size();
		
		if((int) node_ic.size() != num_buckets) node_ic.resize(num_buckets);
		
		const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
		const std::vector< VectorXi >& pressure_neighbors = scene.getPressureNeighbors();
        const std::vector< VectorXs >& node_fluid_vol_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_fluid_vol_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_fluid_vol_z = scene.getNodeFluidVolZ();
        const std::vector< VectorXs >& node_elasto_vol_x = scene.getNodeVolX();
        const std::vector< VectorXs >& node_elasto_vol_y = scene.getNodeVolY();
        const std::vector< VectorXs >& node_elasto_vol_z = scene.getNodeVolZ();
        const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
        const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
        const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
		const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
		const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
		const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
		const std::vector< VectorXs >& node_solid_vel_x = scene.getNodeSolidVelX();
		const std::vector< VectorXs >& node_solid_vel_y = scene.getNodeSolidVelY();
		const std::vector< VectorXs >& node_solid_vel_z = scene.getNodeSolidVelZ();
        const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
		
		const scalar dx = scene.getCellSize();

		buckets.for_each_bucket([&] (int bucket_idx) {
			VectorXs& bucket_node_rhs = node_ic[bucket_idx];
			const VectorXs& bucket_liquid_phi = node_liquid_phi[bucket_idx];
			const VectorXi& bucket_pn = pressure_neighbors[bucket_idx];
            const VectorXs& bucket_pos_p = scene.getNodePosP(bucket_idx);
			
			const int num_rhs = bucket_liquid_phi.size();
			
			bucket_node_rhs.resize(num_rhs);
			
			for(int i = 0; i < num_rhs; ++i)
			{
				bucket_node_rhs(i) = 0.0;
				
				const scalar center_phi = bucket_liquid_phi(i);
				if(center_phi > 0.0) continue;

				// left neighbor
				const int bucket_idx_left = bucket_pn[i * 12 + 0];
				const int node_idx_left = bucket_pn[i * 12 + 1];
				if(bucket_idx_left >= 0 && node_idx_left >= 0) {
					const scalar& w = node_solid_weight_x[bucket_idx_left][node_idx_left];
					bucket_node_rhs(i) += ((node_fluid_vel_x[bucket_idx_left][node_idx_left] *
											(1.0 - node_psi_x[bucket_idx_left][node_idx_left]) +
											node_elasto_vel_x[bucket_idx_left][node_idx_left] *
											node_psi_x[bucket_idx_left][node_idx_left]) * w +
										   node_solid_vel_x[bucket_idx_left][node_idx_left] * (1.0 - w)) / dx;
                }
				
				// right neighbor
				const int bucket_idx_right = bucket_pn[i * 12 + 2];
				const int node_idx_right = bucket_pn[i * 12 + 3];
				if(bucket_idx_right >= 0 && node_idx_right >= 0) {
					const scalar& w = node_solid_weight_x[bucket_idx_right][node_idx_right];
					bucket_node_rhs(i) -= ((node_fluid_vel_x[bucket_idx_right][node_idx_right] *
											(1.0 - node_psi_x[bucket_idx_right][node_idx_right]) +
											node_elasto_vel_x[bucket_idx_right][node_idx_right] *
											node_psi_x[bucket_idx_right][node_idx_right]) * w +
										   node_solid_vel_x[bucket_idx_right][node_idx_right] * (1.0 - w)) / dx;
                }

				// bottom neighbor
				const int bucket_idx_bottom = bucket_pn[i * 12 + 4];
				const int node_idx_bottom = bucket_pn[i * 12 + 5];
				if(bucket_idx_bottom >= 0 && node_idx_bottom >= 0) {
					const scalar& w = node_solid_weight_y[bucket_idx_bottom][node_idx_bottom];
					bucket_node_rhs(i) += ((node_fluid_vel_y[bucket_idx_bottom][node_idx_bottom] *
											(1.0 - node_psi_y[bucket_idx_bottom][node_idx_bottom]) +
											node_elasto_vel_y[bucket_idx_bottom][node_idx_bottom] *
											node_psi_y[bucket_idx_bottom][node_idx_bottom]) * w +
										   node_solid_vel_y[bucket_idx_bottom][node_idx_bottom] * (1.0 - w)) / dx;
                }
				
				// top neighbor
				const int bucket_idx_top = bucket_pn[i * 12 + 6];
				const int node_idx_top = bucket_pn[i * 12 + 7];
				if(bucket_idx_top >= 0 && node_idx_top >= 0) {
					const scalar& w = node_solid_weight_y[bucket_idx_top][node_idx_top];
					bucket_node_rhs(i) -= ((node_fluid_vel_y[bucket_idx_top][node_idx_top] *
											(1.0 - node_psi_y[bucket_idx_top][node_idx_top]) +
											node_elasto_vel_y[bucket_idx_top][node_idx_top] *
											node_psi_y[bucket_idx_top][node_idx_top]) * w +
										   node_solid_vel_y[bucket_idx_top][node_idx_top] * (1.0 - w)) / dx;
                }

				// near neighbor
				const int bucket_idx_near = bucket_pn[i * 12 + 8];
				const int node_idx_near = bucket_pn[i * 12 + 9];
				if(bucket_idx_near >= 0 && node_idx_near >= 0) {
					const scalar& w = node_solid_weight_z[bucket_idx_near][node_idx_near];
					bucket_node_rhs(i) += ((node_fluid_vel_z[bucket_idx_near][node_idx_near] *
											(1.0 - node_psi_z[bucket_idx_near][node_idx_near]) +
											node_elasto_vel_z[bucket_idx_near][node_idx_near] *
											node_psi_z[bucket_idx_near][node_idx_near]) * w +
										   node_solid_vel_z[bucket_idx_near][node_idx_near] * (1.0 - w)) / dx;
                }
				
				// far neighbor
				const int bucket_idx_far = bucket_pn[i * 12 + 10];
				const int node_idx_far = bucket_pn[i * 12 + 11];
				if(bucket_idx_far >= 0 && node_idx_far >= 0) {
					const scalar& w = node_solid_weight_z[bucket_idx_far][node_idx_far];
					bucket_node_rhs(i) -= ((node_fluid_vel_z[bucket_idx_far][node_idx_far] *
											(1.0 - node_psi_z[bucket_idx_far][node_idx_far]) +
											node_elasto_vel_z[bucket_idx_far][node_idx_far] *
											node_psi_z[bucket_idx_far][node_idx_far]) * w +
										   node_solid_vel_z[bucket_idx_far][node_idx_far] * (1.0 - w)) / dx;
                }
			}
		});
	}
	
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
                           int maxiters )
    {
        const Sorter& buckets = scene.getParticleBuckets();
        const int bucket_num_cell = scene.getDefaultNumNodes();
        const int ni = buckets.ni * bucket_num_cell;
        const int nj = buckets.nj * bucket_num_cell;
        const int nk = buckets.nk * bucket_num_cell;
        
        const std::vector< VectorXi >& node_indices_p = scene.getNodeIndicesP();
        const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
        const std::vector< VectorXi >& pressure_neighbors = scene.getPressureNeighbors();
        const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
        const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
        const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
        const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
        const std::vector< VectorXs >& node_fluid_vol_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_fluid_vol_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_fluid_vol_z = scene.getNodeFluidVolZ();
        const std::vector< VectorXs >& node_elasto_vol_x = scene.getNodeVolX();
        const std::vector< VectorXs >& node_elasto_vol_y = scene.getNodeVolY();
        const std::vector< VectorXs >& node_elasto_vol_z = scene.getNodeVolZ();
        const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
        const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
        const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
        const std::vector< VectorXs >& node_solid_vel_x = scene.getNodeSolidVelX();
        const std::vector< VectorXs >& node_solid_vel_y = scene.getNodeSolidVelY();
        const std::vector< VectorXs >& node_solid_vel_z = scene.getNodeSolidVelZ();
        
        std::vector<int> num_effective_nodes( buckets.size() );
        // assign global indices to nodes
        buckets.for_each_bucket([&] (int bucket_idx) {
            const int num_nodes_p = node_liquid_phi[bucket_idx].size();
            const VectorXi& bucket_pn = pressure_neighbors[bucket_idx];
            
            int k = 0;
            for(int i = 0; i < num_nodes_p; ++i)
            {
                node_global_indices[bucket_idx][i] = -1;
                
                if(node_liquid_phi[bucket_idx][i] < 0) {
                    const int bucket_idx_left = bucket_pn[i * 12 + 0];
                    const int node_idx_left = bucket_pn[i * 12 + 1];
                    const int bucket_idx_right = bucket_pn[i * 12 + 2];
                    const int node_idx_right = bucket_pn[i * 12 + 3];
                    const int bucket_idx_bottom = bucket_pn[i * 12 + 4];
                    const int node_idx_bottom = bucket_pn[i * 12 + 5];
                    const int bucket_idx_top = bucket_pn[i * 12 + 6];
                    const int node_idx_top = bucket_pn[i * 12 + 7];
                    const int bucket_idx_near = bucket_pn[i * 12 + 8];
                    const int node_idx_near = bucket_pn[i * 12 + 9];
                    const int bucket_idx_far = bucket_pn[i * 12 + 10];
                    const int node_idx_far = bucket_pn[i * 12 + 11];
                    
                    if(bucket_idx_left < 0 || node_idx_left < 0 || node_solid_weight_x[bucket_idx_left][node_idx_left] > 0 ||
                       bucket_idx_right < 0 || node_idx_right < 0 || node_solid_weight_x[bucket_idx_right][node_idx_right] > 0 ||
                       bucket_idx_bottom < 0 || node_idx_bottom < 0 || node_solid_weight_y[bucket_idx_bottom][node_idx_bottom] > 0 ||
                       bucket_idx_top < 0 || node_idx_top < 0 || node_solid_weight_y[bucket_idx_top][node_idx_top] > 0 ||
                       bucket_idx_near < 0 || node_idx_near < 0 || node_solid_weight_z[bucket_idx_near][node_idx_near] > 0 ||
                       bucket_idx_far < 0 || node_idx_far < 0 || node_solid_weight_z[bucket_idx_far][node_idx_far] > 0)
                    {
                        node_global_indices[bucket_idx][i] = k++;
                    }
                }
            }
            
            num_effective_nodes[bucket_idx] = k;
        });
        
        std::partial_sum(num_effective_nodes.begin(), num_effective_nodes.end(), num_effective_nodes.begin());
        
        const int total_num_nodes = num_effective_nodes[num_effective_nodes.size() - 1];
        if(total_num_nodes == 0) return;
        
        std::vector< Vector2i > effective_node_indices(total_num_nodes);
        std::vector< Vector3i > dof_ijk(total_num_nodes);
        std::vector< double > result(total_num_nodes);
      
        result.assign(total_num_nodes, 0.0);
        buckets.for_each_bucket([&] (int bucket_idx) {
            pressure[bucket_idx].setZero();
        });
        
        if((int) rhs.size() != total_num_nodes) {
            rhs.resize(total_num_nodes);
            matrix.resize(total_num_nodes);
        }

        matrix.zero();
      
        buckets.for_each_bucket([&] (int bucket_idx) {
            const int num_nodes_p = node_liquid_phi[bucket_idx].size();
            const VectorXi& bucket_node_indices_p = node_indices_p[bucket_idx];
            const Vector3i handle = buckets.bucket_handle(bucket_idx);
            
            for(int i = 0; i < num_nodes_p; ++i)
            {
                if(node_global_indices[bucket_idx][i] >= 0) {
                    if(bucket_idx > 0) {
                        node_global_indices[bucket_idx][i] += num_effective_nodes[bucket_idx - 1];
                    }
                    
                    const int global_idx = node_global_indices[bucket_idx][i];
                    effective_node_indices[global_idx] = Vector2i(bucket_idx, i);
                    
                    const Vector3i& local_handle = bucket_node_indices_p.segment<3>(i * 3);
                    dof_ijk[global_idx] = Vector3i(handle(0) * bucket_num_cell + local_handle(0),
                                                   handle(1) * bucket_num_cell + local_handle(1),
                                                   handle(2) * bucket_num_cell + local_handle(2)
                                                   );
                }
            }
        });
        
        const scalar dx = scene.getCellSize();
        const scalar dV = dx * dx * dx;
        const scalar coeff = dt / (dx * dx);
        
        threadutils::for_each(0, total_num_nodes, [&] (int dof_idx) {
            const Vector2i& dof_loc = effective_node_indices[dof_idx];
            const int bucket_idx = dof_loc[0];
            const int node_idx = dof_loc[1];
            const VectorXs& bucket_pos_p = scene.getNodePosP(bucket_idx);
            
            const scalar center_phi = node_liquid_phi[bucket_idx][node_idx];
            rhs[dof_idx] = 0.0;
            
            const VectorXi& bucket_pn = pressure_neighbors[bucket_idx];
            
            const int bucket_idx_left = bucket_pn[node_idx * 12 + 0];
            const int node_idx_left = bucket_pn[node_idx * 12 + 1];
            
            if(bucket_idx_left >= 0 && node_idx_left >= 0) {
                const scalar& w = node_solid_weight_x[bucket_idx_left][node_idx_left];
                
                rhs[dof_idx] += ((node_fluid_vel_x[bucket_idx_left][node_idx_left] *
                                        node_psi_fs_x[bucket_idx_left][node_idx_left] +
                                        node_elasto_vel_x[bucket_idx_left][node_idx_left] *
                                        node_psi_sf_x[bucket_idx_left][node_idx_left]) * w +
                                       node_solid_vel_x[bucket_idx_left][node_idx_left] * (1.0 - w)) / dx;
                
                if(w > 0.0) {
                    scalar term = w * coeff;
                    
                    term *= (1.0 - node_psi_x[bucket_idx_left][node_idx_left]) * node_inv_C_x[bucket_idx_left][node_idx_left] *
                    (node_fluid_vol_x[bucket_idx_left][node_idx_left] + node_elasto_vol_x[bucket_idx_left][node_idx_left] * node_mshdvm_hdvm_x[bucket_idx_left][node_idx_left])
                    + node_psi_x[bucket_idx_left][node_idx_left] * node_inv_Cs_x[bucket_idx_left][node_idx_left] *
                    (node_elasto_vol_x[bucket_idx_left][node_idx_left] + node_fluid_vol_x[bucket_idx_left][node_idx_left] * node_mfhdvm_hdvm_x[bucket_idx_left][node_idx_left]);
                    
                    const int bucket_pressure_left = node_pressure_index_x[bucket_idx_left][node_idx_left * 4 + 0];
                    const int node_pressure_left = node_pressure_index_x[bucket_idx_left][node_idx_left * 4 + 1];
                    
                    if(bucket_pressure_left >= 0 && node_pressure_left >= 0) {
                        const scalar left_phi = node_liquid_phi[bucket_pressure_left][node_pressure_left];
                        
                        if(left_phi < 0.0) {
                            const int dof_left = node_global_indices[bucket_pressure_left][node_pressure_left];
                            assert(dof_left >= 0);
                            
                            if(dof_left >= 0) {
                                matrix.add_to_element(dof_idx, dof_idx, term);
                                matrix.add_to_element(dof_idx, dof_left, -term);
                            }
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    } else { // coarse grid has zero pressure in fine level
                        const Vector3s& center_pos = bucket_pos_p.segment<3>(node_idx * 3);
                        const scalar left_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( dx, 0., 0. ));
                        
                        if(left_phi < 0.0) {
                            matrix.add_to_element(dof_idx, dof_idx, term);
                            rhs[dof_idx] += term * scene.interpolateBucketPressure(center_pos - Vector3s( dx, 0., 0. ));
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    }
                }
            } else {
                const Vector3s& np = bucket_pos_p.segment<3>(node_idx * 3);
                Vector3s nnp = np - Vector3s(0.5 * dx, 0.0, 0.0);
                const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                rhs[dof_idx] += (scene.interpolateBucketFluidVelX(nnp) * w + scene.interpolateBucketSolidVelX(nnp) * (1.0 - w)) / dx;
                
                if(w > 0.0) {
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar left_phi = scene.interpolateBucketLiquidPhi(np - Vector3s( dx, 0., 0. ));
                    
                    if(left_phi < 0.0) {
                        matrix.add_to_element(dof_idx, dof_idx, term);
                        rhs[dof_idx] += term * scene.interpolateBucketPressure(np - Vector3s( dx, 0., 0. ));
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
                        matrix.add_to_element(dof_idx, dof_idx, term / theta);
                    }
                }
            }
            
            const int bucket_idx_right = bucket_pn[node_idx * 12 + 2];
            const int node_idx_right = bucket_pn[node_idx * 12 + 3];
            
            if(bucket_idx_right >= 0 && node_idx_right >= 0) {
                const scalar& w = node_solid_weight_x[bucket_idx_right][node_idx_right];
                rhs[dof_idx] -= ((node_fluid_vel_x[bucket_idx_right][node_idx_right] *
                                        node_psi_fs_x[bucket_idx_right][node_idx_right] +
                                        node_elasto_vel_x[bucket_idx_right][node_idx_right] *
                                        node_psi_sf_x[bucket_idx_right][node_idx_right]) * w +
                                       node_solid_vel_x[bucket_idx_right][node_idx_right] * (1.0 - w)) / dx;
                
                if(w > 0.0) {
                    scalar term = w * coeff;
                    
                    term *= (1.0 - node_psi_x[bucket_idx_right][node_idx_right]) * node_inv_C_x[bucket_idx_right][node_idx_right] *
                    (node_fluid_vol_x[bucket_idx_right][node_idx_right] + node_elasto_vol_x[bucket_idx_right][node_idx_right] * node_mshdvm_hdvm_x[bucket_idx_right][node_idx_right])
                    + node_psi_x[bucket_idx_right][node_idx_right] * node_inv_Cs_x[bucket_idx_right][node_idx_right] *
                    (node_elasto_vol_x[bucket_idx_right][node_idx_right] + node_fluid_vol_x[bucket_idx_right][node_idx_right] * node_mfhdvm_hdvm_x[bucket_idx_right][node_idx_right]);
                    
                    const int bucket_pressure_right = node_pressure_index_x[bucket_idx_right][node_idx_right * 4 + 2];
                    const int node_pressure_right = node_pressure_index_x[bucket_idx_right][node_idx_right * 4 + 3];
                    
                    if(bucket_pressure_right >= 0 && node_pressure_right >= 0) {
                        const scalar right_phi = node_liquid_phi[bucket_pressure_right][node_pressure_right];
                        
                        if(right_phi < 0.0) {
                            const int dof_right = node_global_indices[bucket_pressure_right][node_pressure_right];
                            assert(dof_right >= 0);
                            
                            if(dof_right >= 0) {
                                matrix.add_to_element(dof_idx, dof_idx, term);
                                matrix.add_to_element(dof_idx, dof_right, -term);
                            }
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    } else { // coarse grid has zero pressure in fine level
                        const Vector3s& center_pos = bucket_pos_p.segment<3>(node_idx * 3);
                        const scalar right_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( dx, 0., 0. ));
                        
                        if(right_phi < 0.0) {
                            matrix.add_to_element(dof_idx, dof_idx, term);
                            rhs[dof_idx] += term * scene.interpolateBucketPressure(center_pos + Vector3s( dx, 0., 0. ));
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    }
                }
            } else {
                const Vector3s& np = bucket_pos_p.segment<3>(node_idx * 3);
                Vector3s nnp = np + Vector3s(0.5 * dx, 0.0, 0.0);
                const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                
                if(w > 0.0) {
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    rhs[dof_idx] -= (scene.interpolateBucketFluidVelX(nnp) * w + scene.interpolateBucketSolidVelX(nnp) * (1.0 - w)) / dx;

                    const scalar right_phi = scene.interpolateBucketLiquidPhi(np + Vector3s( dx, 0., 0. ));
                    
                    if(right_phi < 0.0) {
                        matrix.add_to_element(dof_idx, dof_idx, term);
                        rhs[dof_idx] += term * scene.interpolateBucketPressure(np + Vector3s( dx, 0., 0. ));
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
                        matrix.add_to_element(dof_idx, dof_idx, term / theta);
                    }
                }
            }
            
            const int bucket_idx_bottom = bucket_pn[node_idx * 12 + 4];
            const int node_idx_bottom = bucket_pn[node_idx * 12 + 5];
            
            if(bucket_idx_bottom >= 0 && node_idx_bottom >= 0) {
                const scalar& w = node_solid_weight_y[bucket_idx_bottom][node_idx_bottom];
                rhs[dof_idx] += ((node_fluid_vel_y[bucket_idx_bottom][node_idx_bottom] *
                                        node_psi_fs_y[bucket_idx_bottom][node_idx_bottom] +
                                        node_elasto_vel_y[bucket_idx_bottom][node_idx_bottom] *
                                        node_psi_sf_y[bucket_idx_bottom][node_idx_bottom]) * w +
                                       node_solid_vel_y[bucket_idx_bottom][node_idx_bottom] * (1.0 - w)) / dx;
                
                if(w > 0.0) {
                    scalar term = w * coeff;
                    
                    term *= (1.0 - node_psi_y[bucket_idx_bottom][node_idx_bottom]) * node_inv_C_y[bucket_idx_bottom][node_idx_bottom] *
                    (node_fluid_vol_y[bucket_idx_bottom][node_idx_bottom] + node_elasto_vol_y[bucket_idx_bottom][node_idx_bottom] * node_mshdvm_hdvm_y[bucket_idx_bottom][node_idx_bottom])
                    + node_psi_y[bucket_idx_bottom][node_idx_bottom] * node_inv_Cs_y[bucket_idx_bottom][node_idx_bottom] *
                    (node_elasto_vol_y[bucket_idx_bottom][node_idx_bottom] + node_fluid_vol_y[bucket_idx_bottom][node_idx_bottom] * node_mfhdvm_hdvm_y[bucket_idx_bottom][node_idx_bottom]);
                    
                    const int bucket_pressure_bottom = node_pressure_index_y[bucket_idx_bottom][node_idx_bottom * 4 + 0];
                    const int node_pressure_bottom = node_pressure_index_y[bucket_idx_bottom][node_idx_bottom * 4 + 1];
                    
                    if(bucket_pressure_bottom >= 0 && node_pressure_bottom >= 0) {
                        const scalar bottom_phi = node_liquid_phi[bucket_pressure_bottom][node_pressure_bottom];
                        
                        if(bottom_phi < 0.0) {
                            const int dof_bottom = node_global_indices[bucket_pressure_bottom][node_pressure_bottom];
                            assert(dof_bottom >= 0);
                            
                            if(dof_bottom >= 0) {
                                matrix.add_to_element(dof_idx, dof_idx, term);
                                matrix.add_to_element(dof_idx, dof_bottom, -term);
                            }
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    } else { // coarse grid has zero pressure in fine level
                        const Vector3s& center_pos = bucket_pos_p.segment<3>(node_idx * 3);
                        const scalar bottom_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., dx, 0. ));
                        
                        if(bottom_phi < 0.0) {
                            matrix.add_to_element(dof_idx, dof_idx, term);
                            rhs[dof_idx] += term * scene.interpolateBucketPressure(center_pos - Vector3s(0, dx, 0));
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    }
                }
            } else {
                const Vector3s& np = bucket_pos_p.segment<3>(node_idx * 3);
                Vector3s nnp = np - Vector3s(0.0, 0.5 * dx, 0.0);
                const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                
                if(w > 0.0) {
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    const scalar bottom_phi = scene.interpolateBucketLiquidPhi(np - Vector3s( 0., dx, 0. ));
                    rhs[dof_idx] += (scene.interpolateBucketFluidVelY(nnp) * w + scene.interpolateBucketSolidVelY(nnp) * (1.0 - w)) / dx;
                    
                    if(bottom_phi < 0.0) {
                        matrix.add_to_element(dof_idx, dof_idx, term);
                        rhs[dof_idx] += term * scene.interpolateBucketPressure(np - Vector3s(0, dx, 0));
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
                        matrix.add_to_element(dof_idx, dof_idx, term / theta);
                    }
                }
            }
            
            const int bucket_idx_top = bucket_pn[node_idx * 12 + 6];
            const int node_idx_top = bucket_pn[node_idx * 12 + 7];
            
            if(bucket_idx_top >= 0 && node_idx_top >= 0) {
                const scalar& w = node_solid_weight_y[bucket_idx_top][node_idx_top];
                rhs[dof_idx] -= ((node_fluid_vel_y[bucket_idx_top][node_idx_top] *
                                        node_psi_fs_y[bucket_idx_top][node_idx_top] +
                                        node_elasto_vel_y[bucket_idx_top][node_idx_top] *
                                        node_psi_sf_y[bucket_idx_top][node_idx_top]) * w +
                                       node_solid_vel_y[bucket_idx_top][node_idx_top] * (1.0 - w)) / dx;
                
                if(w > 0.0) {
                    scalar term = w * coeff;
                    
                    term *= (1.0 - node_psi_y[bucket_idx_top][node_idx_top]) * node_inv_C_y[bucket_idx_top][node_idx_top] *
                    (node_fluid_vol_y[bucket_idx_top][node_idx_top] + node_elasto_vol_y[bucket_idx_top][node_idx_top] * node_mshdvm_hdvm_y[bucket_idx_top][node_idx_top])
                    + node_psi_y[bucket_idx_top][node_idx_top] * node_inv_Cs_y[bucket_idx_top][node_idx_top] *
                    (node_elasto_vol_y[bucket_idx_top][node_idx_top] + node_fluid_vol_y[bucket_idx_top][node_idx_top] * node_mfhdvm_hdvm_y[bucket_idx_top][node_idx_top]);
                    
                    const int bucket_pressure_top = node_pressure_index_y[bucket_idx_top][node_idx_top * 4 + 2];
                    const int node_pressure_top = node_pressure_index_y[bucket_idx_top][node_idx_top * 4 + 3];
                    
                    if(bucket_pressure_top >= 0 && node_pressure_top >= 0) {
                        const scalar top_phi = node_liquid_phi[bucket_pressure_top][node_pressure_top];
                        
                        if(top_phi < 0.0) {
                            const int dof_top = node_global_indices[bucket_pressure_top][node_pressure_top];
                            assert(dof_top >= 0);
                            
                            if(dof_top >= 0) {
                                matrix.add_to_element(dof_idx, dof_idx, term);
                                matrix.add_to_element(dof_idx, dof_top, -term);
                            }
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    } else { // coarse grid has zero pressure in fine level
                        const Vector3s& center_pos = bucket_pos_p.segment<3>(node_idx * 3);
                        const scalar top_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., dx, 0. ));
                        
                        if(top_phi < 0.0) {
                            rhs[dof_idx] += term * scene.interpolateBucketPressure(center_pos + Vector3s(0, dx, 0));
                            matrix.add_to_element(dof_idx, dof_idx, term);
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    }
                }
            } else {
                const Vector3s& np = bucket_pos_p.segment<3>(node_idx * 3);
                Vector3s nnp = np + Vector3s(0.0, 0.5 * dx, 0.0);
                const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                
                if(w > 0.0) {
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    const scalar top_phi = scene.interpolateBucketLiquidPhi(np + Vector3s( 0., dx, 0. ));
                    rhs[dof_idx] -= (scene.interpolateBucketFluidVelY(nnp) * w + scene.interpolateBucketSolidVelY(nnp) * (1.0 - w)) / dx;
                    
                    if(top_phi < 0.0) {
                        rhs[dof_idx] += term * scene.interpolateBucketPressure(np + Vector3s(0, dx, 0));
                        matrix.add_to_element(dof_idx, dof_idx, term);
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
                        matrix.add_to_element(dof_idx, dof_idx, term / theta);
                    }
                }
            }
            
            const int bucket_idx_near = bucket_pn[node_idx * 12 + 8];
            const int node_idx_near = bucket_pn[node_idx * 12 + 9];
            
            if(bucket_idx_near >= 0 && node_idx_near >= 0) {
                const scalar& w = node_solid_weight_z[bucket_idx_near][node_idx_near];
                rhs[dof_idx] += ((node_fluid_vel_z[bucket_idx_near][node_idx_near] *
                                        node_psi_fs_z[bucket_idx_near][node_idx_near] +
                                        node_elasto_vel_z[bucket_idx_near][node_idx_near] *
                                        node_psi_sf_z[bucket_idx_near][node_idx_near]) * w +
                                       node_solid_vel_z[bucket_idx_near][node_idx_near] * (1.0 - w)) / dx;
                
                if(w > 0.0) {
                    scalar term = w * coeff;
                    
                    term *= (1.0 - node_psi_z[bucket_idx_near][node_idx_near]) * node_inv_C_z[bucket_idx_near][node_idx_near] *
                    (node_fluid_vol_z[bucket_idx_near][node_idx_near] + node_elasto_vol_z[bucket_idx_near][node_idx_near] * node_mshdvm_hdvm_z[bucket_idx_near][node_idx_near])
                    + node_psi_z[bucket_idx_near][node_idx_near] * node_inv_Cs_z[bucket_idx_near][node_idx_near] *
                    (node_elasto_vol_z[bucket_idx_near][node_idx_near] + node_fluid_vol_z[bucket_idx_near][node_idx_near] * node_mfhdvm_hdvm_z[bucket_idx_near][node_idx_near]);
                    
                    const int bucket_pressure_near = node_pressure_index_z[bucket_idx_near][node_idx_near * 4 + 0];
                    const int node_pressure_near = node_pressure_index_z[bucket_idx_near][node_idx_near * 4 + 1];
                    
                    if(bucket_pressure_near >= 0 && node_pressure_near >= 0) {
                        const scalar near_phi = node_liquid_phi[bucket_pressure_near][node_pressure_near];
                        
                        if(near_phi < 0.0) {
                            const int dof_near = node_global_indices[bucket_pressure_near][node_pressure_near];
                            assert(dof_near >= 0);
                            
                            if(dof_near >= 0) {
                                matrix.add_to_element(dof_idx, dof_idx, term);
                                matrix.add_to_element(dof_idx, dof_near, -term);
                            }
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    } else { // coarse grid has zero pressure in fine level
                        const Vector3s& center_pos = bucket_pos_p.segment<3>(node_idx * 3);
                        const scalar near_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., 0., dx ));
                        
                        if(near_phi < 0.0) {
                            matrix.add_to_element(dof_idx, dof_idx, term);
                            rhs[dof_idx] += term * scene.interpolateBucketPressure(center_pos - Vector3s(0, 0, dx));
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    }
                }
            } else {
                const Vector3s& np = bucket_pos_p.segment<3>(node_idx * 3);
                Vector3s nnp = np - Vector3s(0.0, 0.0, 0.5 * dx);
                const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                
                if(w > 0.0) {
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    const scalar near_phi = scene.interpolateBucketLiquidPhi(np - Vector3s( 0., 0., dx ));
                    rhs[dof_idx] += (scene.interpolateBucketFluidVelZ(nnp) * w + scene.interpolateBucketSolidVelZ(nnp) * (1.0 - w)) / dx;
                    
                    if(near_phi < 0.0) {
                        matrix.add_to_element(dof_idx, dof_idx, term);
                        rhs[dof_idx] += term * scene.interpolateBucketPressure(np - Vector3s(0, 0, dx));
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
                        matrix.add_to_element(dof_idx, dof_idx, term / theta);
                    }
                }
            }
            
            const int bucket_idx_far = bucket_pn[node_idx * 12 + 10];
            const int node_idx_far = bucket_pn[node_idx * 12 + 11];
            
            if(bucket_idx_far >= 0 && node_idx_far >= 0) {
                const scalar& w = node_solid_weight_z[bucket_idx_far][node_idx_far];
                rhs[dof_idx] -= ((node_fluid_vel_z[bucket_idx_far][node_idx_far] *
                                        node_psi_fs_z[bucket_idx_far][node_idx_far] +
                                        node_elasto_vel_z[bucket_idx_far][node_idx_far] *
                                        node_psi_sf_z[bucket_idx_far][node_idx_far]) * w +
                                       node_solid_vel_z[bucket_idx_far][node_idx_far] * (1.0 - w)) / dx;
                
                if(w > 0.0) {
                    scalar term = w * coeff;
                    
                    term *= (1.0 - node_psi_z[bucket_idx_far][node_idx_far]) * node_inv_C_z[bucket_idx_far][node_idx_far] *
                    (node_fluid_vol_z[bucket_idx_far][node_idx_far] + node_elasto_vol_z[bucket_idx_far][node_idx_far] * node_mshdvm_hdvm_z[bucket_idx_far][node_idx_far])
                    + node_psi_z[bucket_idx_far][node_idx_far] * node_inv_Cs_z[bucket_idx_far][node_idx_far] *
                    (node_elasto_vol_z[bucket_idx_far][node_idx_far] + node_fluid_vol_z[bucket_idx_far][node_idx_far] * node_mfhdvm_hdvm_z[bucket_idx_far][node_idx_far]);
                    
                    const int bucket_pressure_far = node_pressure_index_z[bucket_idx_far][node_idx_far * 4 + 2];
                    const int node_pressure_far = node_pressure_index_z[bucket_idx_far][node_idx_far * 4 + 3];
                    
                    if(bucket_pressure_far >= 0 && node_pressure_far >= 0) {
                        const scalar far_phi = node_liquid_phi[bucket_pressure_far][node_pressure_far];
                        
                        if(far_phi < 0.0) {
                            const int dof_far = node_global_indices[bucket_pressure_far][node_pressure_far];
                            assert(dof_far >= 0);
                            
                            if(dof_far >= 0) {
                                matrix.add_to_element(dof_idx, dof_idx, term);
                                matrix.add_to_element(dof_idx, dof_far, -term);
                            }
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    } else { // coarse grid has zero pressure in fine level
                        const Vector3s& center_pos = bucket_pos_p.segment<3>(node_idx * 3);
                        const scalar far_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., 0., dx ));
                        
                        if(far_phi < 0.0) {
                            matrix.add_to_element(dof_idx, dof_idx, term);
                            rhs[dof_idx] += term * scene.interpolateBucketPressure(center_pos + Vector3s(0, 0, dx));
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
                            matrix.add_to_element(dof_idx, dof_idx, term / theta);
                        }
                    }
                }
            } else {
                const Vector3s& np = bucket_pos_p.segment<3>(node_idx * 3);
                Vector3s nnp = np + Vector3s(0.0, 0.0, 0.5 * dx);
                const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                
                if(w > 0.0) {
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    rhs[dof_idx] -= (scene.interpolateBucketFluidVelZ(nnp) * w + scene.interpolateBucketSolidVelZ(nnp) * (1.0 - w)) / dx;
                    
                    const scalar far_phi = scene.interpolateBucketLiquidPhi(np + Vector3s( 0., 0., dx ));
                    
                    if(far_phi < 0.0) {
                        matrix.add_to_element(dof_idx, dof_idx, term);
                        rhs[dof_idx] += term * scene.interpolateBucketPressure(np + Vector3s(0, 0, dx));
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
                        matrix.add_to_element(dof_idx, dof_idx, term / theta);
                    }
                }
            }
        });
        
        bool success = false;
        scalar tolerance = 0.0;
        int iterations = 0;
        
        success = AMGPCGSolveSparse(matrix, rhs, result, dof_ijk, criterion, maxiters, tolerance, iterations, ni, nj, nk);
        
        std::cout << "[amg pcg total iter: " << iterations << ", res: " << tolerance << "]" << std::endl;
        
        if(!success) {
            std::cout << "WARNING: AMG PCG solve failed!" << std::endl;
            
            std::cout << "rhs=[";
            for(scalar s : rhs) {
                std::cout << s << "; ";
            }
            std::cout << "];" << std::endl;
        }
		
#ifdef CHECK_AMGPCG_RESULT
		std::vector<scalar> tmp = rhs;
		
		robertbridson::multiply_and_subtract(matrix, result, tmp);
		
		scalar residual = Eigen::Map<VectorXs>((scalar*)&tmp[0], tmp.size()).norm();
		
		scalar len_rhs = Eigen::Map<VectorXs>((scalar*)&rhs[0], rhs.size()).norm();
		
		std::cout << "[amg pcg check result: " << residual << ", " << len_rhs << ", " << (residual / len_rhs) << "]" << std::endl;
#endif
		
        threadutils::for_each(0, total_num_nodes, [&] (int dof_idx) {
            const Vector2i& dof_loc = effective_node_indices[dof_idx];
            pressure[dof_loc[0]][dof_loc[1]] = result[dof_idx];
        });
    }
    
	void multiplyPressureMatrix( const TwoDScene& scene, const std::vector< VectorXs >& node_vec, std::vector< VectorXs >& out_node_vec, const std::vector< VectorXs >& node_inv_mdv_x, const std::vector< VectorXs >& node_inv_mdv_y, const std::vector< VectorXs >& node_inv_mdv_z, const std::vector< VectorXs >& node_inv_mdvs_x, const std::vector< VectorXs >& node_inv_mdvs_y, const std::vector< VectorXs >& node_inv_mdvs_z, const scalar& dt )
	{
		const Sorter& buckets = scene.getParticleBuckets();
		
		const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
		const std::vector< VectorXi >& pressure_neighbors = scene.getPressureNeighbors();
		const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
		const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
		const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
		const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
		const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
		const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
		const std::vector< VectorXs >& node_fluid_vol_x = scene.getNodeFluidVolX();
		const std::vector< VectorXs >& node_fluid_vol_y = scene.getNodeFluidVolY();
		const std::vector< VectorXs >& node_fluid_vol_z = scene.getNodeFluidVolZ();
		const std::vector< VectorXs >& node_elasto_vol_x = scene.getNodeVolX();
		const std::vector< VectorXs >& node_elasto_vol_y = scene.getNodeVolY();
		const std::vector< VectorXs >& node_elasto_vol_z = scene.getNodeVolZ();
		const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
		const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
		const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
		
		const scalar dx = scene.getCellSize();
		const scalar dV = dx * dx * dx;
		const scalar coeff = dt / (dx * dx);
		
		buckets.for_each_bucket([&] (int bucket_idx) {
			const VectorXs& bucket_liquid_phi = node_liquid_phi[bucket_idx];
			const VectorXi& bucket_pn = pressure_neighbors[bucket_idx];;
            const VectorXs& bucket_pos_p = scene.getNodePosP(bucket_idx);
            
			const int num_pressure = out_node_vec[bucket_idx].size();
			
			assert(out_node_vec[bucket_idx].size() == bucket_liquid_phi.size());
			
			VectorXs& bucket_out = out_node_vec[bucket_idx];
			
			for(int i = 0; i < num_pressure; ++i)
			{
				const scalar center_phi = bucket_liquid_phi(i);
				
				bucket_out(i) = 0.0;
				
				if(center_phi > 0.0) continue;
                const Vector3s center_pos = bucket_pos_p.segment<3>(i * 3);
				
				const int bucket_idx_left = bucket_pn[i * 12 + 0];
				const int node_idx_left = bucket_pn[i * 12 + 1];
				
				if(bucket_idx_left >= 0 && node_idx_left >= 0) {
					const scalar& w = node_solid_weight_x[bucket_idx_left][node_idx_left];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_x[bucket_idx_left][node_idx_left] * (1.0 - node_psi_x[bucket_idx_left][node_idx_left]) * node_inv_mdv_x[bucket_idx_left][node_idx_left] + node_elasto_vol_x[bucket_idx_left][node_idx_left] * node_psi_x[bucket_idx_left][node_idx_left] * node_inv_mdvs_x[bucket_idx_left][node_idx_left];
					
					const int bucket_pressure_left = node_pressure_index_x[bucket_idx_left][node_idx_left * 4 + 0];
					const int node_pressure_left = node_pressure_index_x[bucket_idx_left][node_idx_left * 4 + 1];
					
					if(bucket_pressure_left >= 0 && node_pressure_left >= 0) {
						const scalar left_phi = node_liquid_phi[bucket_pressure_left][node_pressure_left];
						
						if(left_phi < 0.0) {
							bucket_out(i) += term * (node_vec[bucket_idx][i] - node_vec[bucket_pressure_left][node_pressure_left]);
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
							bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
						}
                    } else { // coarse grid has zero pressure in fine level
                        const scalar left_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( dx, 0., 0. ));
                        
                        if(left_phi < 0.0) {
                            bucket_out(i) += term * node_vec[bucket_idx][i];
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
                            bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.5 * dx, 0.0, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar left_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( dx, 0., 0. ));
                    
                    if(left_phi < 0.0) {
                        bucket_out(i) += term * node_vec[bucket_idx][i];
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
                        bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                    }
                }
				
				const int bucket_idx_right = bucket_pn[i * 12 + 2];
				const int node_idx_right = bucket_pn[i * 12 + 3];
				
				if(bucket_idx_right >= 0 && node_idx_right >= 0) {
					const scalar& w = node_solid_weight_x[bucket_idx_right][node_idx_right];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_x[bucket_idx_right][node_idx_right] * (1.0 - node_psi_x[bucket_idx_right][node_idx_right]) * node_inv_mdv_x[bucket_idx_right][node_idx_right] + node_elasto_vol_x[bucket_idx_right][node_idx_right] * node_psi_x[bucket_idx_right][node_idx_right] * node_inv_mdvs_x[bucket_idx_right][node_idx_right];
					
					const int bucket_pressure_right = node_pressure_index_x[bucket_idx_right][node_idx_right * 4 + 2];
					const int node_pressure_right = node_pressure_index_x[bucket_idx_right][node_idx_right * 4 + 3];
					
					if(bucket_pressure_right >= 0 && node_pressure_right >= 0) {
						const scalar right_phi = node_liquid_phi[bucket_pressure_right][node_pressure_right];
						
						if(right_phi < 0.0) {
							bucket_out(i) += term * (node_vec[bucket_idx][i] - node_vec[bucket_pressure_right][node_pressure_right]);
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
							bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
						}
                    } else { // coarse grid has zero pressure in fine level
                        const scalar right_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( dx, 0., 0. ));
                        
                        if(right_phi < 0.0) {
                            bucket_out(i) += term * node_vec[bucket_idx][i];
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
                            bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.5 * dx, 0.0, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar right_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( dx, 0., 0. ));
                    
                    if(right_phi < 0.0) {
                        bucket_out(i) += term * node_vec[bucket_idx][i];
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
                        bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                    }
                }
				
				const int bucket_idx_bottom = bucket_pn[i * 12 + 4];
				const int node_idx_bottom = bucket_pn[i * 12 + 5];
				
				if(bucket_idx_bottom >= 0 && node_idx_bottom >= 0) {
					const scalar& w = node_solid_weight_y[bucket_idx_bottom][node_idx_bottom];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_y[bucket_idx_bottom][node_idx_bottom] * (1.0 - node_psi_y[bucket_idx_bottom][node_idx_bottom]) * node_inv_mdv_y[bucket_idx_bottom][node_idx_bottom] + node_elasto_vol_y[bucket_idx_bottom][node_idx_bottom] * node_psi_y[bucket_idx_bottom][node_idx_bottom] * node_inv_mdvs_y[bucket_idx_bottom][node_idx_bottom];
					
					const int bucket_pressure_bottom = node_pressure_index_y[bucket_idx_bottom][node_idx_bottom * 4 + 0];
					const int node_pressure_bottom = node_pressure_index_y[bucket_idx_bottom][node_idx_bottom * 4 + 1];
					
					if(bucket_pressure_bottom >= 0 && node_pressure_bottom >= 0) {
						const scalar bottom_phi = node_liquid_phi[bucket_pressure_bottom][node_pressure_bottom];
						
						if(bottom_phi < 0.0) {
							bucket_out(i) += term * (node_vec[bucket_idx][i] - node_vec[bucket_pressure_bottom][node_pressure_bottom]);
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
							bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
						}
                    } else { // coarse grid has zero pressure in fine level
                        const scalar bottom_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., dx, 0. ));
                        
                        if(bottom_phi < 0.0) {
                            bucket_out(i) += term * node_vec[bucket_idx][i];
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
                            bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.0, 0.5 * dx, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar bottom_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., dx, 0. ));
                    
                    if(bottom_phi < 0.0) {
                        bucket_out(i) += term * node_vec[bucket_idx][i];
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
                        bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                    }
                }
				
				const int bucket_idx_top = bucket_pn[i * 12 + 6];
				const int node_idx_top = bucket_pn[i * 12 + 7];
				
				if(bucket_idx_top >= 0 && node_idx_top >= 0) {
					const scalar& w = node_solid_weight_y[bucket_idx_top][node_idx_top];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_y[bucket_idx_top][node_idx_top] * (1.0 - node_psi_y[bucket_idx_top][node_idx_top]) * node_inv_mdv_y[bucket_idx_top][node_idx_top] + node_elasto_vol_y[bucket_idx_top][node_idx_top] * node_psi_y[bucket_idx_top][node_idx_top] * node_inv_mdvs_y[bucket_idx_top][node_idx_top];
					
					const int bucket_pressure_top = node_pressure_index_y[bucket_idx_top][node_idx_top * 4 + 2];
					const int node_pressure_top = node_pressure_index_y[bucket_idx_top][node_idx_top * 4 + 3];
					
					if(bucket_pressure_top >= 0 && node_pressure_top >= 0) {
						const scalar top_phi = node_liquid_phi[bucket_pressure_top][node_pressure_top];
						
						if(top_phi < 0.0) {
							bucket_out(i) += term * (node_vec[bucket_idx][i] - node_vec[bucket_pressure_top][node_pressure_top]);
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
							bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
						}
                    } else { // coarse grid has zero pressure in fine level
                        const scalar top_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., dx, 0. ));
                        
                        if(top_phi < 0.0) {
                            bucket_out(i) += term * node_vec[bucket_idx][i];
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
                            bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.0, 0.5 * dx, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar top_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., dx, 0. ));
                    
                    if(top_phi < 0.0) {
                        bucket_out(i) += term * node_vec[bucket_idx][i];
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
                        bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                    }
                }
				
				const int bucket_idx_near = bucket_pn[i * 12 + 8];
				const int node_idx_near = bucket_pn[i * 12 + 9];
				
				if(bucket_idx_near >= 0 && node_idx_near >= 0) {
					const scalar& w = node_solid_weight_z[bucket_idx_near][node_idx_near];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_z[bucket_idx_near][node_idx_near] * (1.0 - node_psi_z[bucket_idx_near][node_idx_near]) * node_inv_mdv_z[bucket_idx_near][node_idx_near] + node_elasto_vol_z[bucket_idx_near][node_idx_near] * node_psi_z[bucket_idx_near][node_idx_near] * node_inv_mdvs_z[bucket_idx_near][node_idx_near];
					
					const int bucket_pressure_near = node_pressure_index_z[bucket_idx_near][node_idx_near * 4 + 0];
					const int node_pressure_near = node_pressure_index_z[bucket_idx_near][node_idx_near * 4 + 1];
					
					if(bucket_pressure_near >= 0 && node_pressure_near >= 0) {
						const scalar near_phi = node_liquid_phi[bucket_pressure_near][node_pressure_near];
						
						if(near_phi < 0.0) {
							bucket_out(i) += term * (node_vec[bucket_idx][i] - node_vec[bucket_pressure_near][node_pressure_near]);
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
							bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
						}
                    } else { // coarse grid has zero pressure in fine level
                        const scalar near_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., 0., dx ));
                        
                        if(near_phi < 0.0) {
                            bucket_out(i) += term * node_vec[bucket_idx][i];
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
                            bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.0, 0.0, 0.5 * dx);
                    const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar near_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., 0., dx ));
                    
                    if(near_phi < 0.0) {
                        bucket_out(i) += term * node_vec[bucket_idx][i];
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
                        bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                    }
                }
				
				const int bucket_idx_far = bucket_pn[i * 12 + 10];
				const int node_idx_far = bucket_pn[i * 12 + 11];
				
				if(bucket_idx_far >= 0 && node_idx_far >= 0) {
					const scalar& w = node_solid_weight_z[bucket_idx_far][node_idx_far];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_z[bucket_idx_far][node_idx_far] * (1.0 - node_psi_z[bucket_idx_far][node_idx_far]) * node_inv_mdv_z[bucket_idx_far][node_idx_far] + node_elasto_vol_z[bucket_idx_far][node_idx_far] * node_psi_z[bucket_idx_far][node_idx_far] * node_inv_mdvs_z[bucket_idx_far][node_idx_far];
					
					const int bucket_pressure_far = node_pressure_index_z[bucket_idx_far][node_idx_far * 4 + 2];
					const int node_pressure_far = node_pressure_index_z[bucket_idx_far][node_idx_far * 4 + 3];
					
					if(bucket_pressure_far >= 0 && node_pressure_far >= 0) {
						const scalar far_phi = node_liquid_phi[bucket_pressure_far][node_pressure_far];
						
						if(far_phi < 0.0) {
							bucket_out(i) += term * (node_vec[bucket_idx][i] - node_vec[bucket_pressure_far][node_pressure_far]);
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
							bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
						}
                    } else { // coarse grid has zero pressure in fine level
                        const scalar far_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., 0., dx ));
                        
                        if(far_phi < 0.0) {
                            bucket_out(i) += term * node_vec[bucket_idx][i];
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
                            bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.0, 0.0, 0.5 * dx);
                    const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar far_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., 0., dx ));
                    
                    if(far_phi < 0.0) {
                        bucket_out(i) += term * node_vec[bucket_idx][i];
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
                        bucket_out(i) += term * node_vec[bucket_idx][i] / theta;
                    }
                }
			}
		});
	}
	
	void allocateNodes( const TwoDScene& scene, std::vector< VectorXs >& out_node_vec )
	{
		const Sorter& buckets = scene.getParticleBuckets();
		const std::vector< VectorXi >& pressure_neighbors = scene.getPressureNeighbors();
		out_node_vec.resize(buckets.size());
		buckets.for_each_bucket([&] (int bucket_idx) {
			const VectorXi& bucket_pn = pressure_neighbors[bucket_idx];;
			const int num_pressure = bucket_pn.size() / 12;
			VectorXs& bucket_out = out_node_vec[bucket_idx];
			bucket_out.resize(num_pressure);
		});
	}
	
	void localSolveJacobi( const std::vector< VectorXs >& node_vec, std::vector< VectorXs >& out_node_vec, const std::vector< VectorXs >& jacobi_precond)
	{
		const int num_buckets = node_vec.size();
		threadutils::for_each(0, num_buckets, [&] (int bucket_idx) {
			out_node_vec[bucket_idx] = node_vec[bucket_idx].array() * jacobi_precond[bucket_idx].array();
		});
	}
	
    void applyPressureGradsElastoRHSBiCGSTAB( const TwoDScene& scene,
                                     const std::vector< VectorXs >& pressure_vec,
                                     std::vector< VectorXs >& rhs_vec_x,
                                     std::vector< VectorXs >& rhs_vec_y,
                                     std::vector< VectorXs >& rhs_vec_z,
                                     const scalar& dt )
    {
        const Sorter& buckets = scene.getParticleBuckets();
        
        const scalar dx = scene.getCellSize();
        const scalar dV = dx * dx * dx;
        
        const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
        
        const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
        const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
        const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
        const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
        
        const std::vector< VectorXs >& node_vol_x = scene.getNodeVolX();
        const std::vector< VectorXs >& node_vol_y = scene.getNodeVolY();
        const std::vector< VectorXs >& node_vol_z = scene.getNodeVolZ();
        
        const std::vector< VectorXs >& node_vol_fluid_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_vol_fluid_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_vol_fluid_z = scene.getNodeFluidVolZ();
        
        const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
        const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
        const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
        
        buckets.for_each_bucket([&] (int bucket_idx) {
            VectorXs& bucket_rhs_x = rhs_vec_x[bucket_idx];
            VectorXs& bucket_rhs_y = rhs_vec_y[bucket_idx];
            VectorXs& bucket_rhs_z = rhs_vec_z[bucket_idx];
            
            const VectorXs& bucket_weight_x = node_solid_weight_x[bucket_idx];
            const VectorXs& bucket_weight_y = node_solid_weight_y[bucket_idx];
            const VectorXs& bucket_weight_z = node_solid_weight_z[bucket_idx];
            
            const VectorXs& bucket_pos_x = scene.getNodePosX(bucket_idx);
            const VectorXs& bucket_pos_y = scene.getNodePosY(bucket_idx);
            const VectorXs& bucket_pos_z = scene.getNodePosZ(bucket_idx);
            
            const int num_rhs_x = bucket_rhs_x.size();
            const int num_rhs_y = bucket_rhs_y.size();
            const int num_rhs_z = bucket_rhs_z.size();
            
            for(int i = 0; i < num_rhs_x; ++i) {
                const Vector3s& pos_x = bucket_pos_x.segment<3>(i * 3);
                
                if(bucket_weight_x(i) > 0.) {
                    const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);
                    
                    scalar left_phi, right_phi;
                    if(indices[0] < 0 || indices[1] < 0) left_phi = scene.interpolateBucketLiquidPhi(pos_x - Vector3s(0.5 * dx, 0.0, 0.0));
                    else left_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) right_phi = scene.interpolateBucketLiquidPhi(pos_x + Vector3s(0.5 * dx, 0.0, 0.0));
                    else right_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(left_phi < 0.0 || right_phi < 0.0) {
                        scalar theta = 1.;
                        if(right_phi >= 0.0 || left_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_x - Vector3s(dx * 0.5, 0, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_x + Vector3s(dx * 0.5, 0, 0));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt * node_vol_x[bucket_idx][i];
                        
                        bucket_rhs_x(i) += term;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_y; ++i) {
                const Vector3s& pos_y = bucket_pos_y.segment<3>(i * 3);
                
                if(bucket_weight_y(i) > 0.) {
                    const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);
                    
                    scalar bottom_phi, top_phi;
                    if(indices[0] < 0 || indices[1] < 0) bottom_phi = scene.interpolateBucketLiquidPhi(pos_y - Vector3s(0.0, 0.5 * dx, 0.0));
                    else bottom_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) top_phi = scene.interpolateBucketLiquidPhi(pos_y + Vector3s(0.0, 0.5 * dx, 0.0));
                    else top_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(bottom_phi < 0.0 || top_phi < 0.0) {
                        scalar theta = 1.;
                        if(top_phi >= 0.0 || bottom_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(bottom_phi, top_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_y - Vector3s(0, dx * 0.5, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_y + Vector3s(0, dx * 0.5, 0));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt * node_vol_y[bucket_idx][i];
                        
                        bucket_rhs_y(i) += term;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_z; ++i) {
                const Vector3s& pos_z = bucket_pos_z.segment<3>(i * 3);
                
                if(bucket_weight_z(i) > 0.) {
                    const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);
                    
                    scalar near_phi, far_phi;
                    if(indices[0] < 0 || indices[1] < 0) near_phi = scene.interpolateBucketLiquidPhi(pos_z - Vector3s(0.0, 0.0, 0.5 * dx));
                    else near_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) far_phi = scene.interpolateBucketLiquidPhi(pos_z + Vector3s(0.0, 0.0, 0.5 * dx));
                    else far_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(near_phi < 0.0 || far_phi < 0.0) {
                        scalar theta = 1.;
                        if(far_phi >= 0.0 || near_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(near_phi, far_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_z - Vector3s(0, 0, dx * 0.5));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_z + Vector3s(0, 0, dx * 0.5));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt * node_vol_z[bucket_idx][i];
                        
                        bucket_rhs_z(i) += term;
                    }
                }
            }
        });
    }
    
    void applyPressureGradsFluidRHS( TwoDScene& scene,
                                     const std::vector< VectorXs >& pressure_vec,
                                     std::vector< VectorXs >& rhs_vec_x,
                                     std::vector< VectorXs >& rhs_vec_y,
                                     std::vector< VectorXs >& rhs_vec_z,
                                     const std::vector< VectorXs >& node_mshdvm_hdvm_x, // (M_f+hDVm)^{-1}hDVm
                                     const std::vector< VectorXs >& node_mshdvm_hdvm_y,
                                     const std::vector< VectorXs >& node_mshdvm_hdvm_z,
                                     const scalar& dt )
    {
        const Sorter& buckets = scene.getParticleBuckets();
        
        const scalar dx = scene.getCellSize();
        const scalar dV = dx * dx * dx;
        
        std::vector< VectorXuc >& node_liquid_valid_x = scene.getNodeLiquidValidX();
        std::vector< VectorXuc >& node_liquid_valid_y = scene.getNodeLiquidValidY();
        std::vector< VectorXuc >& node_liquid_valid_z = scene.getNodeLiquidValidZ();
        
        const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
        
        const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
        const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
        const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
        const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
        
        const std::vector< VectorXs >& node_vol_x = scene.getNodeVolX();
        const std::vector< VectorXs >& node_vol_y = scene.getNodeVolY();
        const std::vector< VectorXs >& node_vol_z = scene.getNodeVolZ();
        
        const std::vector< VectorXs >& node_vol_fluid_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_vol_fluid_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_vol_fluid_z = scene.getNodeFluidVolZ();
        
        const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
        const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
        const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
        
        buckets.for_each_bucket([&] (int bucket_idx) {
            VectorXs& bucket_rhs_x = rhs_vec_x[bucket_idx];
            VectorXs& bucket_rhs_y = rhs_vec_y[bucket_idx];
            VectorXs& bucket_rhs_z = rhs_vec_z[bucket_idx];
            
            const VectorXs& bucket_weight_x = node_solid_weight_x[bucket_idx];
            const VectorXs& bucket_weight_y = node_solid_weight_y[bucket_idx];
            const VectorXs& bucket_weight_z = node_solid_weight_z[bucket_idx];
            
            VectorXuc& bucket_valid_x = node_liquid_valid_x[bucket_idx];
            VectorXuc& bucket_valid_y = node_liquid_valid_y[bucket_idx];
            VectorXuc& bucket_valid_z = node_liquid_valid_z[bucket_idx];
            
            const VectorXs& bucket_pos_x = scene.getNodePosX(bucket_idx);
            const VectorXs& bucket_pos_y = scene.getNodePosY(bucket_idx);
            const VectorXs& bucket_pos_z = scene.getNodePosZ(bucket_idx);
            
            const int num_rhs_x = bucket_rhs_x.size();
            const int num_rhs_y = bucket_rhs_y.size();
            const int num_rhs_z = bucket_rhs_z.size();
            
            for(int i = 0; i < num_rhs_x; ++i) {
                bucket_valid_x(i) = 0U;
                const Vector3s& pos_x = bucket_pos_x.segment<3>(i * 3);
                
                if(bucket_weight_x(i) > 0. ) {
                    const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);
                    
                    scalar left_phi, right_phi;
                    if(indices[0] < 0 || indices[1] < 0) left_phi = scene.interpolateBucketLiquidPhi(pos_x - Vector3s(0.5 * dx, 0.0, 0.0));
                    else left_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) right_phi = scene.interpolateBucketLiquidPhi(pos_x + Vector3s(0.5 * dx, 0.0, 0.0));
                    else right_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(left_phi < 0.0 || right_phi < 0.0) {
                        scalar theta = 1.;
                        if(right_phi >= 0.0 || left_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_x - Vector3s(dx * 0.5, 0, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_x + Vector3s(dx * 0.5, 0, 0));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt * (node_vol_fluid_x[bucket_idx][i] + node_vol_x[bucket_idx][i] * node_mshdvm_hdvm_x[bucket_idx][i]);
                        
                        bucket_rhs_x(i) -= term;
                        bucket_valid_x(i) = 1U;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_y; ++i) {
                bucket_valid_y(i) = 0U;
                const Vector3s& pos_y = bucket_pos_y.segment<3>(i * 3);
                
                if(bucket_weight_y(i) > 0. ) {
                    const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);
                    
                    scalar bottom_phi, top_phi;
                    if(indices[0] < 0 || indices[1] < 0) bottom_phi = scene.interpolateBucketLiquidPhi(pos_y - Vector3s(0.0, 0.5 * dx, 0.0));
                    else bottom_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) top_phi = scene.interpolateBucketLiquidPhi(pos_y + Vector3s(0.0, 0.5 * dx, 0.0));
                    else top_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(bottom_phi < 0.0 || top_phi < 0.0) {
                        scalar theta = 1.;
                        if(top_phi >= 0.0 || bottom_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(bottom_phi, top_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_y - Vector3s(0, dx * 0.5, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_y + Vector3s(0, dx * 0.5, 0));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt * (node_vol_fluid_y[bucket_idx][i] + node_vol_y[bucket_idx][i] * node_mshdvm_hdvm_y[bucket_idx][i]);
                        
                        bucket_rhs_y(i) -= term;
                        bucket_valid_y(i) = 1U;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_z; ++i) {
                bucket_valid_z(i) = 0U;
                const Vector3s& pos_z = bucket_pos_z.segment<3>(i * 3);
                
                if(bucket_weight_z(i) > 0. ) {
                    const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);
                    
                    scalar near_phi, far_phi;
                    if(indices[0] < 0 || indices[1] < 0) near_phi = scene.interpolateBucketLiquidPhi(pos_z - Vector3s(0.0, 0.0, 0.5 * dx));
                    else near_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) far_phi = scene.interpolateBucketLiquidPhi(pos_z + Vector3s(0.0, 0.0, 0.5 * dx));
                    else far_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(near_phi < 0.0 || far_phi < 0.0) {
                        scalar theta = 1.;
                        if(far_phi >= 0.0 || near_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(near_phi, far_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_z - Vector3s(0, 0, dx * 0.5));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_z + Vector3s(0, 0, dx * 0.5));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt * (node_vol_fluid_z[bucket_idx][i] + node_vol_z[bucket_idx][i] * node_mshdvm_hdvm_z[bucket_idx][i]);
                        
                        bucket_rhs_z(i) -= term;
                        bucket_valid_z(i) = 1U;
                    }
                }
            }
        });
    }
    
	void applyPressureGradsElastoRHS( TwoDScene& scene,
                                     const std::vector< VectorXs >& pressure_vec,
                                     std::vector< VectorXs >& rhs_vec_x,
                                     std::vector< VectorXs >& rhs_vec_y,
                                     std::vector< VectorXs >& rhs_vec_z,
                                     const std::vector< VectorXs >& node_mfhdvm_hdvm_x, // (M_f+hDVm)^{-1}hDVm
                                     const std::vector< VectorXs >& node_mfhdvm_hdvm_y,
                                     const std::vector< VectorXs >& node_mfhdvm_hdvm_z,
                                     const scalar& dt )
	{
		const Sorter& buckets = scene.getParticleBuckets();
		
		const scalar dx = scene.getCellSize();
		const scalar dV = dx * dx * dx;
		
		const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
		const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
		const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
		
		const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
		const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
		const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
		const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
		
		const std::vector< VectorXs >& node_vol_x = scene.getNodeVolX();
		const std::vector< VectorXs >& node_vol_y = scene.getNodeVolY();
		const std::vector< VectorXs >& node_vol_z = scene.getNodeVolZ();
        
        const std::vector< VectorXs >& node_vol_fluid_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_vol_fluid_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_vol_fluid_z = scene.getNodeFluidVolZ();
		
		const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
		const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
		const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
		
		buckets.for_each_bucket([&] (int bucket_idx) {
			VectorXs& bucket_rhs_x = rhs_vec_x[bucket_idx];
			VectorXs& bucket_rhs_y = rhs_vec_y[bucket_idx];
			VectorXs& bucket_rhs_z = rhs_vec_z[bucket_idx];
			
			const VectorXs& bucket_weight_x = node_solid_weight_x[bucket_idx];
			const VectorXs& bucket_weight_y = node_solid_weight_y[bucket_idx];
			const VectorXs& bucket_weight_z = node_solid_weight_z[bucket_idx];
			
            const VectorXs& bucket_pos_x = scene.getNodePosX(bucket_idx);
            const VectorXs& bucket_pos_y = scene.getNodePosY(bucket_idx);
            const VectorXs& bucket_pos_z = scene.getNodePosZ(bucket_idx);
            
			const int num_rhs_x = bucket_rhs_x.size();
			const int num_rhs_y = bucket_rhs_y.size();
			const int num_rhs_z = bucket_rhs_z.size();
			
			for(int i = 0; i < num_rhs_x; ++i) {
                const Vector3s& pos_x = bucket_pos_x.segment<3>(i * 3);
                
				if(bucket_weight_x(i) > 0.) {
					const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);
                    
                    scalar left_phi, right_phi;
                    if(indices[0] < 0 || indices[1] < 0) left_phi = scene.interpolateBucketLiquidPhi(pos_x - Vector3s(0.5 * dx, 0.0, 0.0));
                    else left_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) right_phi = scene.interpolateBucketLiquidPhi(pos_x + Vector3s(0.5 * dx, 0.0, 0.0));
                    else right_phi = node_liquid_phi[indices[2]][indices[3]];
                    
					if(left_phi < 0.0 || right_phi < 0.0) {
						scalar theta = 1.;
						if(right_phi >= 0.0 || left_phi >= 0.0)
							theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
						
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_x - Vector3s(dx * 0.5, 0, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_x + Vector3s(dx * 0.5, 0, 0));
                        
                        pressure_grad /= (dx * theta);
						
						scalar term = pressure_grad * dt * (node_vol_x[bucket_idx][i] + node_vol_fluid_x[bucket_idx][i] * node_mfhdvm_hdvm_x[bucket_idx][i]);
						
						bucket_rhs_x(i) -= term;
					}
				}
			}
			
			for(int i = 0; i < num_rhs_y; ++i) {
                const Vector3s& pos_y = bucket_pos_y.segment<3>(i * 3);
                
				if(bucket_weight_y(i) > 0.) {
					const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);

                    scalar bottom_phi, top_phi;
                    if(indices[0] < 0 || indices[1] < 0) bottom_phi = scene.interpolateBucketLiquidPhi(pos_y - Vector3s(0.0, 0.5 * dx, 0.0));
                    else bottom_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) top_phi = scene.interpolateBucketLiquidPhi(pos_y + Vector3s(0.0, 0.5 * dx, 0.0));
                    else top_phi = node_liquid_phi[indices[2]][indices[3]];
                    
					if(bottom_phi < 0.0 || top_phi < 0.0) {
						scalar theta = 1.;
						if(top_phi >= 0.0 || bottom_phi >= 0.0)
							theta = std::max(mathutils::fraction_inside(bottom_phi, top_phi), theta_criterion);
						
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_y - Vector3s(0, dx * 0.5, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_y + Vector3s(0, dx * 0.5, 0));
                        
                        pressure_grad /= (dx * theta);
						
						scalar term = pressure_grad * dt * (node_vol_y[bucket_idx][i] + node_vol_fluid_y[bucket_idx][i] * node_mfhdvm_hdvm_y[bucket_idx][i]);
						
						bucket_rhs_y(i) -= term;
					}
				}
			}
			
			for(int i = 0; i < num_rhs_z; ++i) {
                const Vector3s& pos_z = bucket_pos_z.segment<3>(i * 3);
                
				if(bucket_weight_z(i) > 0.) {
					const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);

                    scalar near_phi, far_phi;
                    if(indices[0] < 0 || indices[1] < 0) near_phi = scene.interpolateBucketLiquidPhi(pos_z - Vector3s(0.0, 0.0, 0.5 * dx));
                    else near_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) far_phi = scene.interpolateBucketLiquidPhi(pos_z + Vector3s(0.0, 0.0, 0.5 * dx));
                    else far_phi = node_liquid_phi[indices[2]][indices[3]];
                    
					if(near_phi < 0.0 || far_phi < 0.0) {
						scalar theta = 1.;
						if(far_phi >= 0.0 || near_phi >= 0.0)
							theta = std::max(mathutils::fraction_inside(near_phi, far_phi), theta_criterion);
						
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_z - Vector3s(0, 0, dx * 0.5));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_z + Vector3s(0, 0, dx * 0.5));
                        
                        pressure_grad /= (dx * theta);
						
						scalar term = pressure_grad * dt * (node_vol_z[bucket_idx][i] + node_vol_fluid_z[bucket_idx][i] * node_mfhdvm_hdvm_z[bucket_idx][i]);
						
						bucket_rhs_z(i) -= term;
					}
				}
			}
		});
	}
	
	void applySurfTensionFluid( TwoDScene& scene, const std::vector< VectorXs >& surf_tension_vec, std::vector< VectorXs >& rhs_vec_x, std::vector< VectorXs >& rhs_vec_y, std::vector< VectorXs >& rhs_vec_z, const std::vector< VectorXs >& node_vol_x, const std::vector< VectorXs >& node_vol_y, const std::vector< VectorXs >& node_vol_z)
	{
		const Sorter& buckets = scene.getParticleBuckets();
		const scalar dx = scene.getCellSize();
		
		const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
		const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
		const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
		
		const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
		const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
		const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
		const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
		
		const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
		const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
		const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
		
		buckets.for_each_bucket([&] (int bucket_idx) {
			VectorXs& bucket_rhs_x = rhs_vec_x[bucket_idx];
			VectorXs& bucket_rhs_y = rhs_vec_y[bucket_idx];
			VectorXs& bucket_rhs_z = rhs_vec_z[bucket_idx];
			
			const VectorXs& bucket_weight_x = node_solid_weight_x[bucket_idx];
			const VectorXs& bucket_weight_y = node_solid_weight_y[bucket_idx];
			const VectorXs& bucket_weight_z = node_solid_weight_z[bucket_idx];
			
			const int num_rhs_x = bucket_rhs_x.size();
			const int num_rhs_y = bucket_rhs_y.size();
			const int num_rhs_z = bucket_rhs_z.size();
			
			for(int i = 0; i < num_rhs_x; ++i) {
				if(bucket_weight_x(i) > 0.) {
					const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);
					if((indices[0] < 0 || indices[1] < 0) && (indices[2] < 0 || indices[3] < 0)) continue;
					
					const scalar left_phi = (indices[0] < 0 || indices[1] < 0) ? node_liquid_phi[indices[2]][indices[3]] : node_liquid_phi[indices[0]][indices[1]];
					const scalar right_phi = (indices[2] < 0 || indices[3] < 0) ? node_liquid_phi[indices[0]][indices[1]] : node_liquid_phi[indices[2]][indices[3]];
					
					if(left_phi * right_phi <= 0.0) {
						scalar theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
						
                        scalar surf = 0.0;
                        scalar w = 0.0;
                        
                        if(indices[0] >= 0 && indices[1] >= 0) {
                            surf += surf_tension_vec[ indices[0] ][ indices[1] ];
                            w += 1.0;
                        }
                        
                        if(indices[2] >= 0 && indices[3] >= 0) {
                            surf += surf_tension_vec[ indices[2] ][ indices[3] ];
                            w += 1.0;
                        }
                        
                        surf /= w;
						
						scalar term = surf / dx / theta * node_vol_x[bucket_idx][i];
						
						if(left_phi <= 0.0) term *= -1.0;
						
						bucket_rhs_x(i) -= term;
					}
				}
			}
			
			for(int i = 0; i < num_rhs_y; ++i) {
				if(bucket_weight_y(i) > 0.) {
					const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);
					if(indices[0] < 0 || indices[1] < 0 || indices[2] < 0 || indices[3] < 0) continue;
					
					const scalar left_phi = (indices[1] >= node_liquid_phi[indices[0]].size()) ? 0.5 * dx : node_liquid_phi[indices[0]][indices[1]];
					const scalar right_phi = (indices[3] >= node_liquid_phi[indices[2]].size()) ? 0.5 * dx : node_liquid_phi[indices[2]][indices[3]];
					
					if(left_phi * right_phi <= 0.0) {
						scalar theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
						
                        scalar surf = 0.0;
                        scalar w = 0.0;
                        
                        if(indices[0] >= 0 && indices[1] >= 0) {
                            surf += surf_tension_vec[ indices[0] ][ indices[1] ];
                            w += 1.0;
                        }
                        
                        if(indices[2] >= 0 && indices[3] >= 0) {
                            surf += surf_tension_vec[ indices[2] ][ indices[3] ];
                            w += 1.0;
                        }
                        
                        surf /= w;
                        
                        scalar term = surf / dx / theta * node_vol_y[bucket_idx][i];
						
						if(left_phi <= 0.0) term *= -1.0;
						
						bucket_rhs_y(i) -= term;
					}
				}
			}
			
			for(int i = 0; i < num_rhs_z; ++i) {
				if(bucket_weight_z(i) > 0.) {
					const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);
					if(indices[0] < 0 || indices[1] < 0 || indices[2] < 0 || indices[3] < 0) continue;
					
					const scalar left_phi = (indices[1] >= node_liquid_phi[indices[0]].size()) ? 0.5 * dx : node_liquid_phi[indices[0]][indices[1]];
					const scalar right_phi = (indices[3] >= node_liquid_phi[indices[2]].size()) ? 0.5 * dx : node_liquid_phi[indices[2]][indices[3]];
					
					if(left_phi * right_phi <= 0.0) {
						scalar theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
						
                        scalar surf = 0.0;
                        scalar w = 0.0;
                        
                        if(indices[0] >= 0 && indices[1] >= 0) {
                            surf += surf_tension_vec[ indices[0] ][ indices[1] ];
                            w += 1.0;
                        }
                        
                        if(indices[2] >= 0 && indices[3] >= 0) {
                            surf += surf_tension_vec[ indices[2] ][ indices[3] ];
                            w += 1.0;
                        }
                        
                        surf /= w;
                        
                        scalar term = surf / dx / theta * node_vol_z[bucket_idx][i];
						
						if(left_phi <= 0.0) term *= -1.0;
						
						bucket_rhs_z(i) -= term;
					}
				}
			}
		});
	}
	
    void applyPressureGradsFluidRHSBiCGSTAB( const TwoDScene& scene,
                                            const std::vector< VectorXs >& pressure_vec,
                                            std::vector< VectorXs >& rhs_vec_x,
                                            std::vector< VectorXs >& rhs_vec_y,
                                            std::vector< VectorXs >& rhs_vec_z,
                                            const scalar& dt )
    {
        const Sorter& buckets = scene.getParticleBuckets();
        
        const scalar dx = scene.getCellSize();
        const scalar dV = dx * dx * dx;
        
        const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
        
        const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
        const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
        const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
        const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
        
        const std::vector< VectorXs >& node_vol_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_vol_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_vol_z = scene.getNodeFluidVolZ();
        
        const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
        const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
        const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
        
        buckets.for_each_bucket([&] (int bucket_idx) {
            VectorXs& bucket_rhs_x = rhs_vec_x[bucket_idx];
            VectorXs& bucket_rhs_y = rhs_vec_y[bucket_idx];
            VectorXs& bucket_rhs_z = rhs_vec_z[bucket_idx];
            
            const VectorXs& bucket_weight_x = node_solid_weight_x[bucket_idx];
            const VectorXs& bucket_weight_y = node_solid_weight_y[bucket_idx];
            const VectorXs& bucket_weight_z = node_solid_weight_z[bucket_idx];
            
            const VectorXs& bucket_pos_x = scene.getNodePosX(bucket_idx);
            const VectorXs& bucket_pos_y = scene.getNodePosY(bucket_idx);
            const VectorXs& bucket_pos_z = scene.getNodePosZ(bucket_idx);
            
            const int num_rhs_x = bucket_rhs_x.size();
            const int num_rhs_y = bucket_rhs_y.size();
            const int num_rhs_z = bucket_rhs_z.size();
            
            for(int i = 0; i < num_rhs_x; ++i) {
                const Vector3s& pos_x = bucket_pos_x.segment<3>(i * 3);
                
                if(bucket_weight_x(i) > 0.) {
                    const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);
                    
                    scalar left_phi, right_phi;
                    if(indices[0] < 0 || indices[1] < 0) left_phi = scene.interpolateBucketLiquidPhi(pos_x - Vector3s(0.5 * dx, 0.0, 0.0));
                    else left_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) right_phi = scene.interpolateBucketLiquidPhi(pos_x + Vector3s(0.5 * dx, 0.0, 0.0));
                    else right_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(left_phi < 0.0 || right_phi < 0.0) {
                        scalar theta = 1.;
                        if(right_phi >= 0.0 || left_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_x - Vector3s(dx * 0.5, 0, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_x + Vector3s(dx * 0.5, 0, 0));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt;
                        
                        term *= node_vol_x[bucket_idx][i];
                        
                        bucket_rhs_x(i) += term;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_y; ++i) {
                const Vector3s& pos_y = bucket_pos_y.segment<3>(i * 3);
                
                if(bucket_weight_y(i) > 0.) {
                    const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);
                    
                    scalar bottom_phi, top_phi;
                    if(indices[0] < 0 || indices[1] < 0) bottom_phi = scene.interpolateBucketLiquidPhi(pos_y - Vector3s(0.0, 0.5 * dx, 0.0));
                    else bottom_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) top_phi = scene.interpolateBucketLiquidPhi(pos_y + Vector3s(0.0, 0.5 * dx, 0.0));
                    else top_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(bottom_phi < 0.0 || top_phi < 0.0) {
                        scalar theta = 1.;
                        if(top_phi >= 0.0 || bottom_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(bottom_phi, top_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_y - Vector3s(0, dx * 0.5, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_y + Vector3s(0, dx * 0.5, 0));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt;
                        
                        term *= node_vol_y[bucket_idx][i];
                        
                        bucket_rhs_y(i) += term;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_z; ++i) {
                const Vector3s& pos_z = bucket_pos_z.segment<3>(i * 3);
                
                if(bucket_weight_z(i) > 0.) {
                    const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);
                    
                    scalar near_phi, far_phi;
                    if(indices[0] < 0 || indices[1] < 0) near_phi = scene.interpolateBucketLiquidPhi(pos_z - Vector3s(0.0, 0.0, 0.5 * dx));
                    else near_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) far_phi = scene.interpolateBucketLiquidPhi(pos_z + Vector3s(0.0, 0.0, 0.5 * dx));
                    else far_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(near_phi < 0.0 || far_phi < 0.0) {
                        scalar theta = 1.;
                        if(far_phi >= 0.0 || near_phi >= 0.0)
                            theta = std::max(mathutils::fraction_inside(near_phi, far_phi), theta_criterion);
                        
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_z - Vector3s(0, 0, dx * 0.5));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_z + Vector3s(0, 0, dx * 0.5));
                        
                        pressure_grad /= (dx * theta);
                        
                        scalar term = pressure_grad * dt;
                        
                        term *= node_vol_z[bucket_idx][i];
                        
                        bucket_rhs_z(i) += term;
                    }
                }
            }
        });
    }
    
    void identifyValid( TwoDScene& scene )
    {
        const Sorter& buckets = scene.getParticleBuckets();
        
        const scalar dx = scene.getCellSize();
        const scalar dV = dx * dx * dx;
        
        std::vector< VectorXuc >& node_liquid_valid_x = scene.getNodeLiquidValidX();
        std::vector< VectorXuc >& node_liquid_valid_y = scene.getNodeLiquidValidY();
        std::vector< VectorXuc >& node_liquid_valid_z = scene.getNodeLiquidValidZ();
        
        const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
        
        const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
        const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
        const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
        const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
        
        const std::vector< VectorXs >& node_vol_x = scene.getNodeFluidVolX();
        const std::vector< VectorXs >& node_vol_y = scene.getNodeFluidVolY();
        const std::vector< VectorXs >& node_vol_z = scene.getNodeFluidVolZ();
        
        const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
        const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
        const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
        
        buckets.for_each_bucket([&] (int bucket_idx) {
            const VectorXs& bucket_weight_x = node_solid_weight_x[bucket_idx];
            const VectorXs& bucket_weight_y = node_solid_weight_y[bucket_idx];
            const VectorXs& bucket_weight_z = node_solid_weight_z[bucket_idx];
            
            VectorXuc& bucket_valid_x = node_liquid_valid_x[bucket_idx];
            VectorXuc& bucket_valid_y = node_liquid_valid_y[bucket_idx];
            VectorXuc& bucket_valid_z = node_liquid_valid_z[bucket_idx];
            
            const VectorXs& bucket_pos_x = scene.getNodePosX(bucket_idx);
            const VectorXs& bucket_pos_y = scene.getNodePosY(bucket_idx);
            const VectorXs& bucket_pos_z = scene.getNodePosZ(bucket_idx);
            
            const int num_rhs_x = bucket_valid_x.size();
            const int num_rhs_y = bucket_valid_y.size();
            const int num_rhs_z = bucket_valid_z.size();
            
            for(int i = 0; i < num_rhs_x; ++i) {
                bucket_valid_x(i) = 0U;
                const Vector3s& pos_x = bucket_pos_x.segment<3>(i * 3);
                
                if(bucket_weight_x(i) > 0. && node_vol_x[bucket_idx][i] > 0.) {
                    const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);
                    
                    scalar left_phi, right_phi;
                    if(indices[0] < 0 || indices[1] < 0) left_phi = scene.interpolateBucketLiquidPhi(pos_x - Vector3s(0.5 * dx, 0.0, 0.0));
                    else left_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) right_phi = scene.interpolateBucketLiquidPhi(pos_x + Vector3s(0.5 * dx, 0.0, 0.0));
                    else right_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(left_phi < 0.0 || right_phi < 0.0) {
                        bucket_valid_x(i) = 1U;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_y; ++i) {
                bucket_valid_y(i) = 0U;
                const Vector3s& pos_y = bucket_pos_y.segment<3>(i * 3);
                
                if(bucket_weight_y(i) > 0. && node_vol_y[bucket_idx][i] > 0.) {
                    const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);
                    
                    scalar bottom_phi, top_phi;
                    if(indices[0] < 0 || indices[1] < 0) bottom_phi = scene.interpolateBucketLiquidPhi(pos_y - Vector3s(0.0, 0.5 * dx, 0.0));
                    else bottom_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) top_phi = scene.interpolateBucketLiquidPhi(pos_y + Vector3s(0.0, 0.5 * dx, 0.0));
                    else top_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(bottom_phi < 0.0 || top_phi < 0.0) {
                        bucket_valid_y(i) = 1U;
                    }
                }
            }
            
            for(int i = 0; i < num_rhs_z; ++i) {
                bucket_valid_z(i) = 0U;
                const Vector3s& pos_z = bucket_pos_z.segment<3>(i * 3);
                
                if(bucket_weight_z(i) > 0. && node_vol_z[bucket_idx][i] > 0.) {
                    const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);
                    
                    scalar near_phi, far_phi;
                    if(indices[0] < 0 || indices[1] < 0) near_phi = scene.interpolateBucketLiquidPhi(pos_z - Vector3s(0.0, 0.0, 0.5 * dx));
                    else near_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) far_phi = scene.interpolateBucketLiquidPhi(pos_z + Vector3s(0.0, 0.0, 0.5 * dx));
                    else far_phi = node_liquid_phi[indices[2]][indices[3]];
                    
                    if(near_phi < 0.0 || far_phi < 0.0) {
                        bucket_valid_z(i) = 1U;
                    }
                }
            }
        });
    }
    
	void applyPressureGradsFluid( TwoDScene& scene, const std::vector< VectorXs >& pressure_vec, std::vector< VectorXs >& rhs_vec_x, std::vector< VectorXs >& rhs_vec_y, std::vector< VectorXs >& rhs_vec_z, const std::vector< VectorXs >& node_inv_mdv_x, const std::vector< VectorXs >& node_inv_mdv_y, const std::vector< VectorXs >& node_inv_mdv_z, const scalar& dt )
	{
		const Sorter& buckets = scene.getParticleBuckets();
		
		const scalar dx = scene.getCellSize();
		const scalar dV = dx * dx * dx;
		
		std::vector< VectorXuc >& node_liquid_valid_x = scene.getNodeLiquidValidX();
		std::vector< VectorXuc >& node_liquid_valid_y = scene.getNodeLiquidValidY();
		std::vector< VectorXuc >& node_liquid_valid_z = scene.getNodeLiquidValidZ();
		
		const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
		const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
		const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
		
		const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
		const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
		const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
		const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
		
		const std::vector< VectorXs >& node_vol_x = scene.getNodeFluidVolX();
		const std::vector< VectorXs >& node_vol_y = scene.getNodeFluidVolY();
		const std::vector< VectorXs >& node_vol_z = scene.getNodeFluidVolZ();
		
		const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
		const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
		const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
	
		buckets.for_each_bucket([&] (int bucket_idx) {
			VectorXs& bucket_rhs_x = rhs_vec_x[bucket_idx];
			VectorXs& bucket_rhs_y = rhs_vec_y[bucket_idx];
			VectorXs& bucket_rhs_z = rhs_vec_z[bucket_idx];
			
			const VectorXs& bucket_weight_x = node_solid_weight_x[bucket_idx];
			const VectorXs& bucket_weight_y = node_solid_weight_y[bucket_idx];
			const VectorXs& bucket_weight_z = node_solid_weight_z[bucket_idx];
			
			VectorXuc& bucket_valid_x = node_liquid_valid_x[bucket_idx];
			VectorXuc& bucket_valid_y = node_liquid_valid_y[bucket_idx];
			VectorXuc& bucket_valid_z = node_liquid_valid_z[bucket_idx];
            
            const VectorXs& bucket_pos_x = scene.getNodePosX(bucket_idx);
            const VectorXs& bucket_pos_y = scene.getNodePosY(bucket_idx);
            const VectorXs& bucket_pos_z = scene.getNodePosZ(bucket_idx);
			
			const int num_rhs_x = bucket_rhs_x.size();
			const int num_rhs_y = bucket_rhs_y.size();
			const int num_rhs_z = bucket_rhs_z.size();
			
			for(int i = 0; i < num_rhs_x; ++i) {
				bucket_valid_x(i) = 0U;
                const Vector3s& pos_x = bucket_pos_x.segment<3>(i * 3);
                
				if(bucket_weight_x(i) > 0.) {
					const Vector4i& indices = node_pressure_index_x[bucket_idx].segment<4>(i * 4);

                    scalar left_phi, right_phi;
                    if(indices[0] < 0 || indices[1] < 0) left_phi = scene.interpolateBucketLiquidPhi(pos_x - Vector3s(0.5 * dx, 0.0, 0.0));
                    else left_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) right_phi = scene.interpolateBucketLiquidPhi(pos_x + Vector3s(0.5 * dx, 0.0, 0.0));
                    else right_phi = node_liquid_phi[indices[2]][indices[3]];
                    
					if(left_phi < 0.0 || right_phi < 0.0) {
						scalar theta = 1.;
						if(right_phi >= 0.0 || left_phi >= 0.0)
							theta = std::max(mathutils::fraction_inside(left_phi, right_phi), theta_criterion);
						
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_x - Vector3s(dx * 0.5, 0, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_x + Vector3s(dx * 0.5, 0, 0));
                        
                        pressure_grad /= (dx * theta);
                        
						scalar term = pressure_grad * dt;
						
						term *= node_vol_x[bucket_idx][i] * node_inv_mdv_x[bucket_idx](i);
						
						bucket_rhs_x(i) -= term;
						bucket_valid_x(i) = 1U;
					}
				}
			}
			
			for(int i = 0; i < num_rhs_y; ++i) {
				bucket_valid_y(i) = 0U;
				const Vector3s& pos_y = bucket_pos_y.segment<3>(i * 3);
                
				if(bucket_weight_y(i) > 0.) {
					const Vector4i& indices = node_pressure_index_y[bucket_idx].segment<4>(i * 4);
					
                    scalar bottom_phi, top_phi;
                    if(indices[0] < 0 || indices[1] < 0) bottom_phi = scene.interpolateBucketLiquidPhi(pos_y - Vector3s(0.0, 0.5 * dx, 0.0));
                    else bottom_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) top_phi = scene.interpolateBucketLiquidPhi(pos_y + Vector3s(0.0, 0.5 * dx, 0.0));
                    else top_phi = node_liquid_phi[indices[2]][indices[3]];
                    
					if(bottom_phi < 0.0 || top_phi < 0.0) {
						scalar theta = 1.;
						if(top_phi >= 0.0 || bottom_phi >= 0.0)
							theta = std::max(mathutils::fraction_inside(bottom_phi, top_phi), theta_criterion);
						
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_y - Vector3s(0, dx * 0.5, 0));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_y + Vector3s(0, dx * 0.5, 0));
                        
                        pressure_grad /= (dx * theta);
                        
						scalar term = pressure_grad * dt;
						
						term *= node_vol_y[bucket_idx][i] * node_inv_mdv_y[bucket_idx](i);
						
						bucket_rhs_y(i) -= term;
						bucket_valid_y(i) = 1U;
					}
				}
			}
			
			for(int i = 0; i < num_rhs_z; ++i) {
				bucket_valid_z(i) = 0U;
				const Vector3s& pos_z = bucket_pos_z.segment<3>(i * 3);
                
				if(bucket_weight_z(i) > 0.) {
					const Vector4i& indices = node_pressure_index_z[bucket_idx].segment<4>(i * 4);
                    
                    scalar near_phi, far_phi;
                    if(indices[0] < 0 || indices[1] < 0) near_phi = scene.interpolateBucketLiquidPhi(pos_z - Vector3s(0.0, 0.0, 0.5 * dx));
                    else near_phi = node_liquid_phi[indices[0]][indices[1]];
                    
                    if(indices[2] < 0 || indices[3] < 0) far_phi = scene.interpolateBucketLiquidPhi(pos_z + Vector3s(0.0, 0.0, 0.5 * dx));
                    else far_phi = node_liquid_phi[indices[2]][indices[3]];
                    
					if(near_phi < 0.0 || far_phi < 0.0) {
						scalar theta = 1.;
						if(far_phi >= 0.0 || near_phi >= 0.0)
							theta = std::max(mathutils::fraction_inside(near_phi, far_phi), theta_criterion);
						
                        scalar pressure_grad = 0;
                        if(indices[0] >= 0 && indices[1] >= 0) pressure_grad -= pressure_vec[indices[0]][indices[1]];
                        else pressure_grad -= scene.interpolateBucketPressure(pos_z - Vector3s(0, 0, dx * 0.5));
                        
                        if(indices[2] >= 0 && indices[3] >= 0) pressure_grad += pressure_vec[indices[2]][indices[3]];
                        else pressure_grad += scene.interpolateBucketPressure(pos_z + Vector3s(0, 0, dx * 0.5));
                        
                        pressure_grad /= (dx * theta);
                        
						scalar term = pressure_grad * dt;
						
						term *= node_vol_z[bucket_idx][i] * node_inv_mdv_z[bucket_idx](i);
						
						bucket_rhs_z(i) -= term;
						bucket_valid_z(i) = 1U;
					}
				}
			}
		});
	}
	
    void extrapolateCoarse( const TwoDScene& scene, Array3uc& valid_x, Array3uc& valid_y, Array3uc& valid_z,
                           Array3s& vel_x, Array3s& vel_y, Array3s& vel_z,
                           Array3s& fluid_pg_x,
                           Array3s& fluid_pg_y,
                           Array3s& fluid_pg_z)
    {
        const Sorter& buckets = scene.getParticleBuckets();
        const Vector3i offsets[] = {
            Vector3i(-1, 0, 0),
            Vector3i(1, 0, 0),
            Vector3i(0, -1, 0),
            Vector3i(0, 1, 0),
            Vector3i(0, 0, -1),
            Vector3i(0, 0, 1)
        };
        
        Array3s vel_copy_x;
        Array3uc valid_copy_x;
        Array3s vel_copy_y;
        Array3uc valid_copy_y;
        Array3s vel_copy_z;
        Array3uc valid_copy_z;
        Array3s fluid_pg_copy_x;
        Array3s fluid_pg_copy_y;
        Array3s fluid_pg_copy_z;
        
        for(int layers = 0; layers < 2; ++layers) {
            vel_copy_x = vel_x;
            valid_copy_x = valid_x;
            fluid_pg_copy_x = fluid_pg_x;
            
            buckets.for_each_bucket_x([&] (int bucket_idx) {
                Vector3i handle = buckets.bucket_handle_x(bucket_idx);
                
                if(valid_copy_x(handle)) return;
                
                scalar sum = 0;
                scalar sum_fpg = 0;
                int count = 0;
                
                for(int r = 0; r < 6; ++r)
                {
                    Vector3i neigh_idx = handle + offsets[r];
                    if(!buckets.has_bucket_x(neigh_idx) || !valid_copy_x(neigh_idx)) continue;
                    
                    sum += vel_copy_x(neigh_idx);
                    sum_fpg += fluid_pg_copy_x(neigh_idx);
                    count++;
                }
                
                if(count) {
                    vel_x(handle) = sum / (scalar) count;
                    fluid_pg_x(handle) = sum_fpg / (scalar) count;
                    valid_x(handle) = 1U;
                }
            });
            
            vel_copy_y = vel_y;
            valid_copy_y = valid_y;
            fluid_pg_copy_y = fluid_pg_y;
            
            buckets.for_each_bucket_y([&] (int bucket_idx) {
                Vector3i handle = buckets.bucket_handle_y(bucket_idx);
                
                if(valid_copy_y(handle)) return;
                
                scalar sum = 0;
                scalar sum_fpg = 0;
                int count = 0;
                
                for(int r = 0; r < 6; ++r)
                {
                    Vector3i neigh_idx = handle + offsets[r];
                    if(!buckets.has_bucket_y(neigh_idx) || !valid_copy_y(neigh_idx)) continue;
                    
                    sum += vel_copy_y(neigh_idx);
                    sum_fpg += fluid_pg_copy_y(neigh_idx);
                    count++;
                }
                
                if(count) {
                    vel_y(handle) = sum / (scalar) count;
                    fluid_pg_y(handle) = sum_fpg / (scalar) count;
                    valid_y(handle) = 1U;
                }
            });
            
            vel_copy_z = vel_z;
            valid_copy_z = valid_z;
            fluid_pg_copy_z = fluid_pg_z;
            
            buckets.for_each_bucket_z([&] (int bucket_idx) {
                Vector3i handle = buckets.bucket_handle_z(bucket_idx);
                
                if(valid_copy_z(handle)) return;
                
                scalar sum = 0;
                scalar sum_fpg = 0;
                int count = 0;
                
                for(int r = 0; r < 6; ++r)
                {
                    Vector3i neigh_idx = handle + offsets[r];
                    if(!buckets.has_bucket_z(neigh_idx) || !valid_copy_z(neigh_idx)) continue;
                    
                    sum += vel_copy_z(neigh_idx);
                    sum_fpg += fluid_pg_copy_z(neigh_idx);
                    count++;
                }
                
                if(count) {
                    vel_z(handle) = sum / (scalar) count;
                    fluid_pg_z(handle) = sum_fpg / (scalar) count;
                    valid_z(handle) = 1U;
                }
            });
        }
    }
    
	void extrapolate( const TwoDScene& scene, std::vector< VectorXuc >& node_valid, const std::vector< VectorXi >& node_cpidx,
					 std::vector< VectorXs >& node_vel)
	{
		const Sorter& buckets = scene.getParticleBuckets();
		const int num_nodes = scene.getDefaultNumNodes();
		
		const Vector3i offsets[] = {
			Vector3i(-1, 0, 0),
			Vector3i(1, 0, 0),
			Vector3i(0, -1, 0),
			Vector3i(0, 1, 0),
			Vector3i(0, 0, -1),
			Vector3i(0, 0, 1)
		};
		
		for(int layers = 0; layers < 2; ++layers) {
			std::vector< VectorXs > node_vel_copy = node_vel;
			std::vector< VectorXuc > node_valid_copy = node_valid;
			
			buckets.for_each_bucket([&] (int bucket_idx) {
				if(node_cpidx[bucket_idx].size() == 0) return;
				
				const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
				
				const VectorXuc& bucket_valid = node_valid_copy[bucket_idx];
				const VectorXi& bucket_cpidx = node_cpidx[bucket_idx];
				
				for(int k = 0; k < num_nodes; ++k) for(int j = 0; j < num_nodes; ++j) for(int i = 0; i < num_nodes; ++i)
				{
					const int center_cell_idx = k * num_nodes * num_nodes + j * num_nodes + i;
					
					// ignore empty center cells
					if(bucket_cpidx[center_cell_idx] < 0) continue;
					
					// ignore valid center cells
					const int mapped_center_idx = bucket_cpidx[center_cell_idx];
					if(bucket_valid[mapped_center_idx]) continue;
					
					scalar sum = 0;
					int count = 0;
				
					for(int r = 0; r < 6; ++r) {
						Vector3i neigh_idx = Vector3i(i, j, k) + offsets[r];
						
						Vector3i node_bucket_handle = bucket_handle;
						for(int r = 0; r < 3; ++r) {
							if(neigh_idx(r) < 0) {
								node_bucket_handle(r)--;
								neigh_idx(r) += num_nodes;
							}
							if(neigh_idx(r) >= num_nodes) {
								node_bucket_handle(r)++;
								neigh_idx(r) -= num_nodes;
							}
							
							assert( neigh_idx(r) >= 0 && neigh_idx(r) < num_nodes );
						}
						
						if( node_bucket_handle(0) < 0 || node_bucket_handle(0) >= buckets.dim_size(0) ||
						    node_bucket_handle(1) < 0 || node_bucket_handle(1) >= buckets.dim_size(1) ||
						    node_bucket_handle(2) < 0 || node_bucket_handle(2) >= buckets.dim_size(2)
						   ) continue;
						
						const int node_bucket_idx = buckets.bucket_index(node_bucket_handle);
						if(node_cpidx[node_bucket_idx].size() == 0) continue;
						
						const int neigh_cell_idx = neigh_idx(2) * num_nodes * num_nodes + neigh_idx(1) * num_nodes + neigh_idx(0);
						
						const int mapped_neigh_idx = node_cpidx[node_bucket_idx][neigh_cell_idx];
						
						if(mapped_neigh_idx < 0) continue;
						
						if(!node_valid_copy[node_bucket_idx][mapped_neigh_idx]) continue;
						
						sum += node_vel_copy[node_bucket_idx][mapped_neigh_idx];
						++count;
					}
					
					if(count > 0) {
						node_vel[bucket_idx][mapped_center_idx] = sum / (scalar) count;
						node_valid[bucket_idx][mapped_center_idx] = 1U;
					}
				}
			});
		}
	}
	
	void constructJacobiPreconditioner( const TwoDScene& scene, std::vector< VectorXs >& out_node_vec, const std::vector< VectorXs >& node_inv_mdv_x, const std::vector< VectorXs >& node_inv_mdv_y, const std::vector< VectorXs >& node_inv_mdv_z, const std::vector< VectorXs >& node_inv_mdvs_x, const std::vector< VectorXs >& node_inv_mdvs_y, const std::vector< VectorXs >& node_inv_mdvs_z, const scalar& dt )
	{
		const Sorter& buckets = scene.getParticleBuckets();
		
		const std::vector< VectorXs >& node_liquid_phi = scene.getNodeLiquidPhi();
		const std::vector< VectorXi >& pressure_neighbors = scene.getPressureNeighbors();
		const std::vector< VectorXs >& node_solid_weight_x = scene.getNodeSolidWeightX();
		const std::vector< VectorXs >& node_solid_weight_y = scene.getNodeSolidWeightY();
		const std::vector< VectorXs >& node_solid_weight_z = scene.getNodeSolidWeightZ();
		const std::vector< VectorXi >& node_pressure_index_x = scene.getNodePressureIndexX();
		const std::vector< VectorXi >& node_pressure_index_y = scene.getNodePressureIndexY();
		const std::vector< VectorXi >& node_pressure_index_z = scene.getNodePressureIndexZ();
		const std::vector< VectorXs >& node_fluid_vol_x = scene.getNodeFluidVolX();
		const std::vector< VectorXs >& node_fluid_vol_y = scene.getNodeFluidVolY();
		const std::vector< VectorXs >& node_fluid_vol_z = scene.getNodeFluidVolZ();
		const std::vector< VectorXs >& node_elasto_vol_x = scene.getNodeVolX();
		const std::vector< VectorXs >& node_elasto_vol_y = scene.getNodeVolY();
		const std::vector< VectorXs >& node_elasto_vol_z = scene.getNodeVolZ();
		const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
		const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
		const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();
		
		const scalar dx = scene.getCellSize();
		const scalar coeff = dt / (dx * dx);
		const scalar dV = dx * dx * dx;
		
		buckets.for_each_bucket([&] (int bucket_idx) {
			const VectorXs& bucket_liquid_phi = node_liquid_phi[bucket_idx];
			const VectorXi& bucket_pn = pressure_neighbors[bucket_idx];;
            const VectorXs& bucket_pos_p = scene.getNodePosP(bucket_idx);
            
			const int num_pressure = out_node_vec[bucket_idx].size();
			
			assert(out_node_vec[bucket_idx].size() == bucket_liquid_phi.size());
			
			VectorXs& bucket_out = out_node_vec[bucket_idx];
			
			for(int i = 0; i < num_pressure; ++i)
			{
				const scalar center_phi = bucket_liquid_phi(i);
				
				if(center_phi > 0.0) {
					bucket_out(i) = 1.0;
					continue;
				}
                
                const Vector3s center_pos = bucket_pos_p.segment<3>(i * 3);
                
				scalar diag_term = 0.0;
				
				const int bucket_idx_left = bucket_pn[i * 12 + 0];
				const int node_idx_left = bucket_pn[i * 12 + 1];
				
				if(bucket_idx_left >= 0 && node_idx_left >= 0) {
					const scalar& w = node_solid_weight_x[bucket_idx_left][node_idx_left];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_x[bucket_idx_left][node_idx_left] * (1.0 - node_psi_x[bucket_idx_left][node_idx_left]) * node_inv_mdv_x[bucket_idx_left][node_idx_left] + node_elasto_vol_x[bucket_idx_left][node_idx_left] * node_psi_x[bucket_idx_left][node_idx_left] * node_inv_mdvs_x[bucket_idx_left][node_idx_left];
					
					const int bucket_pressure_left = node_pressure_index_x[bucket_idx_left][node_idx_left * 4 + 0];
					const int node_pressure_left = node_pressure_index_x[bucket_idx_left][node_idx_left * 4 + 1];
					
					if(bucket_pressure_left >= 0 && node_pressure_left >= 0) {
						const scalar left_phi = node_liquid_phi[bucket_pressure_left][node_pressure_left];
						
						if(left_phi < 0.0) {
							diag_term += term;
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
							diag_term += term / theta;
						}
                    } else {
                        diag_term += term;
                        
                        const scalar left_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( dx, 0., 0. ));
                        
                        if(left_phi < 0.0) {
                            diag_term += term;
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
                            diag_term += term / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.5 * dx, 0.0, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar left_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( dx, 0., 0. ));
                    
                    if(left_phi < 0.0) {
                        diag_term += term;
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, left_phi), theta_criterion);
                        diag_term += term / theta;
                    }
                }
				
				const int bucket_idx_right = bucket_pn[i * 12 + 2];
				const int node_idx_right = bucket_pn[i * 12 + 3];
				
				if(bucket_idx_right >= 0 && node_idx_right >= 0) {
					const scalar& w = node_solid_weight_x[bucket_idx_right][node_idx_right];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_x[bucket_idx_right][node_idx_right] * (1.0 - node_psi_x[bucket_idx_right][node_idx_right]) * node_inv_mdv_x[bucket_idx_right][node_idx_right] + node_elasto_vol_x[bucket_idx_right][node_idx_right] * node_psi_x[bucket_idx_right][node_idx_right] * node_inv_mdvs_x[bucket_idx_right][node_idx_right];
					
					const int bucket_pressure_right = node_pressure_index_x[bucket_idx_right][node_idx_right * 4 + 2];
					const int node_pressure_right = node_pressure_index_x[bucket_idx_right][node_idx_right * 4 + 3];
					
					if(bucket_pressure_right >= 0 && node_pressure_right >= 0) {
						const scalar right_phi = node_liquid_phi[bucket_pressure_right][node_pressure_right];
						
						if(right_phi < 0.0) {
							diag_term += term;
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
							diag_term += term / theta;
						}
                    } else {
                        const scalar right_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( dx, 0., 0. ));
                        
                        if(right_phi < 0.0) {
                            diag_term += term;
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
                            diag_term += term / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.5 * dx, 0.0, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightX(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar right_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( dx, 0., 0. ));
                    
                    if(right_phi < 0.0) {
                        diag_term += term;
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, right_phi), theta_criterion);
                        diag_term += term / theta;
                    }
                }
				
				const int bucket_idx_bottom = bucket_pn[i * 12 + 4];
				const int node_idx_bottom = bucket_pn[i * 12 + 5];
				
				if(bucket_idx_bottom >= 0 && node_idx_bottom >= 0) {
					const scalar& w = node_solid_weight_y[bucket_idx_bottom][node_idx_bottom];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_y[bucket_idx_bottom][node_idx_bottom] * (1.0 - node_psi_y[bucket_idx_bottom][node_idx_bottom]) * node_inv_mdv_y[bucket_idx_bottom][node_idx_bottom] + node_elasto_vol_y[bucket_idx_bottom][node_idx_bottom] * node_psi_y[bucket_idx_bottom][node_idx_bottom] * node_inv_mdvs_y[bucket_idx_bottom][node_idx_bottom];
					
					const int bucket_pressure_bottom = node_pressure_index_y[bucket_idx_bottom][node_idx_bottom * 4 + 0];
					const int node_pressure_bottom = node_pressure_index_y[bucket_idx_bottom][node_idx_bottom * 4 + 1];
					
					if(bucket_pressure_bottom >= 0 && node_pressure_bottom >= 0) {
						const scalar bottom_phi = node_liquid_phi[bucket_pressure_bottom][node_pressure_bottom];
						
						if(bottom_phi < 0.0) {
							diag_term += term;
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
							diag_term += term / theta;
						}
                    } else {
                        const scalar bottom_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., dx, 0. ));
                        
                        if(bottom_phi < 0.0) {
                            diag_term += term;
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
                            diag_term += term / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.0, 0.5 * dx, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar bottom_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., dx, 0. ));
                    
                    if(bottom_phi < 0.0) {
                        diag_term += term;
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, bottom_phi), theta_criterion);
                        diag_term += term / theta;
                    }
                }
				
				const int bucket_idx_top = bucket_pn[i * 12 + 6];
				const int node_idx_top = bucket_pn[i * 12 + 7];
				
				if(bucket_idx_top >= 0 && node_idx_top >= 0) {
					const scalar& w = node_solid_weight_y[bucket_idx_top][node_idx_top];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_y[bucket_idx_top][node_idx_top] * (1.0 - node_psi_y[bucket_idx_top][node_idx_top]) * node_inv_mdv_y[bucket_idx_top][node_idx_top] + node_elasto_vol_y[bucket_idx_top][node_idx_top] * node_psi_y[bucket_idx_top][node_idx_top] * node_inv_mdvs_y[bucket_idx_top][node_idx_top];
					
					const int bucket_pressure_top = node_pressure_index_y[bucket_idx_top][node_idx_top * 4 + 2];
					const int node_pressure_top = node_pressure_index_y[bucket_idx_top][node_idx_top * 4 + 3];
					
					if(bucket_pressure_top >= 0 && node_pressure_top >= 0) {
						const scalar top_phi = node_liquid_phi[bucket_pressure_top][node_pressure_top];
						
						if(top_phi < 0.0) {
							diag_term += term;
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
							diag_term += term;
						}
                    } else {
                        const scalar top_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., dx, 0. ));
                        
                        if(top_phi < 0.0) {
                            diag_term += term;
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
                            diag_term += term / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.0, 0.5 * dx, 0.0);
                    const scalar w = scene.interpolateBucketSolidWeightY(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar top_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., dx, 0. ));
                    
                    if(top_phi < 0.0) {
                        diag_term += term;
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, top_phi), theta_criterion);
                        diag_term += term / theta;
                    }
                }
				
				const int bucket_idx_near = bucket_pn[i * 12 + 8];
				const int node_idx_near = bucket_pn[i * 12 + 9];
				
				if(bucket_idx_near >= 0 && node_idx_near >= 0) {
					const scalar& w = node_solid_weight_z[bucket_idx_near][node_idx_near];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_z[bucket_idx_near][node_idx_near] * (1.0 - node_psi_z[bucket_idx_near][node_idx_near]) * node_inv_mdv_z[bucket_idx_near][node_idx_near] + node_elasto_vol_z[bucket_idx_near][node_idx_near] * node_psi_z[bucket_idx_near][node_idx_near] * node_inv_mdvs_z[bucket_idx_near][node_idx_near];
					
					const int bucket_pressure_near = node_pressure_index_z[bucket_idx_near][node_idx_near * 4 + 0];
					const int node_pressure_near = node_pressure_index_z[bucket_idx_near][node_idx_near * 4 + 1];
					
					if(bucket_pressure_near >= 0 && node_pressure_near >= 0) {
						const scalar near_phi = node_liquid_phi[bucket_pressure_near][node_pressure_near];
						
						if(near_phi < 0.0) {
							diag_term += term;
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
							diag_term += term / theta;
						}
                    } else {
                        const scalar near_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., 0., dx ));
                        
                        if(near_phi < 0.0) {
                            diag_term += term;
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
                            diag_term += term / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np - Vector3s(0.0, 0.0, 0.5 * dx);
                    const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar near_phi = scene.interpolateBucketLiquidPhi(center_pos - Vector3s( 0., 0., dx ));
                    
                    if(near_phi < 0.0) {
                        diag_term += term;
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, near_phi), theta_criterion);
                        diag_term += term / theta;
                    }
                }
				
				const int bucket_idx_far = bucket_pn[i * 12 + 10];
				const int node_idx_far = bucket_pn[i * 12 + 11];
				
				if(bucket_idx_far >= 0 && node_idx_far >= 0) {
					const scalar& w = node_solid_weight_z[bucket_idx_far][node_idx_far];
					scalar term = w * coeff;
					
					term *= node_fluid_vol_z[bucket_idx_far][node_idx_far] * (1.0 - node_psi_z[bucket_idx_far][node_idx_far]) * node_inv_mdv_z[bucket_idx_far][node_idx_far] + node_elasto_vol_z[bucket_idx_far][node_idx_far] * node_psi_z[bucket_idx_far][node_idx_far] * node_inv_mdvs_z[bucket_idx_far][node_idx_far];
					
					const int bucket_pressure_far = node_pressure_index_z[bucket_idx_far][node_idx_far * 4 + 2];
					const int node_pressure_far = node_pressure_index_z[bucket_idx_far][node_idx_far * 4 + 3];
					
					if(bucket_pressure_far >= 0 && node_pressure_far >= 0) {
						const scalar far_phi = node_liquid_phi[bucket_pressure_far][node_pressure_far];
						
						if(far_phi < 0.0) {
							diag_term += term;
						} else {
							const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
							diag_term += term / theta;
						}
                    } else {
                        const scalar far_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., 0., dx ));
                        
                        if(far_phi < 0.0) {
                            diag_term += term;
                        } else {
                            const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
                            diag_term += term / theta;
                        }
                    }
                } else {
                    const Vector3s& np = bucket_pos_p.segment<3>(i * 3);
                    Vector3s nnp = np + Vector3s(0.0, 0.0, 0.5 * dx);
                    const scalar w = scene.interpolateBucketSolidWeightZ(nnp);
                    scalar term = w * coeff / scene.getLiquidInfo().liquid_density; // coarse grid has no solid
                    
                    const scalar far_phi = scene.interpolateBucketLiquidPhi(center_pos + Vector3s( 0., 0., dx ));
                    
                    if(far_phi < 0.0) {
                        diag_term += term;
                    } else {
                        const scalar theta = std::max(mathutils::fraction_inside(center_phi, far_phi), theta_criterion);
                        diag_term += term / theta;
                    }
                }
				
				if(diag_term == 0.0) {
					bucket_out(i) = 1.0;
				} else {
					bucket_out(i) = 1.0 / diag_term;
				}
			}
		});
	}
};

