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


#include "Viscosity.h"
#include "TwoDScene.h"
#include "ThreadUtils.h"
#include "MathUtilities.h"
#include "AlgebraicMultigrid.h"

namespace viscosity {
    void applyNodeViscositySolidRHS( const TwoDScene& scene,
                                    const std::vector< VectorXs >& node_fv_0_x,
                                    const std::vector< VectorXs >& node_fv_0_y,
                                    const std::vector< VectorXs >& node_fv_0_z,
                                    const std::vector< VectorXs >& node_fv_1_x,
                                    const std::vector< VectorXs >& node_fv_1_y,
                                    const std::vector< VectorXs >& node_fv_1_z,
                                    std::vector< VectorXs >& node_rhs_x,
                                    std::vector< VectorXs >& node_rhs_y,
                                    std::vector< VectorXs >& node_rhs_z)
    {
        const std::vector< VectorXs >& node_vol_x = scene.getNodeVolX();
        const std::vector< VectorXs >& node_vol_y = scene.getNodeVolY();
        const std::vector< VectorXs >& node_vol_z = scene.getNodeVolZ();
        
        const Sorter& buckets = scene.getParticleBuckets();
        const scalar rho = scene.getLiquidInfo().liquid_density;
        
        buckets.for_each_bucket([&] (int bucket_idx) {
            node_rhs_x[bucket_idx] += VectorXs((node_fv_1_x[bucket_idx] - node_fv_0_x[bucket_idx]).array() * node_vol_x[bucket_idx].array()) * rho;
            node_rhs_y[bucket_idx] += VectorXs((node_fv_1_y[bucket_idx] - node_fv_0_y[bucket_idx]).array() * node_vol_y[bucket_idx].array()) * rho;
            node_rhs_z[bucket_idx] += VectorXs((node_fv_1_z[bucket_idx] - node_fv_0_z[bucket_idx]).array() * node_vol_z[bucket_idx].array()) * rho;
        });
    }
    
    void applyNodeViscosityExplicit( const TwoDScene& scene,
                                    const std::vector< VectorXs >& node_vel_src_x,
                                    const std::vector< VectorXs >& node_vel_src_y,
                                    const std::vector< VectorXs >& node_vel_src_z,
                                    std::vector< VectorXs >& node_vel_x,
                                    std::vector< VectorXs >& node_vel_y,
                                    std::vector< VectorXs >& node_vel_z,
                                    const scalar& dt )
    {
        const std::vector< VectorXi >& node_cpidx_x = scene.getNodeCompressedIndexX();
        const std::vector< VectorXi >& node_cpidx_y = scene.getNodeCompressedIndexY();
        const std::vector< VectorXi >& node_cpidx_z = scene.getNodeCompressedIndexZ();
        const std::vector< VectorXi >& node_cpidx_p = scene.getNodeCompressedIndexP();
        
        const std::vector< VectorXi >& node_indices_x = scene.getNodeIndicesX();
        const std::vector< VectorXi >& node_indices_y = scene.getNodeIndicesY();
        const std::vector< VectorXi >& node_indices_z = scene.getNodeIndicesZ();
        
        const std::vector< VectorXi >& node_index_ex = scene.getNodeIndexEdgeX();
        const std::vector< VectorXi >& node_index_ey = scene.getNodeIndexEdgeY();
        const std::vector< VectorXi >& node_index_ez = scene.getNodeIndexEdgeZ();
        
        const std::vector< VectorXi >& node_index_p_x = scene.getNodePressureIndexX();
        const std::vector< VectorXi >& node_index_p_y = scene.getNodePressureIndexY();
        const std::vector< VectorXi >& node_index_p_z = scene.getNodePressureIndexZ();
        
        const std::vector< VectorXs >& node_liquid_c_vf = scene.getNodeLiquidVolFracCentre();
        const std::vector< VectorXs >& node_liquid_u_vf = scene.getNodeLiquidVolFracU();
        const std::vector< VectorXs >& node_liquid_v_vf = scene.getNodeLiquidVolFracV();
        const std::vector< VectorXs >& node_liquid_w_vf = scene.getNodeLiquidVolFracW();
        const std::vector< VectorXs >& node_liquid_ex_vf = scene.getNodeLiquidVolFracEX();
        const std::vector< VectorXs >& node_liquid_ey_vf = scene.getNodeLiquidVolFracEY();
        const std::vector< VectorXs >& node_liquid_ez_vf = scene.getNodeLiquidVolFracEZ();
        
        const int num_nodes = scene.getDefaultNumNodes();
        const Sorter& buckets = scene.getParticleBuckets();
        
        auto get_value = [&] (const std::vector< VectorXi >& node_cpidx,
                             const std::vector< VectorXs >& data,
                             const Vector3i& bucket_handle,
                             const Vector3i& node_handle,
                             const scalar& default_val) -> scalar
        {
            Vector3i cur_bucket_handle = bucket_handle;
            Vector3i cur_node_handle = node_handle;
            
            // check boundary
            for(int r = 0; r < 3; ++r) {
                while(cur_node_handle(r) < 0) {
                    cur_node_handle(r) += num_nodes;
                    --cur_bucket_handle(r);
                }
                
                while(cur_node_handle(r) >= num_nodes) {
                    cur_node_handle(r) -= num_nodes;
                    ++cur_bucket_handle(r);
                }
            }
            
            if(cur_bucket_handle(0) < 0 || cur_bucket_handle(0) >= buckets.ni ||
               cur_bucket_handle(1) < 0 || cur_bucket_handle(1) >= buckets.nj ||
               cur_bucket_handle(2) < 0 || cur_bucket_handle(2) >= buckets.nk )
            {
                return default_val;
            }
            
            const int bucket_idx = buckets.bucket_index(cur_bucket_handle);
            if(data[bucket_idx].size() == 0) return default_val;
            
            const int node_idx = cur_node_handle(2) * num_nodes * num_nodes + cur_node_handle(1) * num_nodes + cur_node_handle(0);
            const int mapped_idx = node_cpidx[bucket_idx][node_idx];
            if(mapped_idx == -1) return default_val;
            
            return data[bucket_idx][mapped_idx];
        };
        
        auto get_value_fast = [&] (const std::vector< VectorXi >& node_index_p,
                                     const std::vector< VectorXs >& data,
                                     int bucket_idx, int node_idx, int num_sel, int sel, const scalar& default_val) -> scalar
        {
            if(node_index_p[bucket_idx].size() == 0)
                return default_val;
            
            const int np_bucket_idx = node_index_p[bucket_idx]((node_idx * num_sel + sel) * 2 + 0);
            const int mapped_idx = node_index_p[bucket_idx]((node_idx * num_sel + sel) * 2 + 1);
            
            if(np_bucket_idx == -1 || mapped_idx == -1)
                return default_val;
            
            return data[np_bucket_idx][mapped_idx];
        };
        
        const scalar dx = scene.getCellSize();
        const scalar coeff = dt / (dx * dx) * scene.getLiquidInfo().viscosity / scene.getLiquidInfo().liquid_density;
        
        buckets.for_each_bucket([&] (int bucket_idx) {
            const int num_node_x = node_indices_x[bucket_idx].size() / 3;
            const int num_node_y = node_indices_y[bucket_idx].size() / 3;
            const int num_node_z = node_indices_z[bucket_idx].size() / 3;
            
            const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
            
            for(int i = 0; i < num_node_x; ++i) {
                const Vector3i& node_handle = node_indices_x[bucket_idx].segment<3>(i * 3);
                const scalar factor = coeff;
                
                const scalar vol_left = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
                const scalar vol_right = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
                
                const scalar vol_back = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, i, 4, 0, 0.0);
                const scalar vol_front = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, i, 4, 1, 0.0);
                
                const scalar vol_bottom = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, i, 4, 2, 0.0);
                const scalar vol_top = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, i, 4, 3, 0.0);
                
                const scalar centre_vel = node_vel_src_x[bucket_idx](i);
                scalar vel = centre_vel;
                
                vel += 2.0 * factor * vol_right * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), centre_vel) - centre_vel);
                vel += 2.0 * factor * vol_left * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), centre_vel) - centre_vel);
                vel += factor * vol_top * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), centre_vel) - centre_vel);
                vel += factor * vol_bottom * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), centre_vel) - centre_vel);
                vel += factor * vol_front * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), centre_vel) - centre_vel);
                vel += factor * vol_back * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), centre_vel) - centre_vel);
                
                vel += factor * vol_top * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), 0.0) -
                                           get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(-1, 1, 0)), 0.0));
                vel += factor * vol_bottom * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), 0.0) -
                                              get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0));
                
                vel += factor * vol_front * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), 0.0) -
                                             get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 1)), 0.0));
                vel += factor * vol_back * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), 0.0) -
                                            get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0));
                
                node_vel_x[bucket_idx](i) = vel;
            }
            
            for(int i = 0; i < num_node_y; ++i) {
                const Vector3i& node_handle = node_indices_y[bucket_idx].segment<3>(i * 3);
                const scalar factor = coeff;
                
                const scalar vol_bottom = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
                const scalar vol_top = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
                
                const scalar vol_back = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, i, 4, 0, 0.0);
                const scalar vol_front = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, i, 4, 1, 0.0);
                
                const scalar vol_left = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, i, 4, 2, 0.0);
                const scalar vol_right = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, i, 4, 3, 0.0);
                
                const scalar centre_vel = node_vel_src_y[bucket_idx](i);
                scalar vel = centre_vel;
                
                vel += factor * vol_right * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), centre_vel) - centre_vel);
                vel += factor * vol_left * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), centre_vel) - centre_vel);
                vel += 2.0 * factor * vol_top * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), centre_vel) - centre_vel);
                vel += 2.0 * factor * vol_bottom * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), centre_vel) - centre_vel);
                vel += factor * vol_front * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), centre_vel) - centre_vel);
                vel += factor * vol_back * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), centre_vel) - centre_vel);
                
                vel += factor * vol_right * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), 0.0) -
                                             get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, -1, 0)), 0.0));
                vel += factor * vol_left * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), 0.0) -
                                            get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0));
                
                vel += factor * vol_front * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), 0.0) -
                                             get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 1)), 0.0));
                vel += factor * vol_back * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), 0.0) -
                                            get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0));
                
                node_vel_y[bucket_idx](i) = vel;
            }
            
            for(int i = 0; i < num_node_z; ++i) {
                const Vector3i& node_handle = node_indices_z[bucket_idx].segment<3>(i * 3);
                const scalar factor = coeff;
                
                const scalar vol_back = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
                const scalar vol_front = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
                
                const scalar vol_bottom = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, i, 4, 0, 0.0);
                const scalar vol_top = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, i, 4, 1, 0.0);
                
                const scalar vol_left = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, i, 4, 2, 0.0);
                const scalar vol_right = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, i, 4, 3, 0.0);
                
                const scalar centre_vel = node_vel_src_z[bucket_idx](i);
                scalar vel = centre_vel;
                
                vel += factor * vol_right * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), centre_vel) - centre_vel);
                vel += factor * vol_left * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), centre_vel) - centre_vel);
                vel += factor * vol_top * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), centre_vel) - centre_vel);
                vel += factor * vol_bottom * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), centre_vel) - centre_vel);
                vel += 2.0 * factor * vol_front * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), centre_vel) - centre_vel);
                vel += 2.0 * factor * vol_back * (get_value(node_cpidx_z, node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), centre_vel) - centre_vel);
                
                vel += factor * vol_right * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), 0.0) -
                                             get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, -1)), 0.0));
                vel += factor * vol_left * (get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), 0.0) -
                                            get_value(node_cpidx_x, node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0));
                
                vel += factor * vol_top * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), 0.0) -
                                           get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, -1)), 0.0));
                vel += factor * vol_bottom * (get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), 0.0) -
                                              get_value(node_cpidx_y, node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0));
                
                node_vel_z[bucket_idx](i) = vel;
            }
        });
    }
}

