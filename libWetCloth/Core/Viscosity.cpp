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
#include "pcgsolver/pcg_solver.h"

#include <iostream>
#include <cstdlib>

using namespace robertbridson;

namespace viscosity {
	template<typename T>
	T get_value
	(const std::vector< Eigen::Matrix<T, Eigen::Dynamic, 1> >& data,
	 const Vector3i& bucket_handle,
	 const Vector3i& node_handle,
	 const T& default_val,
	 const int num_nodes,
	 const Sorter& buckets)
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
		return data[bucket_idx][node_idx];
	};
	
	bool get_node_index(
						const std::vector< unsigned char >& bucket_activated,
						const Vector3i& bucket_handle,
						const Vector3i& node_handle,
						const int num_nodes,
						const Sorter& buckets,
						Vector2i& bucket_node)
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
			bucket_node(0) = bucket_node(1) = -1;
			return false;
		}
		
		const int bucket_idx = buckets.bucket_index(cur_bucket_handle);
		if(!bucket_activated[bucket_idx]) {
			bucket_node(0) = bucket_node(1) = -1;
			return false;
		}
		
		const int node_idx = cur_node_handle(2) * num_nodes * num_nodes + cur_node_handle(1) * num_nodes + cur_node_handle(0);

		bucket_node(0) = bucket_idx;
		bucket_node(1) = node_idx;
		return true;
	};
	
	void applyNodeViscosityImplicit( const TwoDScene& scene,
									const std::vector< VectorXi >& node_global_indices_x,
									const std::vector< VectorXi >& node_global_indices_y,
									const std::vector< VectorXi >& node_global_indices_z,
									int offset_nodes_x,
									int offset_nodes_y,
									int offset_nodes_z,
									const SparseMatrix< scalar >& matrix,
									const std::vector< scalar >& rhs,
									std::vector< scalar >& soln,
									std::vector< VectorXs >& node_vel_x,
									std::vector< VectorXs >& node_vel_y,
									std::vector< VectorXs >& node_vel_z,
									scalar& residual,
									int& iter_out,
									const scalar& criterion,
									int maxiters)
	{
		soln.assign(rhs.size(), 0.0);
		
		PCGSolver<double> solver;
		solver.set_solver_parameters(criterion, maxiters, 0.97, 0.1);
		bool success = false;
		
		success = solver.solve(matrix, rhs, soln, residual, iter_out);
		if(!success) {
			std::cerr << "\n\n\n**********VISCOSITY FAILED**************\n\n\n" << std::endl;
			exit(0);
		}
		
		const std::vector< VectorXuc >& node_state_u = scene.getNodeStateX();
		const std::vector< VectorXuc >& node_state_v = scene.getNodeStateY();
		const std::vector< VectorXuc >& node_state_w = scene.getNodeStateZ();
		
		const std::vector< VectorXs >& node_solid_u = scene.getNodeSolidVelX();
		const std::vector< VectorXs >& node_solid_v = scene.getNodeSolidVelY();
		const std::vector< VectorXs >& node_solid_w = scene.getNodeSolidVelZ();
		
		const Sorter& buckets = scene.getParticleBuckets();
		buckets.for_each_bucket([&] (int bucket_idx) {
			const int num_nodes_x = node_global_indices_x[bucket_idx].size();
			const int num_nodes_y = node_global_indices_y[bucket_idx].size();
			const int num_nodes_z = node_global_indices_z[bucket_idx].size();
			
			for(int i = 0; i < num_nodes_x; ++i)
			{
				if(node_state_u[bucket_idx][i] == NS_FLUID) {
					const int dof_idx = node_global_indices_x[bucket_idx][i];
					if(dof_idx < 0) continue;
					node_vel_x[bucket_idx][i] = soln[dof_idx + offset_nodes_x];
				} else {
					node_vel_x[bucket_idx][i] = node_solid_u[bucket_idx][i];
				}
			}
			
			for(int i = 0; i < num_nodes_y; ++i)
			{
				if(node_state_v[bucket_idx][i] == NS_FLUID) {
					const int dof_idx = node_global_indices_y[bucket_idx][i];
					if(dof_idx < 0) continue;
					node_vel_y[bucket_idx][i] = soln[dof_idx + offset_nodes_y];
				} else {
					node_vel_y[bucket_idx][i] = node_solid_v[bucket_idx][i];
				}
			}
			
			for(int i = 0; i < num_nodes_z; ++i)
			{
				if(node_state_w[bucket_idx][i] == NS_FLUID) {
					const int dof_idx = node_global_indices_z[bucket_idx][i];
					if(dof_idx < 0) continue;
					node_vel_z[bucket_idx][i] = soln[dof_idx + offset_nodes_z];
				} else {
					node_vel_z[bucket_idx][i] = node_solid_w[bucket_idx][i];
				}
			}
		});
	}
	
	void updateViscosityRHS( const TwoDScene& scene,
							const std::vector< VectorXi >& node_global_indices_x,
							const std::vector< VectorXi >& node_global_indices_y,
							const std::vector< VectorXi >& node_global_indices_z,
							const std::vector< Vector2i >& effective_node_indices_x,
							const std::vector< Vector2i >& effective_node_indices_y,
							const std::vector< Vector2i >& effective_node_indices_z,
							const std::vector< VectorXs >& node_vel_src_x,
							const std::vector< VectorXs >& node_vel_src_y,
							const std::vector< VectorXs >& node_vel_src_z,
							std::vector< scalar >& rhs,
							int offset_nodes_x,
							int offset_nodes_y,
							int offset_nodes_z,
							const scalar& dt  )
	{
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
		
		const std::vector< VectorXuc >& node_state_u = scene.getNodeStateX();
		const std::vector< VectorXuc >& node_state_v = scene.getNodeStateY();
		const std::vector< VectorXuc >& node_state_w = scene.getNodeStateZ();
		
		const int num_nodes = scene.getDefaultNumNodes();
		const Sorter& buckets = scene.getParticleBuckets();
		
		const scalar dx = scene.getCellSize();
		const scalar factor = dt / (dx * dx) * scene.getLiquidInfo().viscosity / scene.getLiquidInfo().liquid_density;
		
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
		
		const int total_num_nodes_x = (int) effective_node_indices_x.size();
		const int total_num_nodes_y = (int) effective_node_indices_y.size();
		const int total_num_nodes_z = (int) effective_node_indices_z.size();
		
        const std::vector<unsigned char>& bucket_activated = scene.getBucketActivated();
        
		threadutils::for_each(0, total_num_nodes_x, [&] (int dof_idx) {
			const Vector2i& dof_loc = effective_node_indices_x[dof_idx];
			const int bucket_idx = dof_loc[0];
			const int node_idx = dof_loc[1];
			
			const int index = dof_idx + offset_nodes_x;
			
			rhs[index] = node_liquid_u_vf[bucket_idx][node_idx] * node_vel_src_x[bucket_idx][node_idx];
			
			const scalar vol_left = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, node_idx, 2, 0, 0.0);
			const scalar vol_right = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, node_idx, 2, 1, 0.0);
			
			const scalar vol_back = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, node_idx, 4, 0, 0.0);
			const scalar vol_front = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, node_idx, 4, 1, 0.0);
			
			const scalar vol_bottom = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, node_idx, 4, 2, 0.0);
			const scalar vol_top = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, node_idx, 4, 3, 0.0);
			
			const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
            const Vector3i& node_handle = scene.getNodeHandle(node_idx);
			
			Vector2i cur_node;
			if(vol_right > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
			}
			
			if(vol_left > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
			}
			
			if(vol_top > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_top;
			}
			
			if(vol_bottom > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_bottom;
			}
			
			if(vol_front > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_front;
			}
			
			if(vol_back > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_back;
			}
			
			if(vol_top > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
			}
			
			if(vol_bottom > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
			}
			
			if(vol_front > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
			}
			
			if(vol_back > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
			}
		});
		
		threadutils::for_each(0, total_num_nodes_y, [&] (int dof_idx) {
			const Vector2i& dof_loc = effective_node_indices_y[dof_idx];
			const int bucket_idx = dof_loc[0];
			const int node_idx = dof_loc[1];
			
			const int index = dof_idx + offset_nodes_y;
			
			rhs[index] = node_liquid_v_vf[bucket_idx][node_idx] * node_vel_src_y[bucket_idx][node_idx];
			
			const scalar vol_bottom = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, node_idx, 2, 0, 0.0);
			const scalar vol_top = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, node_idx, 2, 1, 0.0);
			
			const scalar vol_back = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, node_idx, 4, 0, 0.0);
			const scalar vol_front = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, node_idx, 4, 1, 0.0);
			
			const scalar vol_left = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, node_idx, 4, 2, 0.0);
			const scalar vol_right = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, node_idx, 4, 3, 0.0);
			
			const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
            const Vector3i& node_handle = scene.getNodeHandle(node_idx);
			
			Vector2i cur_node;
			if(vol_right > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_right;
			}
			
			if(vol_left > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_left;
			}
			
			if(vol_top > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
			}
			
			if(vol_bottom > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
			}
			
			if(vol_front > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_front;
			}
			
			if(vol_back > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_back;
			}
			
			if(vol_right > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, -1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
			}
			
			if(vol_left > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
			}
			
			if(vol_front > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
			}
			
			if(vol_back > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
			}
		});
		
		threadutils::for_each(0, total_num_nodes_z, [&] (int dof_idx) {
			const Vector2i& dof_loc = effective_node_indices_z[dof_idx];
			const int bucket_idx = dof_loc[0];
			const int node_idx = dof_loc[1];
			
			const int index = dof_idx + offset_nodes_z;
			
			rhs[index] = node_liquid_w_vf[bucket_idx][node_idx] * node_vel_src_z[bucket_idx][node_idx];
			
			const scalar vol_back = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, node_idx, 2, 0, 0.0);
			const scalar vol_front = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, node_idx, 2, 1, 0.0);
			
			const scalar vol_bottom = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, node_idx, 4, 0, 0.0);
			const scalar vol_top = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, node_idx, 4, 1, 0.0);
			
			const scalar vol_left = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, node_idx, 4, 2, 0.0);
			const scalar vol_right = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, node_idx, 4, 3, 0.0);
			
			const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
            const Vector3i& node_handle = scene.getNodeHandle(node_idx);
			
			Vector2i cur_node;
			if(vol_right > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_right;
			}
			
			if(vol_left > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_left;
			}
			
			if(vol_top > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_top;
			}
			
			if(vol_bottom > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_bottom;
			}
			
			if(vol_front > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
			}
			
			if(vol_back > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
			}
			
			if(vol_right > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
			}
			
			if(vol_left > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
			}
			
			if(vol_top > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
			}
			
			if(vol_bottom > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
			}
		});
	}
	
	void constructViscosityMatrixRHS( const TwoDScene& scene,
									 std::vector< VectorXi >& node_global_indices_x,
									 std::vector< VectorXi >& node_global_indices_y,
									 std::vector< VectorXi >& node_global_indices_z,
									 std::vector< Vector2i >& effective_node_indices_x,
									 std::vector< Vector2i >& effective_node_indices_y,
									 std::vector< Vector2i >& effective_node_indices_z,
									 const std::vector< VectorXs >& node_vel_src_x,
									 const std::vector< VectorXs >& node_vel_src_y,
									 const std::vector< VectorXs >& node_vel_src_z,
									 SparseMatrix< scalar >& matrix,
									 std::vector< scalar >& rhs,
									 int& offset_nodes_x,
									 int& offset_nodes_y,
									 int& offset_nodes_z,
									 const scalar& dt )
	{
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
		
		const std::vector< VectorXuc >& node_state_u = scene.getNodeStateX();
		const std::vector< VectorXuc >& node_state_v = scene.getNodeStateY();
		const std::vector< VectorXuc >& node_state_w = scene.getNodeStateZ();
		
		const int num_nodes = scene.getDefaultNumNodes();
		const Sorter& buckets = scene.getParticleBuckets();
		
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
		const scalar factor = dt / (dx * dx) * scene.getLiquidInfo().viscosity / scene.getLiquidInfo().liquid_density;
		
		std::vector<int> num_effective_nodes_x( buckets.size() );
		std::vector<int> num_effective_nodes_y( buckets.size() );
		std::vector<int> num_effective_nodes_z( buckets.size() );
		
		// assign global indices to nodes
		buckets.for_each_bucket([&] (int bucket_idx) {
			const int num_nodes_x = (int) node_vel_src_x[bucket_idx].size();
			const int num_nodes_y = (int) node_vel_src_y[bucket_idx].size();
			const int num_nodes_z = (int) node_vel_src_z[bucket_idx].size();
			
			// filter x
			int k = 0;
			for(int i = 0; i < num_nodes_x; ++i)
			{
				node_global_indices_x[bucket_idx][i] = -1;
				
				if(node_state_u[bucket_idx][i] == NS_FLUID) {
					const scalar vol_u = node_liquid_u_vf[bucket_idx][i];
					if(vol_u > 0.0) {
						node_global_indices_x[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_left = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
					if(vol_left > 0.0) {
						node_global_indices_x[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_right = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
					if(vol_right > 0.0) {
						node_global_indices_x[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_back = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, i, 4, 0, 0.0);
					if(vol_back > 0.0) {
						node_global_indices_x[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_front = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, i, 4, 1, 0.0);
					if(vol_front > 0.0) {
						node_global_indices_x[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_bottom = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, i, 4, 2, 0.0);
					if(vol_bottom > 0.0) {
						node_global_indices_x[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_top = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, i, 4, 3, 0.0);
					if(vol_top > 0.0) {
						node_global_indices_x[bucket_idx][i] = k++;
						continue;
					}
				}
			}
			num_effective_nodes_x[bucket_idx] = k;
			
			// filter y
			k = 0;
			for(int i = 0; i < num_nodes_y; ++i)
			{
				node_global_indices_y[bucket_idx][i] = -1;
				
				if(node_state_v[bucket_idx][i] == NS_FLUID) {
					const scalar vol_v = node_liquid_v_vf[bucket_idx][i];
					if(vol_v > 0.0) {
						node_global_indices_y[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_bottom = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
					if(vol_bottom > 0.0) {
						node_global_indices_y[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_top = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
					if(vol_top > 0.0) {
						node_global_indices_y[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_back = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, i, 4, 0, 0.0);
					if(vol_back > 0.0) {
						node_global_indices_y[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_front = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, i, 4, 1, 0.0);
					if(vol_front > 0.0) {
						node_global_indices_y[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_left = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, i, 4, 2, 0.0);
					if(vol_left > 0.0) {
						node_global_indices_y[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_right = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, i, 4, 3, 0.0);
					if(vol_right > 0.0) {
						node_global_indices_y[bucket_idx][i] = k++;
						continue;
					}
				}
			}
			num_effective_nodes_y[bucket_idx] = k;
			
			k = 0;
			for(int i = 0; i < num_nodes_z; ++i)
			{
				node_global_indices_z[bucket_idx][i] = -1;
				
				if(node_state_w[bucket_idx][i] == NS_FLUID) {
					const scalar vol_w = node_liquid_w_vf[bucket_idx][i];
					if(vol_w > 0.0) {
						node_global_indices_z[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_back = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
					if(vol_back > 0.0) {
						node_global_indices_z[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_front = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
					if(vol_front > 0.0) {
						node_global_indices_z[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_bottom = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, i, 4, 0, 0.0);
					if(vol_bottom > 0.0) {
						node_global_indices_z[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_top = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, i, 4, 1, 0.0);
					if(vol_top > 0.0) {
						node_global_indices_z[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_left = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, i, 4, 2, 0.0);
					if(vol_left > 0.0) {
						node_global_indices_z[bucket_idx][i] = k++;
						continue;
					}
					const scalar vol_right = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, i, 4, 3, 0.0);
					if(vol_right > 0.0) {
						node_global_indices_z[bucket_idx][i] = k++;
						continue;
					}
				}
			}
			num_effective_nodes_z[bucket_idx] = k;
		});
		
		std::partial_sum(num_effective_nodes_x.begin(), num_effective_nodes_x.end(), num_effective_nodes_x.begin());
		std::partial_sum(num_effective_nodes_y.begin(), num_effective_nodes_y.end(), num_effective_nodes_y.begin());
		std::partial_sum(num_effective_nodes_z.begin(), num_effective_nodes_z.end(), num_effective_nodes_z.begin());
		
		const int total_num_nodes_x = num_effective_nodes_x[num_effective_nodes_x.size() - 1];
		const int total_num_nodes_y = num_effective_nodes_y[num_effective_nodes_y.size() - 1];
		const int total_num_nodes_z = num_effective_nodes_z[num_effective_nodes_z.size() - 1];
		
		if(total_num_nodes_x == 0 && total_num_nodes_y == 0 && total_num_nodes_z == 0) return;
		
		offset_nodes_x = 0;
		offset_nodes_y = total_num_nodes_x;
		offset_nodes_z = total_num_nodes_x + total_num_nodes_y;
		
		const int total_num_nodes = total_num_nodes_x + total_num_nodes_y + total_num_nodes_z;
		
		effective_node_indices_x.resize(total_num_nodes_x);
		effective_node_indices_y.resize(total_num_nodes_y);
		effective_node_indices_z.resize(total_num_nodes_z);
		
		if((int) rhs.size() != total_num_nodes) {
			rhs.resize(total_num_nodes);
			matrix.resize(total_num_nodes);
		}
		
		matrix.zero();
		
		// assign indices from reduced to global
		buckets.for_each_bucket([&] (int bucket_idx) {
			const int num_nodes_x = (int) node_vel_src_x[bucket_idx].size();
			const int num_nodes_y = (int) node_vel_src_y[bucket_idx].size();
			const int num_nodes_z = (int) node_vel_src_z[bucket_idx].size();
			
			for(int i = 0; i < num_nodes_x; ++i)
			{
				if(node_global_indices_x[bucket_idx][i] >= 0) {
					if(bucket_idx > 0) {
						node_global_indices_x[bucket_idx][i] += num_effective_nodes_x[bucket_idx - 1];
					}
					
					const int global_idx = node_global_indices_x[bucket_idx][i];
					effective_node_indices_x[global_idx] = Vector2i(bucket_idx, i);
				}
			}
			
			for(int i = 0; i < num_nodes_y; ++i)
			{
				if(node_global_indices_y[bucket_idx][i] >= 0) {
					if(bucket_idx > 0) {
						node_global_indices_y[bucket_idx][i] += num_effective_nodes_y[bucket_idx - 1];
					}
					
					const int global_idx = node_global_indices_y[bucket_idx][i];
					effective_node_indices_y[global_idx] = Vector2i(bucket_idx, i);
				}
			}
			
			for(int i = 0; i < num_nodes_z; ++i)
			{
				if(node_global_indices_z[bucket_idx][i] >= 0) {
					if(bucket_idx > 0) {
						node_global_indices_z[bucket_idx][i] += num_effective_nodes_z[bucket_idx - 1];
					}
					
					const int global_idx = node_global_indices_z[bucket_idx][i];
					effective_node_indices_z[global_idx] = Vector2i(bucket_idx, i);
				}
			}
		});
		
		auto u_ind = [&] (const Vector2i& node) -> int {
			assert(node(0) >= 0 && node(0) < (int)node_global_indices_x.size());
			assert(node(1) >= 0 && node(1) < (int)node_global_indices_x[node(0)].size());
			assert(node_global_indices_x[node(0)][node(1)] >= 0);
			
			return node_global_indices_x[node(0)][node(1)] + offset_nodes_x;
		};
		
		auto v_ind = [&] (const Vector2i& node) -> int {
			assert(node(0) >= 0 && node(0) < (int)node_global_indices_y.size());
			assert(node(1) >= 0 && node(1) < (int)node_global_indices_y[node(0)].size());
			assert(node_global_indices_y[node(0)][node(1)] >= 0);
			
			return node_global_indices_y[node(0)][node(1)] + offset_nodes_y;
		};
		
		auto w_ind = [&] (const Vector2i& node) -> int {
			assert(node(0) >= 0 && node(0) < (int)node_global_indices_z.size());
			assert(node(1) >= 0 && node(1) < (int)node_global_indices_z[node(0)].size());
			assert(node_global_indices_z[node(0)][node(1)] >= 0);
			
			return node_global_indices_z[node(0)][node(1)] + offset_nodes_z;
		};
		
        const std::vector<unsigned char>& bucket_activated = scene.getBucketActivated();
        
		// assamble X-matrix
		threadutils::for_each(0, total_num_nodes_x, [&] (int dof_idx) {
			const Vector2i& dof_loc = effective_node_indices_x[dof_idx];
			const int bucket_idx = dof_loc[0];
			const int node_idx = dof_loc[1];
			
			const int index = dof_idx + offset_nodes_x;
			
			rhs[index] = node_liquid_u_vf[bucket_idx][node_idx] * node_vel_src_x[bucket_idx][node_idx];
			matrix.set_element(index, index, node_liquid_u_vf[bucket_idx][node_idx]);
			
			const scalar vol_left = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, node_idx, 2, 0, 0.0);
			const scalar vol_right = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, node_idx, 2, 1, 0.0);
			
			const scalar vol_back = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, node_idx, 4, 0, 0.0);
			const scalar vol_front = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, node_idx, 4, 1, 0.0);
			
			const scalar vol_bottom = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, node_idx, 4, 2, 0.0);
			const scalar vol_top = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, node_idx, 4, 3, 0.0);
			
			const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
            const Vector3i& node_handle = scene.getNodeHandle(node_idx);
			
			Vector2i cur_node;
			if(vol_right > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, 2 * factor * vol_right);
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, u_ind(cur_node), -2 * factor * vol_right);
				else if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
			}
			
			if(vol_left > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, 2 * factor * vol_left);
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, u_ind(cur_node), -2 * factor * vol_left);
				else if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
			}
			
			if(vol_top > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_top);
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, u_ind(cur_node), -factor * vol_top);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_top;
			}
			
			if(vol_bottom > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_bottom);
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, u_ind(cur_node), -factor * vol_bottom);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_bottom;
			}
			
			if(vol_front > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_front);
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, u_ind(cur_node), -factor * vol_front);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_front;
			}
			
			if(vol_back > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_back);
				const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, u_ind(cur_node), -factor * vol_back);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_back;
			}
			
			if(vol_top > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), -factor * vol_top);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), factor * vol_top);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
			}
			
			if(vol_bottom > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), factor * vol_bottom);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), -factor * vol_bottom);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
			}
			
			if(vol_front > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), -factor * vol_front);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), factor * vol_front);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
			}
			
			if(vol_back > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), factor * vol_back);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), -factor * vol_back);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
			}
		});
		
		threadutils::for_each(0, total_num_nodes_y, [&] (int dof_idx) {
			const Vector2i& dof_loc = effective_node_indices_y[dof_idx];
			const int bucket_idx = dof_loc[0];
			const int node_idx = dof_loc[1];
			
			const int index = dof_idx + offset_nodes_y;
			
			rhs[index] = node_liquid_v_vf[bucket_idx][node_idx] * node_vel_src_y[bucket_idx][node_idx];
			matrix.set_element(index, index, node_liquid_v_vf[bucket_idx][node_idx]);
			
			const scalar vol_bottom = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, node_idx, 2, 0, 0.0);
			const scalar vol_top = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, node_idx, 2, 1, 0.0);
			
			const scalar vol_back = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, node_idx, 4, 0, 0.0);
			const scalar vol_front = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, node_idx, 4, 1, 0.0);
			
			const scalar vol_left = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, node_idx, 4, 2, 0.0);
			const scalar vol_right = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, node_idx, 4, 3, 0.0);
			
			const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
            const Vector3i& node_handle = scene.getNodeHandle(node_idx);
			
			Vector2i cur_node;
			if(vol_right > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_right);
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, v_ind(cur_node), -factor * vol_right);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_right;
			}
			
			if(vol_left > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_left);
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, v_ind(cur_node), -factor * vol_left);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_left;
			}
			
			if(vol_top > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, 2 * factor * vol_top);
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, v_ind(cur_node), -2. * factor * vol_top);
				else if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
			}
			
			if(vol_bottom > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, 2 * factor * vol_bottom);
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, v_ind(cur_node), -2. * factor * vol_bottom);
				else if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
			}
			
			if(vol_front > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_front);
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, v_ind(cur_node), -factor * vol_front);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_front;
			}
			
			if(vol_back > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_back);
				const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, v_ind(cur_node), -factor * vol_back);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_back;
			}
			
			if(vol_right > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), -factor * vol_right);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, -1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), factor * vol_right);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
			}
			
			if(vol_left > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), factor * vol_left);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), -factor * vol_left);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
			}
			
			if(vol_front > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), -factor * vol_front);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), factor * vol_front);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
				}
			}
			
			if(vol_back > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), factor * vol_back);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, w_ind(cur_node), -factor * vol_back);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
				}
			}
		});
		
		threadutils::for_each(0, total_num_nodes_z, [&] (int dof_idx) {
			const Vector2i& dof_loc = effective_node_indices_z[dof_idx];
			const int bucket_idx = dof_loc[0];
			const int node_idx = dof_loc[1];
			
			const int index = dof_idx + offset_nodes_z;
			
			rhs[index] = node_liquid_w_vf[bucket_idx][node_idx] * node_vel_src_z[bucket_idx][node_idx];
			matrix.set_element(index, index, node_liquid_w_vf[bucket_idx][node_idx]);
			
			const scalar vol_back = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, node_idx, 2, 0, 0.0);
			const scalar vol_front = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, node_idx, 2, 1, 0.0);
			
			const scalar vol_bottom = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, node_idx, 4, 0, 0.0);
			const scalar vol_top = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, node_idx, 4, 1, 0.0);
			
			const scalar vol_left = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, node_idx, 4, 2, 0.0);
			const scalar vol_right = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, node_idx, 4, 3, 0.0);
			
			const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
            const Vector3i& node_handle = scene.getNodeHandle(node_idx);
			
			Vector2i cur_node;
			if(vol_right > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_right);
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, w_ind(cur_node), -factor * vol_right);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_right;
			}
			
			if(vol_left > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_left);
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, w_ind(cur_node), -factor * vol_left);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_left;
			}
			
			if(vol_top > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_top);
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, w_ind(cur_node), -factor * vol_top);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_top;
			}
			
			if(vol_bottom > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, factor * vol_bottom);
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, w_ind(cur_node), -factor * vol_bottom);
				else if(state == NS_SOLID)
					rhs[index] -= -node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_bottom;
			}
			
			if(vol_front > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, 2. * factor * vol_front);
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, w_ind(cur_node), -2. * factor * vol_front);
				else if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_front;
			}
			
			if(vol_back > 0. && get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle - Vector3i(0, 0, 1)), num_nodes, buckets, cur_node))
			{
				matrix.add_to_element(index, index, 2. * factor * vol_back);
				const NODE_STATE state = (NODE_STATE) node_state_w[cur_node(0)][cur_node(1)];
				if(state == NS_FLUID)
					matrix.add_to_element(index, w_ind(cur_node), -2. * factor * vol_back);
				else if(state == NS_SOLID)
					rhs[index] -= -2. * node_vel_src_z[cur_node(0)][cur_node(1)] * factor * vol_back;
			}
			
			if(vol_right > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), -factor * vol_right);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), factor * vol_right);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_right;
				}
			}
			
			if(vol_left > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), factor * vol_left);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_u[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, u_ind(cur_node), -factor * vol_left);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_x[cur_node(0)][cur_node(1)] * factor * vol_left;
				}
			}
			
			if(vol_top > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), -factor * vol_top);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), factor * vol_top);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_top;
				}
			}
			
			if(vol_bottom > 0.) {
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), factor * vol_bottom);
					else if(state == NS_SOLID)
						rhs[index] -= node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
				
				if(get_node_index(bucket_activated, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), num_nodes, buckets, cur_node))
				{
					const NODE_STATE state = (NODE_STATE) node_state_v[cur_node(0)][cur_node(1)];
					if(state == NS_FLUID)
						matrix.add_to_element(index, v_ind(cur_node), -factor * vol_bottom);
					else if(state == NS_SOLID)
						rhs[index] -= -node_vel_src_y[cur_node(0)][cur_node(1)] * factor * vol_bottom;
				}
			}
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
		const std::vector< VectorXi >& node_index_ex = scene.getNodeIndexEdgeX();
		const std::vector< VectorXi >& node_index_ey = scene.getNodeIndexEdgeY();
		const std::vector< VectorXi >& node_index_ez = scene.getNodeIndexEdgeZ();
		
		const std::vector< VectorXi >& node_index_p_x = scene.getNodePressureIndexX();
		const std::vector< VectorXi >& node_index_p_y = scene.getNodePressureIndexY();
		const std::vector< VectorXi >& node_index_p_z = scene.getNodePressureIndexZ();
		
		const std::vector< VectorXs >& node_liquid_c_vf = scene.getNodeLiquidVolFracCentre();
		const std::vector< VectorXs >& node_liquid_ex_vf = scene.getNodeLiquidVolFracEX();
		const std::vector< VectorXs >& node_liquid_ey_vf = scene.getNodeLiquidVolFracEY();
		const std::vector< VectorXs >& node_liquid_ez_vf = scene.getNodeLiquidVolFracEZ();
		
		const int num_nodes = scene.getDefaultNumNodes();
		const Sorter& buckets = scene.getParticleBuckets();
		
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
            const int num_node = scene.getNumNodes(bucket_idx);
            
            if(!scene.isBucketActivated(bucket_idx)) return;
			
			const Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);
			
			for(int i = 0; i < num_node; ++i) {
                const Vector3i& node_handle = scene.getNodeHandle(i);
				const scalar factor = coeff;
				
				const scalar vol_left = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
				const scalar vol_right = get_value_fast(node_index_p_x, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
				
				const scalar vol_back = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, i, 4, 0, 0.0);
				const scalar vol_front = get_value_fast(node_index_ex, node_liquid_ey_vf, bucket_idx, i, 4, 1, 0.0);
				
				const scalar vol_bottom = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, i, 4, 2, 0.0);
				const scalar vol_top = get_value_fast(node_index_ex, node_liquid_ez_vf, bucket_idx, i, 4, 3, 0.0);
				
				const scalar centre_vel = node_vel_src_x[bucket_idx](i);
				scalar vel = centre_vel;
				
				vel += 2.0 * factor * vol_right * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += 2.0 * factor * vol_left * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_top * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_bottom * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_front * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_back * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), centre_vel, num_nodes, buckets) - centre_vel);
				
				vel += factor * vol_top * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), 0.0, num_nodes, buckets) -
										   get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(-1, 1, 0)), 0.0, num_nodes, buckets));
				vel += factor * vol_bottom * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), 0.0, num_nodes, buckets) -
											  get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0, num_nodes, buckets));
				
				vel += factor * vol_front * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), 0.0, num_nodes, buckets) -
											 get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 1)), 0.0, num_nodes, buckets));
				vel += factor * vol_back * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), 0.0, num_nodes, buckets) -
											get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0, num_nodes, buckets));
				
				node_vel_x[bucket_idx](i) = vel;
			}
			
			for(int i = 0; i < num_node; ++i) {
				const Vector3i& node_handle = scene.getNodeHandle(i);
				const scalar factor = coeff;
				
				const scalar vol_bottom = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
				const scalar vol_top = get_value_fast(node_index_p_y, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
				
				const scalar vol_back = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, i, 4, 0, 0.0);
				const scalar vol_front = get_value_fast(node_index_ey, node_liquid_ex_vf, bucket_idx, i, 4, 1, 0.0);
				
				const scalar vol_left = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, i, 4, 2, 0.0);
				const scalar vol_right = get_value_fast(node_index_ey, node_liquid_ez_vf, bucket_idx, i, 4, 3, 0.0);
				
				const scalar centre_vel = node_vel_src_y[bucket_idx](i);
				scalar vel = centre_vel;
				
				vel += factor * vol_right * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_left * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += 2.0 * factor * vol_top * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += 2.0 * factor * vol_bottom * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_front * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_back * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), centre_vel, num_nodes, buckets) - centre_vel);
				
				vel += factor * vol_right * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), 0.0, num_nodes, buckets) -
											 get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, -1, 0)), 0.0, num_nodes, buckets));
				vel += factor * vol_left * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), 0.0, num_nodes, buckets) -
											get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0, num_nodes, buckets));
				
				vel += factor * vol_front * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), 0.0, num_nodes, buckets) -
											 get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 1)), 0.0, num_nodes, buckets));
				vel += factor * vol_back * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), 0.0, num_nodes, buckets) -
											get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0, num_nodes, buckets));
				
				node_vel_y[bucket_idx](i) = vel;
			}
			
			for(int i = 0; i < num_node; ++i) {
				const Vector3i& node_handle = scene.getNodeHandle(i);
				const scalar factor = coeff;
				
				const scalar vol_back = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, i, 2, 0, 0.0);
				const scalar vol_front = get_value_fast(node_index_p_z, node_liquid_c_vf, bucket_idx, i, 2, 1, 0.0);
				
				const scalar vol_bottom = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, i, 4, 0, 0.0);
				const scalar vol_top = get_value_fast(node_index_ez, node_liquid_ex_vf, bucket_idx, i, 4, 1, 0.0);
				
				const scalar vol_left = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, i, 4, 2, 0.0);
				const scalar vol_right = get_value_fast(node_index_ez, node_liquid_ey_vf, bucket_idx, i, 4, 3, 0.0);
				
				const scalar centre_vel = node_vel_src_z[bucket_idx](i);
				scalar vel = centre_vel;
				
				vel += factor * vol_right * (get_value( node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_left * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(-1, 0, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_top * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += factor * vol_bottom * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, -1, 0)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += 2.0 * factor * vol_front * (get_value(node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 1)), centre_vel, num_nodes, buckets) - centre_vel);
				vel += 2.0 * factor * vol_back * (get_value( node_vel_src_z, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), centre_vel, num_nodes, buckets) - centre_vel);
				
				vel += factor * vol_right * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, 0)), 0.0, num_nodes, buckets) -
											 get_value( node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(1, 0, -1)), 0.0, num_nodes, buckets));
				vel += factor * vol_left * (get_value(node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), 0.0, num_nodes, buckets) -
											get_value( node_vel_src_x, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0, num_nodes, buckets));
				
				vel += factor * vol_top * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, 0)), 0.0, num_nodes, buckets) -
										   get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 1, -1)), 0.0, num_nodes, buckets));
				vel += factor * vol_bottom * (get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, -1)), 0.0, num_nodes, buckets) -
											  get_value(node_vel_src_y, bucket_handle, Vector3i(node_handle + Vector3i(0, 0, 0)), 0.0, num_nodes, buckets));
				
				node_vel_z[bucket_idx](i) = vel;
			}
		});
	}
}

