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


#include "SceneStepper.h"
#include <numeric>

SceneStepper::~SceneStepper()
{}

bool SceneStepper::advectScene( TwoDScene& scene, scalar dt )
{
	VectorXs& x = scene.getX();
	const VectorXs& v = scene.getV();
	const VectorXs& fv = scene.getFluidV();
	
	assert(!std::isnan(v.sum()));
	assert(!std::isnan(x.sum()));
	
	const int num_parts = scene.getNumParticles();
	const int num_elasto_parts = scene.getNumElastoParticles();
	threadutils::for_each(0, num_parts, [&] (int pidx) {
		bool is_fluid = pidx >= num_elasto_parts;
		
		if(is_fluid) {
			x.segment<3>(pidx * 4) += fv.segment<3>(pidx * 4) * dt;
		} else {
            x.segment<4>(pidx * 4) += v.segment<4>(pidx * 4) * dt;
        }
	});

//	std::cout << "v: \n" << v << std::endl;
//	std::cout << "x: \n" << x << std::endl;

	return true;
}

void SceneStepper::setUseApic(bool apic)
{
	m_apic = apic;
}

bool SceneStepper::useApic() const
{
	return m_apic;
}

void SceneStepper::mapNodeToSoftParticles( const TwoDScene& scene, const std::vector< VectorXs >& node_vec_x, const std::vector< VectorXs >& node_vec_y, const std::vector< VectorXs >& node_vec_z, VectorXs& part_vec ) const
{
	const Sorter& buckets = scene.getParticleBuckets();
	
	part_vec.setZero();
	
    const scalar invD = scene.getInverseDCoeff();
	
	const int kernel_order = scene.getKernelOrder();
	
    const int num_soft_elasto = scene.getNumSoftElastoParticles();
    
    threadutils::for_each(0, num_soft_elasto, [&] (int pidx) {
		auto& indices_x = scene.getParticleNodesX(pidx);
		auto& indices_y = scene.getParticleNodesY(pidx);
		auto& indices_z = scene.getParticleNodesZ(pidx);
		
		auto& grads_x = scene.getParticleGradsX(pidx);
		auto& grads_y = scene.getParticleGradsY(pidx);
		auto& grads_z = scene.getParticleGradsZ(pidx);
		
		auto& weights = scene.getParticleWeights(pidx);
		const Vector3s& part_pos = scene.getX().segment<3>(pidx * 4);
		
        scalar sum_x(0.), sum_y(0.), sum_z(0.);
		for(int i = 0; i < 27; ++i) {
			if(!indices_x(i, 2) || weights(i, 0) == 0.0) continue;
			sum_x += node_vec_x[ indices_x(i, 0) ]( indices_x(i, 1) ) * weights(i, 0);
		}
        part_vec(pidx * 4 + 0) = sum_x;
        
		for(int i = 0; i < 27; ++i) {
			if(!indices_y(i, 2) || weights(i, 1) == 0.0) continue;
			sum_y += node_vec_y[ indices_y(i, 0) ]( indices_y(i, 1) ) * weights(i, 1);
		}
        part_vec(pidx * 4 + 1) = sum_y;
		
		for(int i = 0; i < 27; ++i) {
			if(!indices_z(i, 2) || weights(i, 2) == 0.0) continue;
			sum_z += node_vec_z[ indices_z(i, 0) ]( indices_z(i, 1) ) * weights(i, 2);
		}
        part_vec(pidx * 4 + 2) = sum_z;
	});
}

void SceneStepper::allocateLagrangianVectors( const TwoDScene& scene, VectorXs& vec )
{
    const int num_elasto = scene.getNumSoftElastoParticles();
    vec.resize(num_elasto * 4);
    vec.setZero();
}

void SceneStepper::buildLocalGlobalMapping( const TwoDScene& scene,
                                           std::vector< VectorXi >& node_global_indices_x,
                                           std::vector< VectorXi >& node_global_indices_y,
                                           std::vector< VectorXi >& node_global_indices_z,
                                           std::vector< Vector3i >& effective_node_indices,
                                           std::vector< Vector3i >& dof_ijk)
{
    const Sorter& buckets = scene.getParticleBuckets();
    const int bucket_num_cell = scene.getDefaultNumNodes();
    const int ni = buckets.ni * bucket_num_cell;
    const int nj = buckets.nj * bucket_num_cell;
    const int nk = buckets.nk * bucket_num_cell;
    
    std::vector<int> num_effective_nodes( buckets.size() );
    
    const std::vector< VectorXs >& node_elasto_vol_x = scene.getNodeVolX();
    const std::vector< VectorXs >& node_elasto_vol_y = scene.getNodeVolY();
    const std::vector< VectorXs >& node_elasto_vol_z = scene.getNodeVolZ();
    
    allocateNodeVectors(scene, node_global_indices_x, node_global_indices_y, node_global_indices_z);
    
    buckets.for_each_bucket([&] (int bucket_idx) {
        int k = 0;
        
        const int num_nodes_x = node_elasto_vol_x[bucket_idx].size();
        const int num_nodes_y = node_elasto_vol_y[bucket_idx].size();
        const int num_nodes_z = node_elasto_vol_z[bucket_idx].size();
        
        for(int i = 0; i < num_nodes_x; ++i)
        {
            node_global_indices_x[bucket_idx][i] = -1;
            
            if(node_elasto_vol_x[bucket_idx][i] > 0.0) {
                node_global_indices_x[bucket_idx][i] = k++;
            }
        }
        
        for(int i = 0; i < num_nodes_y; ++i)
        {
            node_global_indices_y[bucket_idx][i] = -1;
            
            if(node_elasto_vol_y[bucket_idx][i] > 0.0) {
                node_global_indices_y[bucket_idx][i] = k++;
            }
        }
        
        for(int i = 0; i < num_nodes_z; ++i)
        {
            node_global_indices_z[bucket_idx][i] = -1;
            
            if(node_elasto_vol_z[bucket_idx][i] > 0.0) {
                node_global_indices_z[bucket_idx][i] = k++;
            }
        }
        
        num_effective_nodes[bucket_idx] = k;
    });
    
    std::partial_sum(num_effective_nodes.begin(), num_effective_nodes.end(), num_effective_nodes.begin());
    
    const int total_num_nodes = num_effective_nodes[num_effective_nodes.size() - 1];
    if(total_num_nodes == 0) return;
    
    effective_node_indices.resize(total_num_nodes);
    dof_ijk.resize(total_num_nodes);
    
    const std::vector< VectorXi >& node_indices_x = scene.getNodeIndicesX();
    const std::vector< VectorXi >& node_indices_y = scene.getNodeIndicesY();
    const std::vector< VectorXi >& node_indices_z = scene.getNodeIndicesZ();
    
    buckets.for_each_bucket([&] (int bucket_idx) {
        const Vector3i handle = buckets.bucket_handle(bucket_idx);
        
        const int num_nodes_x = node_elasto_vol_x[bucket_idx].size();
        const int num_nodes_y = node_elasto_vol_y[bucket_idx].size();
        const int num_nodes_z = node_elasto_vol_z[bucket_idx].size();
        
        for(int i = 0; i < num_nodes_x; ++i)
        {
            if(node_global_indices_x[bucket_idx][i] >= 0) {
                if(bucket_idx > 0) {
                    node_global_indices_x[bucket_idx][i] += num_effective_nodes[bucket_idx - 1];
                }
                
                const int global_idx = node_global_indices_x[bucket_idx][i];
                
                effective_node_indices[global_idx] = Vector3i(bucket_idx, 0, i);
                
                const Vector3i& local_handle = node_indices_x[bucket_idx].segment<3>(i * 3);
                
                dof_ijk[global_idx] = Vector3i(handle(0) * bucket_num_cell * 3 + local_handle(0),
                                               handle(1) * bucket_num_cell + local_handle(1),
                                               handle(2) * bucket_num_cell + local_handle(2)
                                               );
            }
        }
        
        for(int i = 0; i < num_nodes_y; ++i)
        {
            if(node_global_indices_y[bucket_idx][i] >= 0) {
                if(bucket_idx > 0) {
                    node_global_indices_y[bucket_idx][i] += num_effective_nodes[bucket_idx - 1];
                }
                
                const int global_idx = node_global_indices_y[bucket_idx][i];
                
                effective_node_indices[global_idx] = Vector3i(bucket_idx, 1, i);
                
                const Vector3i& local_handle = node_indices_y[bucket_idx].segment<3>(i * 3);
                
                dof_ijk[global_idx] = Vector3i(handle(0) * bucket_num_cell * 3 + bucket_num_cell + local_handle(0),
                                               handle(1) * bucket_num_cell + local_handle(1),
                                               handle(2) * bucket_num_cell + local_handle(2)
                                               );
            }
        }
        
        for(int i = 0; i < num_nodes_z; ++i)
        {
            if(node_global_indices_z[bucket_idx][i] >= 0) {
                if(bucket_idx > 0) {
                    node_global_indices_z[bucket_idx][i] += num_effective_nodes[bucket_idx - 1];
                }
                
                const int global_idx = node_global_indices_z[bucket_idx][i];
                
                effective_node_indices[global_idx] = Vector3i(bucket_idx, 2, i);
                
                const Vector3i& local_handle = node_indices_z[bucket_idx].segment<3>(i * 3);
                
                dof_ijk[global_idx] = Vector3i(handle(0) * bucket_num_cell * 3 + bucket_num_cell * 2 + local_handle(0),
                                               handle(1) * bucket_num_cell + local_handle(1),
                                               handle(2) * bucket_num_cell + local_handle(2)
                                               );
            }
        }
    });
}

void SceneStepper::buildNodeToSoftParticlesMat( const TwoDScene& scene,
                                               const std::vector< VectorXi >& node_global_indices_x,
                                               const std::vector< VectorXi >& node_global_indices_y,
                                               const std::vector< VectorXi >& node_global_indices_z,
                                               const std::vector< Vector3i >& effective_node_indices,
                                               TripletXs& tri_W,
                                               SparseXs& W ) const
{
    const Sorter& buckets = scene.getParticleBuckets();

    const int num_soft_elasto = scene.getNumSoftElastoParticles();

    const int system_size = effective_node_indices.size();
    
    const int num_dof_elasto = num_soft_elasto * 4;
    
    W.resize(num_dof_elasto, system_size);
    
    tri_W.resize(num_soft_elasto * 3 * 27);
    
    const std::vector< VectorXi >& node_indices_x = scene.getNodeIndicesX();
    const std::vector< VectorXi >& node_indices_y = scene.getNodeIndicesY();
    const std::vector< VectorXi >& node_indices_z = scene.getNodeIndicesZ();
    
    threadutils::for_each(0, num_soft_elasto, [&] (int pidx) {
        auto& indices_x = scene.getParticleNodesX(pidx);
        auto& indices_y = scene.getParticleNodesY(pidx);
        auto& indices_z = scene.getParticleNodesZ(pidx);
        
        auto& grads_x = scene.getParticleGradsX(pidx);
        auto& grads_y = scene.getParticleGradsY(pidx);
        auto& grads_z = scene.getParticleGradsZ(pidx);
        
        auto& weights = scene.getParticleWeights(pidx);
        const Vector3s& part_pos = scene.getX().segment<3>(pidx * 4);
        
        for(int i = 0; i < 27; ++i) {
            if(!indices_x(i, 2) || weights(i, 0) == 0.0) {
                tri_W[(pidx * 3 + 0) * 27 + i] = Triplets(0, 0, 0.0);
            } else {
                const int global_idx = node_global_indices_x[indices_x(i, 0)][indices_x(i, 1)];
                
                if(global_idx == -1) {
                    tri_W[(pidx * 3 + 0) * 27 + i] = Triplets(0, 0, 0.0);
                } else {
                    tri_W[(pidx * 3 + 0) * 27 + i] = Triplets(pidx * 4 + 0, global_idx, weights(i, 0));
                }
            }
            
            if(!indices_y(i, 2) || weights(i, 1) == 0.0) {
                tri_W[(pidx * 3 + 1) * 27 + i] = Triplets(0, 0, 0.0);
            } else {
                const int global_idx = node_global_indices_y[indices_y(i, 0)][indices_y(i, 1)];
                
                if(global_idx == -1) {
                    tri_W[(pidx * 3 + 1) * 27 + i] = Triplets(0, 0, 0.0);
                } else {
                    tri_W[(pidx * 3 + 1) * 27 + i] = Triplets(pidx * 4 + 1, global_idx, weights(i, 1));
                }
            }
            
            if(!indices_z(i, 2) || weights(i, 2) == 0.0) {
                tri_W[(pidx * 3 + 2) * 27 + i] = Triplets(0, 0, 0.0);
            } else {
                const int global_idx = node_global_indices_z[indices_z(i, 0)][indices_z(i, 1)];
                
                if(global_idx == -1) {
                    tri_W[(pidx * 3 + 2) * 27 + i] = Triplets(0, 0, 0.0);
                } else {
                    tri_W[(pidx * 3 + 2) * 27 + i] = Triplets(pidx * 4 + 2, global_idx, weights(i, 2));
                }
            }
        }
    });
    
    tri_W.erase(std::remove_if(tri_W.begin(), tri_W.end(), [] (const auto& info) {return info.value() == 0.0;}), tri_W.end());
    
    W.reserve(tri_W.size());
    
    W.setFromTriplets(tri_W.begin(), tri_W.end());
}

void SceneStepper::mapSoftParticlesToNode( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z, const VectorXs& part_vec ) const
{
	const Sorter& buckets = scene.getParticleBuckets();
    const VectorXs& x = scene.getX();
	
	if((int) node_vec_x.size() != buckets.size()) node_vec_x.resize(buckets.size());
	if((int) node_vec_y.size() != buckets.size()) node_vec_y.resize(buckets.size());
	if((int) node_vec_z.size() != buckets.size()) node_vec_z.resize(buckets.size());
	
	const scalar invD = scene.getInverseDCoeff();
	
	const int kernel_order = scene.getKernelOrder();
    
    const std::vector<int>& particle_to_surfels = scene.getParticleToSurfels();
    const int num_elasto = scene.getNumElastoParticles();
    const std::vector< Matrix27x4s >& particle_weights = scene.getParticleWeights();
	
	buckets.for_each_bucket([&] (int bucket_idx) {
		const int num_nodes_x = scene.getNumNodesX(bucket_idx);
		VectorXs& bucket_node_vec_x = node_vec_x[ bucket_idx ];
		
		const int num_nodes_y = scene.getNumNodesY(bucket_idx);
		VectorXs& bucket_node_vec_y = node_vec_y[ bucket_idx ];
		
		const int num_nodes_z = scene.getNumNodesZ(bucket_idx);
		VectorXs& bucket_node_vec_z = node_vec_z[ bucket_idx ];
		
		if(bucket_node_vec_x.size() != num_nodes_x) bucket_node_vec_x.resize( num_nodes_x );
		if(bucket_node_vec_y.size() != num_nodes_y) bucket_node_vec_y.resize( num_nodes_y );
		if(bucket_node_vec_z.size() != num_nodes_z) bucket_node_vec_z.resize( num_nodes_z );
		
		const VectorXs& node_pos_x = scene.getNodePosX( bucket_idx );
		const VectorXs& node_pos_y = scene.getNodePosY( bucket_idx );
		const VectorXs& node_pos_z = scene.getNodePosZ( bucket_idx );
		
		for(int i = 0; i < num_nodes_x; ++i) {
			auto& particle_indices = scene.getNodeParticlePairsX(bucket_idx, i);
			const Vector3s& np = node_pos_x.segment<3>(i * 3);
			
            scalar ret(0.);
			for(auto& pp : particle_indices)
			{
				if(pp.first >= num_elasto || particle_to_surfels[pp.first] >= 0) continue;
				
				const auto& weights = particle_weights[pp.first];

				ret += part_vec( pp.first * 4 + 0 ) * weights(pp.second, 0);
			}
            
            bucket_node_vec_x(i) = ret;
		}
		
		for(int i = 0; i < num_nodes_y; ++i) {
			auto& particle_indices = scene.getNodeParticlePairsY(bucket_idx, i);
			const Vector3s& np = node_pos_y.segment<3>(i * 3);
			
            scalar ret(0.);
			for(auto& pp : particle_indices)
			{
                if(pp.first >= num_elasto || particle_to_surfels[pp.first] >= 0) continue;
				
				const auto& weights = particle_weights[pp.first];

				ret += part_vec( pp.first * 4 + 1 ) * weights(pp.second, 1);
			}
            bucket_node_vec_y(i) = ret;
		}
		
		for(int i = 0; i < num_nodes_z; ++i) {
			auto& particle_indices = scene.getNodeParticlePairsZ(bucket_idx, i);
			const Vector3s& np = node_pos_z.segment<3>(i * 3);
			
            scalar ret(0.);
			for(auto& pp : particle_indices)
			{
                if(pp.first >= num_elasto || particle_to_surfels[pp.first] >= 0) continue;
				
				const auto& weights = particle_weights[pp.first];

				ret += part_vec( pp.first * 4 + 2 ) * weights(pp.second, 2);
			}
            bucket_node_vec_z(i) = ret;
		}
	});
}


void SceneStepper::mapSoftParticlesToNodeSqr( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z, const VectorXs& part_vec ) const
{
    const Sorter& buckets = scene.getParticleBuckets();
    const VectorXs& x = scene.getX();
    
    if((int) node_vec_x.size() != buckets.size()) node_vec_x.resize(buckets.size());
    if((int) node_vec_y.size() != buckets.size()) node_vec_y.resize(buckets.size());
    if((int) node_vec_z.size() != buckets.size()) node_vec_z.resize(buckets.size());
    
    const scalar invD = scene.getInverseDCoeff();
    
    const int kernel_order = scene.getKernelOrder();
    
    const std::vector<int>& particle_to_surfels = scene.getParticleToSurfels();
    const int num_elasto = scene.getNumElastoParticles();
    const std::vector< Matrix27x4s >& particle_weights = scene.getParticleWeights();
    
    buckets.for_each_bucket([&] (int bucket_idx) {
        const int num_nodes_x = scene.getNumNodesX(bucket_idx);
        VectorXs& bucket_node_vec_x = node_vec_x[ bucket_idx ];
        
        const int num_nodes_y = scene.getNumNodesY(bucket_idx);
        VectorXs& bucket_node_vec_y = node_vec_y[ bucket_idx ];
        
        const int num_nodes_z = scene.getNumNodesZ(bucket_idx);
        VectorXs& bucket_node_vec_z = node_vec_z[ bucket_idx ];
        
        if(bucket_node_vec_x.size() != num_nodes_x) bucket_node_vec_x.resize( num_nodes_x );
        if(bucket_node_vec_y.size() != num_nodes_y) bucket_node_vec_y.resize( num_nodes_y );
        if(bucket_node_vec_z.size() != num_nodes_z) bucket_node_vec_z.resize( num_nodes_z );
        
        const VectorXs& node_pos_x = scene.getNodePosX( bucket_idx );
        const VectorXs& node_pos_y = scene.getNodePosY( bucket_idx );
        const VectorXs& node_pos_z = scene.getNodePosZ( bucket_idx );
        
        for(int i = 0; i < num_nodes_x; ++i) {
            auto& particle_indices = scene.getNodeParticlePairsX(bucket_idx, i);
            const Vector3s& np = node_pos_x.segment<3>(i * 3);
            
            scalar ret(0.);
            for(auto& pp : particle_indices)
            {
                if(pp.first >= num_elasto || particle_to_surfels[pp.first] >= 0) continue;
                
                const auto& weights = particle_weights[pp.first];

                scalar w = weights(pp.second, 0);
                ret += part_vec( pp.first * 4 + 0 ) * (w * w);
            }
            
            bucket_node_vec_x(i) = ret;
        }
        
        for(int i = 0; i < num_nodes_y; ++i) {
            auto& particle_indices = scene.getNodeParticlePairsY(bucket_idx, i);
            const Vector3s& np = node_pos_y.segment<3>(i * 3);
            
            scalar ret(0.);
            for(auto& pp : particle_indices)
            {
                if(pp.first >= num_elasto || particle_to_surfels[pp.first] >= 0) continue;
                
                const auto& weights = particle_weights[pp.first];

                scalar w = weights(pp.second, 1);
                ret += part_vec( pp.first * 4 + 1 ) * (w * w);
            }
            bucket_node_vec_y(i) = ret;
        }
        
        for(int i = 0; i < num_nodes_z; ++i) {
            auto& particle_indices = scene.getNodeParticlePairsZ(bucket_idx, i);
            const Vector3s& np = node_pos_z.segment<3>(i * 3);
            
            scalar ret(0.);
            for(auto& pp : particle_indices)
            {
                if(pp.first >= num_elasto || particle_to_surfels[pp.first] >= 0) continue;
                
                const auto& weights = particle_weights[pp.first];

                scalar w = weights(pp.second, 2);
                ret += part_vec( pp.first * 4 + 2 ) * (w * w);
            }
            bucket_node_vec_z(i) = ret;
        }
    });
}


void SceneStepper::mapGaussToNode( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z, const MatrixXs& gauss_vec ) const
{
	const Sorter& g_buckets = scene.getGaussBuckets();
	
	//	std::cout << "gauss_vec: \n" << gauss_vec << std::endl;
	
	g_buckets.for_each_bucket_particles_colored([&] (int pidx, int bucket_idx) {
		auto& indices_x = scene.getGaussNodesX(pidx);
		auto& grads_x = scene.getGaussGradsX(pidx);
		auto& indices_y = scene.getGaussNodesY(pidx);
		auto& grads_y = scene.getGaussGradsY(pidx);
		auto& indices_z = scene.getGaussNodesZ(pidx);
		auto& grads_z = scene.getGaussGradsZ(pidx);
		
		for(int i = 0; i < indices_x.rows(); ++i)
		{
			const int node_bucket_idx = indices_x(i, 0);
			const int node_idx = indices_x(i, 1);
            
            if(!indices_x(i, 2)) continue;
            
			node_vec_x[node_bucket_idx](node_idx) += gauss_vec.block<1, 3>(pidx * 3 + 0, 0).dot(grads_x.row(i).transpose());
		}
		
		for(int i = 0; i < indices_y.rows(); ++i)
		{
			const int node_bucket_idx = indices_y(i, 0);
			const int node_idx = indices_y(i, 1);
            
            if(!indices_y(i, 2)) continue;
            
			node_vec_y[node_bucket_idx](node_idx) += gauss_vec.block<1, 3>(pidx * 3 + 1, 0).dot(grads_y.row(i).transpose());
		}
		
		for(int i = 0; i < indices_z.rows(); ++i)
		{
			const int node_bucket_idx = indices_z(i, 0);
			const int node_idx = indices_z(i, 1);
            
            if(!indices_z(i, 2)) continue;
            
			node_vec_z[node_bucket_idx](node_idx) += gauss_vec.block<1, 3>(pidx * 3 + 2, 0).dot(grads_z.row(i).transpose());
		}
	}, scene.getNumBucketColors());
}

void SceneStepper::allocateCenterNodeVectors( const TwoDScene& scene, std::vector< VectorXi >& node_vec_p ) const
{
    const int num_buckets = scene.getNumBuckets();
    
    if((int) node_vec_p.size() != num_buckets) node_vec_p.resize(num_buckets);
    
    const Sorter& buckets = scene.getParticleBuckets();
    
    buckets.for_each_bucket([&] (int bucket_idx) {
        const int num_node_p = scene.getNumNodesP(bucket_idx);
        node_vec_p[bucket_idx].resize(num_node_p);
        node_vec_p[bucket_idx].setZero();
    });
}

void SceneStepper::allocateCenterNodeVectors( const TwoDScene& scene, std::vector< VectorXs >& node_vec_p ) const
{
    const int num_buckets = scene.getNumBuckets();
    
    if((int) node_vec_p.size() != num_buckets) node_vec_p.resize(num_buckets);
    
    const Sorter& buckets = scene.getParticleBuckets();
    
    buckets.for_each_bucket([&] (int bucket_idx) {
        const int num_node_p = scene.getNumNodesP(bucket_idx);
        node_vec_p[bucket_idx].resize(num_node_p);
        node_vec_p[bucket_idx].setZero();
    });
}

void SceneStepper::allocateNodeVectors( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z ) const
{
	const int num_buckets = scene.getNumBuckets();
	
	if((int) node_vec_x.size() != num_buckets) node_vec_x.resize(num_buckets);
	if((int) node_vec_y.size() != num_buckets) node_vec_y.resize(num_buckets);
	if((int) node_vec_z.size() != num_buckets) node_vec_z.resize(num_buckets);
	
	const Sorter& buckets = scene.getParticleBuckets();
	
	buckets.for_each_bucket([&] (int bucket_idx) {
		const int num_node_x = scene.getNumNodesX(bucket_idx);
		node_vec_x[bucket_idx].resize(num_node_x);
		node_vec_x[bucket_idx].setZero();
		
		const int num_node_y = scene.getNumNodesY(bucket_idx);
		node_vec_y[bucket_idx].resize(num_node_y);
		node_vec_y[bucket_idx].setZero();
		
		const int num_node_z = scene.getNumNodesZ(bucket_idx);
		node_vec_z[bucket_idx].resize(num_node_z);
		node_vec_z[bucket_idx].setZero();
	});
}

void SceneStepper::allocateNodeVectors( const TwoDScene& scene, std::vector< VectorXi >& node_vec_x, std::vector< VectorXi >& node_vec_y, std::vector< VectorXi >& node_vec_z ) const
{
    const int num_buckets = scene.getNumBuckets();
    
    if((int) node_vec_x.size() != num_buckets) node_vec_x.resize(num_buckets);
    if((int) node_vec_y.size() != num_buckets) node_vec_y.resize(num_buckets);
    if((int) node_vec_z.size() != num_buckets) node_vec_z.resize(num_buckets);
    
    const Sorter& buckets = scene.getParticleBuckets();
    
    buckets.for_each_bucket([&] (int bucket_idx) {
        const int num_node_x = scene.getNumNodesX(bucket_idx);
        node_vec_x[bucket_idx].resize(num_node_x);
        node_vec_x[bucket_idx].setZero();
        
        const int num_node_y = scene.getNumNodesY(bucket_idx);
        node_vec_y[bucket_idx].resize(num_node_y);
        node_vec_y[bucket_idx].setZero();
        
        const int num_node_z = scene.getNumNodesZ(bucket_idx);
        node_vec_z[bucket_idx].resize(num_node_z);
        node_vec_z[bucket_idx].setZero();
    });
}

scalar SceneStepper::dotNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az, const std::vector< VectorXs >& node_vec_bx, const std::vector< VectorXs >& node_vec_by, const std::vector< VectorXs >& node_vec_bz ) const
{
	VectorXs bucket_dot(node_vec_ax.size());
	
	assert(node_vec_ax.size() == node_vec_bx.size());
	assert(node_vec_ay.size() == node_vec_by.size());
	assert(node_vec_az.size() == node_vec_bz.size());
	
	const int num_buckets = node_vec_ax.size();
	
	threadutils::for_each(0, num_buckets, [&] (int bucket_idx){
		bucket_dot[bucket_idx] = node_vec_ax[bucket_idx].dot(node_vec_bx[bucket_idx]) + node_vec_ay[bucket_idx].dot(node_vec_by[bucket_idx]) + node_vec_az[bucket_idx].dot(node_vec_bz[bucket_idx]);
	});
	
	return bucket_dot.sum();
}

scalar SceneStepper::dotNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az, const std::vector< VectorXs >& node_vec_bx, const std::vector< VectorXs >& node_vec_by, const std::vector< VectorXs >& node_vec_bz, const VectorXs& twist_vec_a, const VectorXs& twist_vec_b ) const
{
    return dotNodeVectors(node_vec_ax, node_vec_ay, node_vec_az, node_vec_bx, node_vec_by, node_vec_bz) + twist_vec_a.dot(twist_vec_b);
}

scalar SceneStepper::dotNodeVectors( const std::vector< VectorXs >& node_vec_a, const std::vector< VectorXs >& node_vec_b ) const
{
	VectorXs bucket_dot(node_vec_a.size());
	
	assert(node_vec_a.size() == node_vec_b.size());
	
	const int num_buckets = node_vec_a.size();
	
	threadutils::for_each(0, num_buckets, [&] (int bucket_idx){
		bucket_dot[bucket_idx] = node_vec_a[bucket_idx].dot(node_vec_b[bucket_idx]);
	});
	
	return bucket_dot.sum();
}

scalar SceneStepper::lengthNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az, const VectorXs& twist_vec ) const
{
    VectorXs bucket_length(node_vec_ax.size());
    
    const int num_buckets = node_vec_ax.size();
    
    threadutils::for_each(0, num_buckets, [&] (int bucket_idx){
        bucket_length[bucket_idx] = node_vec_ax[bucket_idx].squaredNorm() + node_vec_ay[bucket_idx].squaredNorm() + node_vec_az[bucket_idx].squaredNorm();
    });
    
    return sqrt(bucket_length.sum() + twist_vec.squaredNorm());
}

scalar SceneStepper::lengthNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az ) const
{
	VectorXs bucket_length(node_vec_ax.size());
	
	const int num_buckets = node_vec_ax.size();
	
	threadutils::for_each(0, num_buckets, [&] (int bucket_idx){
		bucket_length[bucket_idx] = node_vec_ax[bucket_idx].squaredNorm() + node_vec_ay[bucket_idx].squaredNorm() + node_vec_az[bucket_idx].squaredNorm();
	});
	
	return sqrt(bucket_length.sum());
}

scalar SceneStepper::lengthNodeVectors( const std::vector< VectorXs >& node_vec ) const
{
	VectorXs bucket_length(node_vec.size());
	
	const int num_buckets = node_vec.size();
	
	threadutils::for_each(0, num_buckets, [&] (int bucket_idx){
		bucket_length[bucket_idx] = node_vec[bucket_idx].squaredNorm();
	});
	
	return sqrt(bucket_length.sum());
}


