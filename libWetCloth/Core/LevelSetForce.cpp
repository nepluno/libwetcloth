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

#include "LevelSetForce.h"
#include "TwoDScene.h"
#include "MathUtilities.h"

LevelSetForce::LevelSetForce( const std::shared_ptr<TwoDScene>& scene, const scalar& l0, const scalar& b )
: Force()
, m_scene(scene)
, m_l0(std::min(scene->getCellSize() * 3.0, l0))
, m_b(b)
{
	assert( m_l0 >= 0.0 );
	assert( m_b >= 0.0 );
}

LevelSetForce::~LevelSetForce()
{}

void LevelSetForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
	assert( x.size() == v.size() );
	assert( x.size()%4 == 0 );
    std::cerr << "NOT IMPLEMENTED!" << std::endl;
}

void LevelSetForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{
	assert( x.size() == v.size() );
	assert( x.size() == gradE.size() );
	assert( x.size()%4 == 0 );
    
    const int num_elasto = m_process_list.size();
    const std::vector< VectorXs >& node_pos = m_scene->getNodePosSolidPhi();
    const std::vector< VectorXs >& solid_phi = m_scene->getNodeSolidPhi();
    
    const VectorXs& vol = m_scene->getVol();
    
    const scalar K = m_scene->getLiquidInfo().levelset_young_modulus;
    
    const scalar iD = m_scene->getInverseDCoeff();
    
    threadutils::for_each(0, num_elasto, [&] (int idx) {
        const int pidx = m_process_list[idx];
        const auto& node_indices_sphi = m_scene->getParticleNodesSolidPhi(pidx);
        const auto& particle_weights = m_scene->getParticleWeights(pidx);
        
        const Vector3s& pos = x.segment<3>(pidx * 4);
        
        scalar phi_ori = 0.0;
        Vector3s grad_phi = Vector3s::Zero();
        
        for(int nidx = 0; nidx < node_indices_sphi.rows(); ++nidx)
        {
            const int bucket_idx = node_indices_sphi(nidx, 0);
            const int node_idx = node_indices_sphi(nidx, 1);
            
            scalar phi;
            if(node_indices_sphi(nidx, 2)) {
                phi = solid_phi[bucket_idx](node_idx);
            } else {
                phi = 3.0 * m_scene->getCellSize();
            }
            
            const scalar w = particle_weights(nidx, 3);
            const Vector3s& np = node_pos[bucket_idx].segment<3>(node_idx * 3);
            
            phi_ori += phi * w;
            grad_phi += phi * iD * w * (np - pos);
        }
        
        if(grad_phi.norm() > 1e-12) grad_phi.normalize();
        
        phi_ori -= m_l0;
        
        if(phi_ori < 0.0) {
            const scalar k = K * pow(vol(pidx), 1. / 3.);
            
            gradE.segment<3>(pidx * 4) += k * phi_ori * grad_phi;
        }
    });
}

bool LevelSetForce::parallelized() const
{
    return true;
}

void LevelSetForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
	assert( x.size() == v.size() );
	assert( x.size() == m.size() );
	assert( x.size()%4 == 0 );
    
    const int num_elasto = m_process_list.size();
    const std::vector< VectorXs >& node_pos = m_scene->getNodePosSolidPhi();
    const std::vector< VectorXs >& solid_phi = m_scene->getNodeSolidPhi();
    
    const VectorXs& vol = m_scene->getVol();
    
    const scalar K = m_scene->getLiquidInfo().levelset_young_modulus;

    const scalar iD = m_scene->getInverseDCoeff();
    
    threadutils::for_each(0, num_elasto, [&] (int idx) {
        const int pidx = m_process_list[idx];
        const auto& node_indices_sphi = m_scene->getParticleNodesSolidPhi(pidx);
        const auto& particle_weights = m_scene->getParticleWeights(pidx);
        
        scalar phi_ori = 0.0;
        Vector3s grad_phi = Vector3s::Zero();
        
        const Vector3s& pos = x.segment<3>(pidx * 4);
        
        for(int nidx = 0; nidx < node_indices_sphi.rows(); ++nidx)
        {
            const int bucket_idx = node_indices_sphi(nidx, 0);
            const int node_idx = node_indices_sphi(nidx, 1);
            
            scalar phi;
            if(node_indices_sphi(nidx, 2)) {
                phi = solid_phi[bucket_idx](node_idx);
            } else {
                phi = 3.0 * m_scene->getCellSize();
            }
            
            const scalar w = particle_weights(nidx, 3);
            const Vector3s& np = node_pos[bucket_idx].segment<3>(node_idx * 3);
            
            phi_ori += phi * w;
            grad_phi += phi * iD * w * (np - pos);
        }
        
        if(grad_phi.norm() > 1e-12) grad_phi.normalize();
        
        phi_ori -= m_l0;
        
        if(phi_ori < 0.0) {
            const scalar k = K * pow(vol(pidx), 1. / 3.);
            
            Matrix3s hess = k * grad_phi * grad_phi.transpose();
            
            for(int s = 0; s < 3; ++s) for(int r = 0; r < 3; ++r) {
                hessE[hessE_index + idx * 9 + s * 3 + r] = Triplets(pidx * 4 + r, pidx * 4 + s, hess(r, s));
            }
        } else {
            for(int s = 0; s < 3; ++s) for(int r = 0; r < 3; ++r) {
                hessE[hessE_index + idx * 9 + s * 3 + r] = Triplets(pidx * 4 + r, pidx * 4 + s, 0.0);
            }
        }
    });
}

void LevelSetForce::updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt )
{

}

void LevelSetForce::preCompute()
{
}

void LevelSetForce::updateStartState()
{
}

int LevelSetForce::numHessX()
{
    const int num_elasto = m_scene->getNumSoftElastoParticles();
    const VectorXs& x = m_scene->getX();
    m_process_list.resize(0);
    m_process_list.reserve(num_elasto);
    
    const scalar dx = m_scene->getCellSize();
    
    for(int i = 0; i < num_elasto; ++i)
    {
		if(m_scene->isFixed(i) & 1) continue;
		
        Vector3s vel;
        scalar phi = m_scene->computePhiVel(x.segment<3>(i * 4), vel);
        
        if(phi < m_l0) m_process_list.push_back(i);
    }
    
	return 9 * (int) m_process_list.size();
}

Force* LevelSetForce::createNewCopy()
{
	return new LevelSetForce(*this);
}

int LevelSetForce::flag() const
{
	return 1;
}
