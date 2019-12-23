//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "SimpleGravityForce.h"
#include "TwoDScene.h"
#include "ThreadUtils.h"
#include "sorter.h"

SimpleGravityForce::SimpleGravityForce( const Vector3s& gravity )
: Force()
, m_gravity(gravity)
{
	assert( (m_gravity.array()==m_gravity.array()).all() );
	assert( (m_gravity.array()!=std::numeric_limits<scalar>::infinity()).all() );
}

SimpleGravityForce::~SimpleGravityForce()
{}

void SimpleGravityForce::addLiquidGradEToNode( const TwoDScene& scene, std::vector< VectorXs >& node_rhs_x, std::vector< VectorXs >& node_rhs_y, std::vector< VectorXs >& node_rhs_z, const scalar& coeff )
{
    const std::vector< VectorXs >& fluid_mass_x = scene.getNodeFluidMassX();
    const std::vector< VectorXs >& fluid_mass_y = scene.getNodeFluidMassY();
    const std::vector< VectorXs >& fluid_mass_z = scene.getNodeFluidMassZ();
    
    const Sorter& buckets = scene.getParticleBuckets();
    buckets.for_each_bucket([&] (int bucket_idx) {
        const int num_nodes_x = node_rhs_x[bucket_idx].size();
        const int num_nodes_y = node_rhs_y[bucket_idx].size();
        const int num_nodes_z = node_rhs_z[bucket_idx].size();
        
        for(int i = 0; i < num_nodes_x; ++i) {
            node_rhs_x[bucket_idx][i] -= fluid_mass_x[bucket_idx][i] * m_gravity(0) * coeff;
        }
        
        for(int i = 0; i < num_nodes_y; ++i) {
            node_rhs_y[bucket_idx][i] -= fluid_mass_y[bucket_idx][i] * m_gravity(1) * coeff;
        }
        
        for(int i = 0; i < num_nodes_z; ++i) {
            node_rhs_z[bucket_idx][i] -= fluid_mass_z[bucket_idx][i] * m_gravity(2) * coeff;
        }
    });
}

void SimpleGravityForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
	assert( x.size() == v.size() );
	assert( x.size() == m.size() );
	assert( x.size()%4 == 0 );
	
	// Assume 0 potential is at origin
	for( int i = 0; i < x.size()/4; ++i ) E -= m(4*i)*m_gravity.dot(x.segment<3>(4*i));
}

void SimpleGravityForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{	
	const int num_elasto = gradE.size()/4;
	threadutils::for_each(0, num_elasto, [&] (int i) {
		gradE.segment<3>(4*i) -= m(4*i)*m_gravity;
	});
}

void SimpleGravityForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
	assert( x.size() == v.size() );
	assert( x.size() == m.size() );

	assert( x.size()%4 == 0 );
	// Nothing to do.
}

void SimpleGravityForce::updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt )
{
	
}

int SimpleGravityForce::numHessX( )
{
	return 0;
}

void SimpleGravityForce::preCompute()
{
}

void SimpleGravityForce::updateStartState()
{
}

Force* SimpleGravityForce::createNewCopy()
{
	return new SimpleGravityForce(*this);
}

int SimpleGravityForce::flag() const
{
	return 3;
}
