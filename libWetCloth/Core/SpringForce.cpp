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


#include "SpringForce.h"

SpringForce::SpringForce( const Vector2iT& endpoints, const scalar& k, const scalar& l0, const scalar& b )
: Force()
, m_endpoints(endpoints)
, m_k(k)
, m_l0(l0)
, m_b(b)
, m_zero_length(l0 == 0.0)
{
	assert( m_k >= 0.0 );
	assert( m_l0 >= 0.0 );
	assert( m_b >= 0.0 );
}

SpringForce::~SpringForce()
{}

void SpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
	assert( x.size() == v.size() );
	assert( x.size()%4 == 0 );
	
	scalar l = (x.segment<3>(4*m_endpoints(1))-x.segment<3>(4*m_endpoints(0))).norm();
	scalar psi_coeff = pow((psi(m_endpoints(1))+psi(m_endpoints(0))) * 0.5, lambda);
	E += 0.5*m_k*(l-m_l0)*(l-m_l0) * psi_coeff;
}

void SpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{
	assert( x.size() == v.size() );
	assert( x.size() == gradE.size() );
	assert( x.size()%4 == 0 );
	scalar psi_coeff = pow((psi(m_endpoints(1))+psi(m_endpoints(0))) * 0.5, lambda);
	
	if(m_zero_length) {
		Vector3s nhat = x.segment<3>(4*m_endpoints(1))-x.segment<3>(4*m_endpoints(0));
		nhat *= m_k * psi_coeff;
		gradE.segment<3>(4 * m_endpoints(0)) -= nhat;
		gradE.segment<3>(4 * m_endpoints(1)) += nhat;
	} else {
	// Compute the elastic component
		Vector3s nhat = x.segment<3>(4*m_endpoints(1))-x.segment<3>(4*m_endpoints(0));
		scalar l = nhat.norm();
		assert( l != 0.0 );
		nhat /= l;
		nhat *= m_k*(l-m_l0)*psi_coeff;
		gradE.segment<3>(4 * m_endpoints(0)) -= nhat;
		gradE.segment<3>(4 * m_endpoints(1)) += nhat;
	}
}

void SpringForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
	assert( x.size() == v.size() );
	assert( x.size() == m.size() );
	assert( x.size()%4 == 0 );
	scalar psi_coeff = pow((psi(m_endpoints(1)) + psi(m_endpoints(0))) * 0.5, lambda);
	
	// Contribution from elastic component
	if(m_zero_length) {
		const scalar coeff = m_k * psi_coeff;
		for(int r = 0; r < 3; ++r) {
			hessE[hessE_index + r * 4 + 0] = Triplets(4 * m_endpoints(0) + r, 4 * m_endpoints(0) + r, coeff);
			hessE[hessE_index + r * 4 + 1] = Triplets(4 * m_endpoints(1) + r, 4 * m_endpoints(1) + r, coeff);
			hessE[hessE_index + r * 4 + 2] = Triplets(4 * m_endpoints(0) + r, 4 * m_endpoints(1) + r, -coeff);
			hessE[hessE_index + r * 4 + 3] = Triplets(4 * m_endpoints(1) + r, 4 * m_endpoints(0) + r, -coeff);
		}
	} else {
		Vector3s nhat = x.segment<3>(4 * m_endpoints(1)) - x.segment<3>(4 * m_endpoints(0));
		scalar l = nhat.norm();
		assert( l != 0 );
		nhat /= l;
		
		Matrix3s hess;
		hess = nhat * nhat.transpose();
		hess += (l - m_l0) * (Matrix3s::Identity() - hess) / l;
		hess *= m_k * psi_coeff;
		
		for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) {
			int element_idx = r * 3 + s;
			hessE[hessE_index + element_idx * 4 + 0] = Triplets(4 * m_endpoints(0) + r, 4 * m_endpoints(0) + s, hess(r, s));
			hessE[hessE_index + element_idx * 4 + 1] = Triplets(4 * m_endpoints(1) + r, 4 * m_endpoints(1) + s, hess(r, s));
			hessE[hessE_index + element_idx * 4 + 2] = Triplets(4 * m_endpoints(0) + r, 4 * m_endpoints(1) + s, -hess(r, s));
			hessE[hessE_index + element_idx * 4 + 3] = Triplets(4 * m_endpoints(1) + r, 4 * m_endpoints(0) + s, -hess(r, s));
		}
	}
}

void SpringForce::updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt )
{
	
}

void SpringForce::preCompute()
{
}

void SpringForce::updateStartState()
{
}

int SpringForce::numHessX()
{
	if(m_zero_length) {
		return 3 * 4;
	} else {
		return 9 * 4;
	}
}

Force* SpringForce::createNewCopy()
{
	return new SpringForce(*this);
}

int SpringForce::flag() const
{
	return 1;
}
