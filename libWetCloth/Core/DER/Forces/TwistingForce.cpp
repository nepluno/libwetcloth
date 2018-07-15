//
// This file is part of the libWetHair open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// The Clear BSD License
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
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

#include "TwistingForce.h"
#include "ViscousOrNotViscous.h"

template<typename ViscousT>
void TwistingForce<ViscousT>::computeLocal( LocalMultiplierType& localL, const StrandForce& strand, const IndexType vtx, const scalar& dt )
{
    const scalar kt = ViscousT::kt( strand, vtx );
    const scalar undefTwist = ViscousT::thetaBar( strand, vtx );
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    const scalar twist = strand.m_strandState->m_twists[vtx];
    const scalar psi_coeff = strand.m_packing_fraction[ vtx ];

    localL = -kt * ilen * psi_coeff * (twist - undefTwist + dt * strand.m_strandState->m_gradTwists[vtx].dot(strand.m_v_plus.segment<11>(4 * (vtx - 1))));    
}

template<typename ViscousT>
scalar TwistingForce<ViscousT>::localEnergy( const StrandForce& strand, const IndexType vtx )
{
    const scalar kt = ViscousT::kt( strand, vtx );
    const scalar undefTwist = ViscousT::thetaBar( strand, vtx );
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    const scalar twist = strand.m_strandState->m_twists[vtx];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];
	
    return 0.5 * kt * square( twist - undefTwist ) * ilen * psi_coeff;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::computeLocal(Eigen::Matrix<scalar, 11, 1>& localF,
        const StrandForce& strand, const IndexType vtx )
{
    const scalar kt = ViscousT::kt( strand, vtx );
    const scalar undefTwist = ViscousT::thetaBar( strand, vtx );
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    const scalar twist = strand.m_strandState->m_twists[vtx];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];
	
    localF = -kt * ilen * ( twist - undefTwist ) * strand.m_strandState->m_gradTwists[vtx] * psi_coeff;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::computeLocal( Eigen::Matrix<scalar, 11, 11>& localJ,
        const StrandForce& strand, const IndexType vtx )
{
    const scalar kt = ViscousT::kt( strand, vtx );
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Mat11& gradTwistSquared = strand.m_strandState->m_gradTwistsSquared[vtx];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];
	
    localJ = -kt * ilen * gradTwistSquared * psi_coeff;
    if( strand.m_requiresExactForceJacobian )
    {
        const Mat11& hessTwist = strand.m_strandState->m_hessTwists[vtx];
        const scalar localL = ViscousT::twistingMultiplier( strand, vtx );

        localJ += localL * hessTwist;
    }
}

template<typename ViscousT>
void TwistingForce<ViscousT>::addInPosition( VecX& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<11>( 4 * ( vtx - 1 ) ) += localForce;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::addInPosition( VecX& globalMultiplier, const IndexType vtx, const LocalMultiplierType& localL )
{
    globalMultiplier(vtx) += localL;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentE( scalar& energy, StrandForce& strand )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        energy += localEnergy( strand, vtx );
    }
}

template<typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentF( VecX& force, StrandForce& strand )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalForceType localF;
        computeLocal( localF, strand, vtx );
        addInPosition( force, vtx, localF );
    }
}

template class TwistingForce<NonViscous>;
template class TwistingForce<Viscous>;
