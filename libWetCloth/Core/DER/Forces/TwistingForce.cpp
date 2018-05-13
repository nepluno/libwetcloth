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
        const scalar undeformedTwist = ViscousT::thetaBar( strand, vtx );
        const scalar twist = strand.m_strandState->m_twists[vtx];
        const Mat11& hessTwist = strand.m_strandState->m_hessTwists[vtx];
        localJ += -kt * ilen * ( twist - undeformedTwist ) * hessTwist * psi_coeff;
    }
}

template<typename ViscousT>
void TwistingForce<ViscousT>::addInPosition( VecX& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<11>( 4 * ( vtx - 1 ) ) += localForce;
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

template<typename ViscousT>
void TwistingForce<ViscousT>::accumulateIntegrationVars( 
    const unsigned& pos_start, const unsigned& j_start, const unsigned& tildek_start, const unsigned& global_start_dof, 
    StrandForce& strand, VectorXs& lambda, TripletXs& J, TripletXs& tildeK, TripletXs& stiffness, VectorXs& Phi, const int& lambda_start )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        unsigned dfirst = global_start_dof + 4 * (vtx-1);
        unsigned dsecond = global_start_dof + 4 * vtx;
        unsigned dthird = global_start_dof + 4 * (vtx+1);

        const scalar twist = strand.m_strandState->m_twists[vtx];
        const scalar undeformedTwist = ViscousT::thetaBar( strand, vtx );
        const scalar kt = ViscousT::kt( strand, vtx );
        const scalar ilen = strand.m_invVoronoiLengths[vtx];

        unsigned idx_pos = pos_start + vtx - 1;
        Phi[ idx_pos ] = twist - undeformedTwist;
        stiffness[ idx_pos ] = Triplets( idx_pos, idx_pos, kt * ilen );

        const Vec11& gt = strand.m_strandState->m_gradTwists[vtx];
        unsigned idx_j = j_start + 11 * (vtx - 1);
        for( int r = 0; r < 4; ++r ){
            J[ idx_j + r ] = Triplets( idx_pos, dfirst + r, gt(r));
            J[ idx_j + 4 + r ] = Triplets( idx_pos, dsecond + r, gt(r + 4));
            if( r < 3 ) J[ idx_j + 8 + r ] = Triplets( idx_pos, dthird + r, gt(r + 8));
        }

        if( std::is_same< ViscousT, NonViscous >::value ){
            scalar weight = -lambda[ lambda_start + (vtx-1) ];
            const Mat11& hessTwist = strand.m_strandState->m_hessTwists[vtx];
            unsigned idx_tildek = tildek_start + 121 * (vtx - 1);
            for(int r = 0; r < 11; ++r) {
                for(int s = 0; s < 11; ++s){
                  if(r == 3 || r == 7 || s == 3 || s == 7) {
                    tildeK[ idx_tildek + r * 11 + s] = Triplets( dfirst + r, dfirst + s, hessTwist(r, s) * kt * ilen * (twist - undeformedTwist));
                  } else {
                    tildeK[ idx_tildek + r * 11 + s] = Triplets( dfirst + r, dfirst + s, hessTwist(r, s) * weight);
                  }
                }
            }
        }
    }
}

template class TwistingForce<NonViscous>;
template class TwistingForce<Viscous>;
