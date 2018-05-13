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

#include "StretchingForce.h"
#include "ViscousOrNotViscous.h"

template<typename ViscousT>
scalar StretchingForce<ViscousT>::localEnergy( const StrandForce& strand, const IndexType vtx )
{
    const scalar ks = ViscousT::ks( strand, vtx );
    const scalar restLength = ViscousT::ellBar( strand, vtx );

    const scalar length = strand.m_strandState->m_lengths[vtx];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];

    return 0.5 * ks * square( length / restLength - 1.0 ) * restLength * psi_coeff;
}

template<typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(Eigen::Matrix<scalar, 6, 1>& localF, const StrandForce& strand, const IndexType vtx )
{
    const scalar ks = ViscousT::ks( strand, vtx );
    const scalar restLength = ViscousT::ellBar( strand, vtx );

    const scalar length = strand.m_strandState->m_lengths[vtx];
    const Vec3& edge = strand.m_strandState->m_tangents[ vtx ];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];

    Vec3 f = ks * ( length / restLength - 1.0 ) * edge * psi_coeff;

    localF.segment<3>( 0 ) = f;
    localF.segment<3>( 3 ) = -f;
}

template<typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(Eigen::Matrix<scalar, 6, 6>& localJ,
        const StrandForce& strand, const IndexType vtx )
{
    const scalar ks = ViscousT::ks( strand, vtx );
    const scalar restLength = ViscousT::ellBar( strand, vtx );

    const scalar length = strand.m_strandState->m_lengths[vtx];
    const Vec3& edge = strand.m_strandState->m_tangents[ vtx ];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];

    bool useApprox = !strand.m_requiresExactForceJacobian && length < restLength;

    Mat3 M ;
    if( useApprox ){
        M = ks * psi_coeff / restLength * ( edge * edge.transpose() );
    }
    else{
		M = ks * psi_coeff / restLength * (edge * edge.transpose() + (length - restLength) * (Mat3::Identity() - edge * edge.transpose()) / length);
    }

    localJ.block<3, 3>( 0, 0 ) = localJ.block<3, 3>( 3, 3 ) = -M;
    localJ.block<3, 3>( 0, 3 ) = localJ.block<3, 3>( 3, 0 ) = M;
}

template<typename ViscousT>
void StretchingForce<ViscousT>::addInPosition( ForceVectorType& globalForce, const IndexType vtx, const LocalForceType& localForce )
{
    globalForce.segment<3>( 4 * vtx ) += localForce.segment<3>( 0 );
    globalForce.segment<3>( 4 * ( vtx + 1 ) ) += localForce.segment<3>( 3 );
}

template<typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentE( scalar& energy, StrandForce& strand )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        energy += localEnergy( strand, vtx );
    }
}

template<typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentF( VecX& force, StrandForce& strand )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalForceType localF;
        computeLocal( localF, strand, vtx );
        addInPosition( force, vtx, localF );
    }
}

template<typename ViscousT>
void StretchingForce<ViscousT>::accumulateIntegrationVars( 
    const unsigned& pos_start, const unsigned& j_start, const unsigned& tildek_start, const unsigned& global_start_dof, 
    StrandForce& strand, VectorXs& lambda, TripletXs& J, TripletXs& tildeK, TripletXs& stiffness, VectorXs& Phi, const int& lambda_start )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        unsigned dfirst = global_start_dof + 4 * vtx;
        unsigned dsecond = global_start_dof + 4 * (vtx+1);

        const scalar length = strand.m_strandState->m_lengths[vtx];
        unsigned idx_pos = pos_start + vtx;
        Phi[ idx_pos ] = length - ViscousT::ellBar( strand, vtx );
        stiffness[ idx_pos ] = Triplets( idx_pos, idx_pos, ViscousT::ks( strand, vtx ) );

        const Vec3& edge = strand.m_strandState->m_tangents[ vtx ];
        unsigned idx_j = j_start + 6 * vtx;
        for( int r = 0; r < 3; ++r ){
            J[ idx_j + r ] = Triplets( idx_pos, dfirst + r, -edge(r));
            J[ idx_j + 3 + r ] = Triplets( idx_pos, dsecond + r, edge(r));
        }


        if( std::is_same< ViscousT, NonViscous >::value )
        { // only fill in tildeK once    
            Mat3 M = ( Mat3::Identity() - edge * edge.transpose()) / length;
            scalar weight = -lambda[ lambda_start + vtx ];
            
            unsigned idx_tildek = tildek_start + 36 * vtx;
            for(int r = 0; r < 3; ++r) {
                for(int s = 0; s < 3; ++s) {
                    tildeK[ idx_tildek + r * 6 + s] = Triplets( dfirst + r, dfirst + s, M(r, s) * weight);
                    tildeK[ idx_tildek + (r + 3) * 6 + 3 + s] = Triplets( dsecond + r, dsecond + s, M(r, s) * weight);

                    tildeK[ idx_tildek + r * 6 + 3 + s] = Triplets( dfirst + r, dsecond + s, -M(r, s) * weight);
                    tildeK[ idx_tildek + (r + 3) * 6 + s] = Triplets( dsecond + r, dfirst + s, -M(r, s) * weight);
                }
            }
        }

    }
}

template class StretchingForce<NonViscous>;
template class StretchingForce<Viscous>;
