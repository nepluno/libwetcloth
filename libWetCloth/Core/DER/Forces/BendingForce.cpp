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

#include "BendingForce.h"
#include "ViscousOrNotViscous.h"
#include "../Dependencies/ElasticStrandUtils.h"

template<typename ViscousT>
void BendingForce<ViscousT>::computeLocal( LocalMultiplierType& localL, const StrandForce& strand, const IndexType vtx, const scalar& dt )
{
    // B.ilen.(phi+hJv)
    const Mat2& B = ViscousT::bendingMatrix( strand, vtx );
    const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec2& kappa = strand.m_strandState->m_kappas[vtx];
    const GradKType& gradKappa = strand.m_strandState->m_gradKappas[vtx];
    const scalar psi_coeff = strand.m_packing_fraction[ vtx ];

    Vec2 Jv = gradKappa.transpose() * strand.m_v_plus.segment<11>(4 * (vtx - 1));

    localL = -ilen * psi_coeff * B * (kappa - kappaBar + dt * Jv);
}

template<typename ViscousT>
scalar BendingForce<ViscousT>::localEnergy( const StrandForce& strand, const IndexType vtx )
{
    const Mat2& B = ViscousT::bendingMatrix( strand, vtx );
    const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec2& kappa = strand.m_strandState->m_kappas[vtx];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];

    return 0.5 * ilen * psi_coeff * ( kappa - kappaBar ).dot( Vec2( B * ( kappa - kappaBar ) ) );
}

template<typename ViscousT>
void BendingForce<ViscousT>::computeLocal( Eigen::Matrix<scalar, 11, 1>& localF,
        const StrandForce& strand, const IndexType vtx )
{
    const Mat2& B = ViscousT::bendingMatrix( strand, vtx );
    const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec2& kappa = strand.m_strandState->m_kappas[vtx];
    const GradKType& gradKappa = strand.m_strandState->m_gradKappas[vtx];
	const scalar psi_coeff = strand.m_packing_fraction[ vtx ];
	
    localF = -ilen * psi_coeff * gradKappa * B * ( kappa - kappaBar );
}

template<typename ViscousT>
void BendingForce<ViscousT>::computeLocal(Eigen::Matrix<scalar, 11, 11>& localJ,
        const StrandForce& strand, const IndexType vtx )
{
    const scalar psi_coeff = strand.m_packing_fraction[ vtx ];
    const scalar ilen = strand.m_invVoronoiLengths[vtx];
    localJ = -ilen * ViscousT::bendingCoefficient( strand, vtx ) * psi_coeff * strand.m_strandState->m_bendingProducts[vtx];

#ifndef USE_APPROX_GRAD_KAPPA
    if( strand.m_requiresExactForceJacobian )
    {
        // const Mat2& bendingMatrixBase = strand.m_strandParams->bendingMatrixBase(vtx);
        // const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
        // const Vec2& kappa = strand.m_strandState->m_kappas[vtx];
        // const std::pair<LocalJacobianType, LocalJacobianType>& hessKappa = strand.m_strandState->m_hessKappas[vtx];
        // const Vec2& temp = bendingMatrixBase * ( kappa - kappaBar );

        // localJ += temp( 0 ) * hessKappa.first + temp( 1 ) * hessKappa.second;

        const std::pair<LocalJacobianType, LocalJacobianType>& hessKappa = strand.m_strandState->m_hessKappas[vtx];
        const Vec2 temp = ViscousT::bendingMultiplier( strand, vtx );
        localJ += temp( 0 ) * hessKappa.first + temp( 1 ) * hessKappa.second;
    }
#endif
}

template<typename ViscousT>
void BendingForce<ViscousT>::addInPosition( VecX& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<11>( 4 * ( vtx - 1 ) ) += localForce;
}

template<typename ViscousT>
void BendingForce<ViscousT>::addInPosition( VecX& globalMultiplier, const IndexType vtx, const LocalMultiplierType& localL )
{
    globalMultiplier.segment<2>( 2 * vtx ) += localL;
}

template<typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentE( scalar& energy, StrandForce& strand )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        energy += localEnergy( strand, vtx );
    }
}

template<typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentF( VecX& force, StrandForce& strand )
{
    for( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalForceType localF;
        computeLocal( localF, strand, vtx );
        addInPosition( force, vtx, localF );
    }
}

template class BendingForce<NonViscous>;
template class BendingForce<Viscous>;
