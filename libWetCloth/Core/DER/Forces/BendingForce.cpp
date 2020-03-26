//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
    if ( strand.m_requiresExactForceJacobian )
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
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        energy += localEnergy( strand, vtx );
    }
}

template<typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentF( VecX& force, StrandForce& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalForceType localF;
        computeLocal( localF, strand, vtx );
        addInPosition( force, vtx, localF );
    }
}

template class BendingForce<NonViscous>;
template class BendingForce<Viscous>;
