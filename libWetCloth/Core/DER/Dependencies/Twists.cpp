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

#include "Twists.h"
#include "ElasticStrandUtils.h"

void Twists::compute()
{
    m_value.resize( m_size );
    const std::vector<scalar>& refTwists = m_refTwists.get();
    const VecX& dofs = m_dofs.get();

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        m_value[vtx] = refTwists[vtx] + dofs[4 * vtx + 3] - dofs[4 * vtx - 1];
    }

    setDependentsDirty();
}

void GradTwists::compute()
{
    m_value.resize( m_size );
    const Vec3Array& curvatureBinormals = m_curvatureBinormals.get();
    const std::vector<scalar>& lengths = m_lengths.get();

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        Vec11& Dtwist = m_value[vtx];
		Dtwist.setZero();

        const Vec3& kb = curvatureBinormals[vtx];

        Dtwist.segment<3>( 0 ) = -0.5 / lengths[vtx - 1] * kb;
        Dtwist.segment<3>( 8 ) = 0.5 / lengths[vtx] * kb;
        Dtwist.segment<3>( 4 ) = -( Dtwist.segment<3>( 0 ) + Dtwist.segment<3>( 8 ) );
        Dtwist( 3 ) = -1;
        Dtwist( 7 ) = 1;
    }

    setDependentsDirty();
}

void GradTwistsSquared::compute()
{
    m_value.resize( m_size );
    const Vec11Array& gradTwists = m_gradTwists.get();

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        const Vec11& gradTwist = gradTwists[vtx];
        m_value[vtx] = gradTwist * gradTwist.transpose();
    }

    setDependentsDirty();
}

void HessTwists::compute()
{
    m_value.resize( m_size );
    const Vec3Array& tangents = m_tangents.get();
    const std::vector<scalar>& lengths = m_lengths.get();
    const Vec3Array& curvatureBinormals = m_curvatureBinormals.get();

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        Mat11& DDtwist = m_value[vtx];

        DDtwist.setZero();

        const Vec3& te = tangents[vtx - 1];
        const Vec3& tf = tangents[vtx];
        const scalar norm_e = lengths[vtx - 1];
        const scalar norm_f = lengths[vtx];
        const Vec3& kb = curvatureBinormals[vtx];

        scalar chi = 1 + te.dot( tf );

        //    assert( chi>0 );
        if( chi <= 0 )
        {
            std::cerr << "StrandState::computeHessTwist for state " << this
                    << " chi = " << chi << " te = " << te << " tf = " << tf << std::endl;
            chi = 1e-12;
        }

        const Vec3& tilde_t = 1.0 / chi * ( te + tf );

        const Mat3& D2mDe2 = -0.25 / square( norm_e )
                * ( outerProd<3>( kb, te + tilde_t ) + outerProd<3>( te + tilde_t, kb ) );
        const Mat3& D2mDf2 = -0.25 / square( norm_f )
                * ( outerProd<3>( kb, tf + tilde_t ) + outerProd<3>( tf + tilde_t, kb ) );
        const Mat3& D2mDeDf = 0.5 / ( norm_e * norm_f )
                * ( 2.0 / chi * crossMat( te ) - outerProd<3>( kb, tilde_t ) );
        const Mat3& D2mDfDe = D2mDeDf.transpose();

        DDtwist.block<3, 3>( 0, 0 ) = D2mDe2;
        DDtwist.block<3, 3>( 0, 4 ) = -D2mDe2 + D2mDeDf;
        DDtwist.block<3, 3>( 4, 0 ) = -D2mDe2 + D2mDfDe;
        DDtwist.block<3, 3>( 4, 4 ) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
        DDtwist.block<3, 3>( 0, 8 ) = -D2mDeDf;
        DDtwist.block<3, 3>( 8, 0 ) = -D2mDfDe;
        DDtwist.block<3, 3>( 8, 4 ) = D2mDfDe - D2mDf2;
        DDtwist.block<3, 3>( 4, 8 ) = D2mDeDf - D2mDf2;
        DDtwist.block<3, 3>( 8, 8 ) = D2mDf2;

        assert( isSymmetric( DDtwist ) );
    }

    setDependentsDirty();
}
