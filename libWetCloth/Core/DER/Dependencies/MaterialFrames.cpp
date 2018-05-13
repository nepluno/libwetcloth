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

#include "MaterialFrames.h"

template<>
const char* MaterialFrames<1>::name() const
{
    return "MaterialFrames1";
}

template<>
const char* MaterialFrames<2>::name() const
{
    return "MaterialFrames2";
}

// GCC did not like this... split into two
//template<int FrameN>
//void MaterialFrames<FrameN>::compute()

template<>
inline Vec3 MaterialFrames<1>::linearMix( const Vec3& u, const Vec3& v, scalar s, scalar c )
{
    return c * u + s * v;
}

template<>
inline Vec3 MaterialFrames<2>::linearMix( const Vec3& u, const Vec3& v, scalar s, scalar c )
{
    return -s * u + c * v;
}

template<>
void MaterialFrames<1>::compute()
{
    m_value.resize( m_size );
    const Vec3Array& referenceFrames1 = m_referenceFrames1.get();
    const Vec3Array& referenceFrames2 = m_referenceFrames2.get();
    const VecX& sinThetas = m_trigThetas.getSines();
    const VecX& cosThetas = m_trigThetas.getCosines();

	for (IndexType vtx = 0; vtx < m_firstValidIndex; ++vtx)
	{
		m_value[vtx].setZero();
	}

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        const Vec3& u = referenceFrames1[vtx];
        const Vec3& v = referenceFrames2[vtx];
        const scalar s = sinThetas[vtx];
        const scalar c = cosThetas[vtx];

        m_value[vtx] = linearMix( u, v, s, c );
    }

    setDependentsDirty();
}

template<>
void MaterialFrames<2>::compute()
{
    m_value.resize( m_size );
    const Vec3Array& referenceFrames1 = m_referenceFrames1.get();
    const Vec3Array& referenceFrames2 = m_referenceFrames2.get();
    const VecX& sinThetas = m_trigThetas.getSines();
    const VecX& cosThetas = m_trigThetas.getCosines();

	for (IndexType vtx = 0; vtx < m_firstValidIndex; ++vtx)
	{
		m_value[vtx].setZero();
	}

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        const Vec3& u = referenceFrames1[vtx];
        const Vec3& v = referenceFrames2[vtx];
        const scalar s = sinThetas[vtx];
        const scalar c = cosThetas[vtx];

        m_value[vtx] = linearMix( u, v, s, c );
    }

    setDependentsDirty();
}
