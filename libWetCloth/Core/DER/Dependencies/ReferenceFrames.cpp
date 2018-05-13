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

#include "ReferenceFrames.h"
#include "ElasticStrandUtils.h"

void ReferenceFrames1::storeInitialFrames( const Vec3& initRefFrame1 )
{
    m_value.resize( m_size );
    const Vec3Array& tangents = m_tangents.get();
    const Vec3& tan0 = tangents[0];
    assert( isApproxUnit(tan0) );

    // Do we have an even approximately valid initial reference frame?
    if( initRefFrame1.squaredNorm() > 0.5 && fabs( initRefFrame1.dot( tan0 ) ) < 0.25 )
    {
        // If so, just project it on the plane normal to the tangent vector
        const Vec3 projectedInitRefFrame1 =
                ( initRefFrame1 - initRefFrame1.dot( tan0 ) * tan0 ).normalized();
        m_value[0] = projectedInitRefFrame1;
    }
    else // If a valid initial first reference frame hasn't been provided, use an arbitrary one
    {
        m_value[0] = findNormal( tan0 ).segment<3>(0);
    }

    // Next initial reference frames are obtained by space-parallel transportation along the rod
    for( IndexType vtx = 1; vtx < size(); ++vtx )
    {
        m_value[vtx] = orthonormalParallelTransport( m_value[vtx - 1], tangents[vtx - 1],
                tangents[vtx] );
        orthoNormalize( m_value[vtx], tangents[vtx] );
    }

    // Store tangents backup for time-parallel transport
    m_previousTangents = tangents;

    setClean();
    setDependentsDirty();
}

void ReferenceFrames1::compute()
{
    m_value.resize( m_size );
    const Vec3Array& tangents = m_tangents.get();

	for (IndexType vtx = 0; vtx < m_firstValidIndex; ++vtx)
	{
		m_value[vtx].setZero();
	}

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        Vec3& previousTangent = m_previousTangents[vtx];
        const Vec3& currentTangent = tangents[vtx];

        m_value[vtx] = orthonormalParallelTransport( m_value[vtx], previousTangent, currentTangent );
        orthoNormalize( m_value[vtx], currentTangent );

        // Store the current tangent for the next time-parallel transportation
        previousTangent = currentTangent;
    }
 
    setDependentsDirty();
}

bool ReferenceFrames1::checkNormality()
{
    bool normal = true;

    const Vec3Array& tangents = m_previousTangents;

    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        if( !isSmall( m_value[vtx].dot( tangents[vtx] ) ) )
        {
            normal = false;
            std::cerr << "ReferenceFrames1::checkNormality() fails at vtx = " << vtx << std::endl;
            break;
        }
    }    
    return normal;
}

void ReferenceFrames2::compute()
{
    m_value.resize( m_size );
    const Vec3Array& tangents = m_tangents.get();
    const Vec3Array& referenceFrames1 = m_referenceFrames1.get();
	for (IndexType vtx = 0; vtx < m_firstValidIndex; ++vtx)
	{
		m_value[vtx].setZero();
	}
    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        m_value[vtx] = tangents[vtx].cross( referenceFrames1[vtx] );
    }

    setDependentsDirty();
}

void ReferenceTwists::compute()
{
    m_value.resize( m_size );
    const Vec3Array& tangents = m_tangents.get();
    const Vec3Array& referenceFrames1 = m_referenceFrames1.get();
	for (IndexType vtx = 0; vtx < m_firstValidIndex; ++vtx)
	{
		m_value[vtx] = 0.0;
	}
    for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        const Vec3& u0 = referenceFrames1[vtx - 1];
        const Vec3& u1 = referenceFrames1[vtx];
        const Vec3& tangent = tangents[vtx];

        // transport reference frame to next edge
        Vec3 ut = orthonormalParallelTransport( u0, tangents[vtx - 1], tangent );

        // rotate by current value of reference twist
        const scalar beforeTwist = m_value[vtx];
        rotateAxisAngle( ut, tangent, beforeTwist );

        // compute increment to reference twist to align reference frames
        m_value[vtx] = beforeTwist + signedAngle( ut, u1, tangent );
    }

    setDependentsDirty();
}
