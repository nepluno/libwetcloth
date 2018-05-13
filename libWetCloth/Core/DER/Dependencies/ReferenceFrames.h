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
#ifndef REFERENCEFRAMES_H_
#define REFERENCEFRAMES_H_

#include "DegreesOfFreedom.h"

/**
 * \brief This class maintains the first reference frame vectors, orthogonal to the edges.
 *
 * Initialized by space-parallel transportation along the rod, the reference frame vectors later
 * evolve by time-parallel transportation.
 *
 * Unit: no dimension
 */
class ReferenceFrames1: public DependencyNode<Vec3Array>
{
public:
    ReferenceFrames1( Tangents& tangents ) :
            DependencyNode<Vec3Array>( 0, tangents.size() ), m_tangents( tangents )
    {
        m_tangents.addDependent( this );

        storeInitialFrames( Vec3() );
    }

    virtual const char* name() const
    {
        return "ReferenceFrames1";
    }

    /**
     * \brief Propagates the initial (orthornormalized) reference frame initRefFrame1 by
     * space-parallel transportation along the rod
     * @param initRefFrame1 If zero or too far from orthogonal to the first edge, is replaced
     * by an arbitrary orthonormal vector.
     */
    void storeInitialFrames( const Vec3& initRefFrame1 );

    /**
     * \brief Access to m_previousTangents.
     *
     * m_previousTangents is an internal parameter caching the tangents from the previous
     * time step so we are able to do time-parallel transportation. The only need we have to
     * access it directly is for serialization/deserialization.
     */
    Vec3Array& getPreviousTangents()
    {
        return m_previousTangents;
    }

    void setPreviousTangents( const Vec3Array& deserialPrevTangents )
    {
        storeInitialFrames( Vec3() );
        m_previousTangents = deserialPrevTangents;
    }

    bool checkNormality();

protected:
    /**
     * \brief Computes new reference frames by time-parallel transportation along the
     * m_previousTangents->m_tangents motion.
     */
    virtual void compute();

    Tangents& m_tangents;
    Vec3Array m_previousTangents;
};

/**
 * \brief The second reference frame is simply maintained as tangent x referenceFrames1
 */
class ReferenceFrames2: public DependencyNode<Vec3Array>
{
public:
    ReferenceFrames2( Tangents& tangents, ReferenceFrames1& referenceFrames1 ) :
            DependencyNode<Vec3Array>( 0, tangents.size() ), m_tangents( tangents ), m_referenceFrames1(
                    referenceFrames1 )
    {
        m_tangents.addDependent( this );
        m_referenceFrames1.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ReferenceFrames2";
    }

protected:
    virtual void compute();

    Tangents& m_tangents;
    ReferenceFrames1& m_referenceFrames1;
};

/**
 * \brief This maintains the reference twist defined at the angle between the (space-parallel
 * transported) Bishop frame and the actual (time-parallel transported) reference frame.
 */
class ReferenceTwists: public DependencyNode<std::vector<scalar> >
{
public:
    ReferenceTwists( Tangents& tangents, ReferenceFrames1& referenceFrames1 ) :
            DependencyNode<std::vector<scalar> >( 1, tangents.size() ), m_tangents( tangents ), m_referenceFrames1(
                    referenceFrames1 )
    {
        m_tangents.addDependent( this );
        m_referenceFrames1.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ReferenceTwists";
    }

protected:
    virtual void compute();

    Tangents& m_tangents;
    ReferenceFrames1& m_referenceFrames1;
};

#endif
