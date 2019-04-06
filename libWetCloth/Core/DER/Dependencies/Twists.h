//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2012 Jean-Marie Aubry
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
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

#ifndef TWISTS_H_
#define TWISTS_H_

#include "ReferenceFrames.h"

/**
 * Unit: no dimension
 */
class Twists: public DependencyNode<std::vector<scalar> >
{
public:
    Twists( ReferenceTwists& refTwists, DOFs& dofs ) :
            DependencyNode<std::vector<scalar> >( 1, dofs.getNumEdges() ), m_refTwists( refTwists ), m_dofs(
                    dofs )
    {
        assert( size()==m_refTwists.size() );

        m_refTwists.addDependent( this );
        m_dofs.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Twists";
    }

protected:
    virtual void compute();

    ReferenceTwists& m_refTwists;
    DOFs& m_dofs;
};

/**
 * Unit: cm^-1
 */
class GradTwists: public DependencyNode<Vec11Array>
{
public:
    GradTwists( Lengths& lengths, CurvatureBinormals& curvatureBinormals ) :
            DependencyNode<Vec11Array>( 1, lengths.size() ), m_lengths( lengths ), m_curvatureBinormals(
                    curvatureBinormals )
    {
        m_lengths.addDependent( this );
        m_curvatureBinormals.addDependent( this );
    }

    virtual const char* name() const
    {
        return "GradTwists";
    }

protected:
    virtual void compute();

    Lengths& m_lengths;
    CurvatureBinormals& m_curvatureBinormals;
};

/**
 * Unit: cm^-2
 */
class GradTwistsSquared: public DependencyNode<Mat11Array>
{
public:
    GradTwistsSquared( GradTwists& gradTwists ) :
            DependencyNode<Mat11Array>( 1, gradTwists.size() ), m_gradTwists( gradTwists )
    {
        m_gradTwists.addDependent(this);
    }

    virtual const char* name() const
    {
        return "GradTwistsSquared";
    }

protected:
    virtual void compute();

    GradTwists& m_gradTwists;
};

/**
 * Unit: cm^-2
 */
class HessTwists: public DependencyNode<Mat11Array>
{
public:
    HessTwists( Tangents&tangents, Lengths& lengths, CurvatureBinormals& curvatureBinormals ) :
            DependencyNode<Mat11Array>( 1, lengths.size() ), m_tangents( tangents ), m_lengths(
                    lengths ), m_curvatureBinormals( curvatureBinormals )
    {
        m_tangents.addDependent( this );
        m_lengths.addDependent( this );
        m_curvatureBinormals.addDependent( this );
    }

    virtual const char* name() const
    {
        return "HessTwists";
    }

protected:
    virtual void compute();

    Tangents& m_tangents;
    Lengths& m_lengths;
    CurvatureBinormals& m_curvatureBinormals;
};

#endif
