//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
