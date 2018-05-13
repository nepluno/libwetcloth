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

#ifndef DEGREESOFFREEDOM_H_
#define DEGREESOFFREEDOM_H_

#include "DependencyNode.h"

/**
 * Unit: cm for position dofs, no dimension for theta
 */
class DOFs: public DependencyNode<VecX>
{
public:
    DOFs( const VecX& dofValues ) :
            DependencyNode<VecX>( dofValues )
    {
        assert( dofValues.size() % 4 == 3 );
        m_numEdges = dofValues.size() / 4;
        setClean();
    }

    // As this class is meant to be the root of the dependency tree, we provide a const getter
    // This means that DOFs are never "computed" but updated by the outer loop algorithm
    const VecX& get() const
    {
        return m_value;
    }

    VecX& get() 
    {
        return m_value;
    }

    Vec3 getVertex( IndexType vtx ) const
    {
        assert( vtx < (m_numEdges + 1) );

        return get().segment<3>( 4 * vtx );
    }

    void setVertex( IndexType vtx, const Vec3& point )
    {
        m_value.segment<3>( 4 * vtx ) = point;
        setDependentsDirty();
    }

    void setDof( IndexType i, const scalar& val )
    {
        m_value[i] = val;
        setDependentsDirty();
    }

    // Accessors to the theta degrees of freedom
    const Eigen::Map<const VecX, Eigen::Unaligned, Eigen::InnerStride<4> > getThetas() const
    {
        return Eigen::Map<const VecX, Eigen::Unaligned, Eigen::InnerStride<4> >(
                m_value.data() + 3, m_numEdges );
    }

    scalar getTheta( IndexType vtx ) const
    {
        assert( vtx < m_numEdges );

        return get()[4 * vtx + 3];
    }

    void setThetas( const VecX& thetas, int numberOfFixedThetas = 0 )
    {
        assert( thetas.size()==m_numEdges );

        Eigen::Map<VecX, Eigen::Unaligned, Eigen::InnerStride<4> >(
                m_value.data() + 4 * numberOfFixedThetas + 3, m_numEdges - numberOfFixedThetas ) =
                thetas.tail( m_numEdges - numberOfFixedThetas );
        setDependentsDirty();
    }

    void setTheta( IndexType vtx, scalar theta )
    {
        m_value[4 * vtx + 3] = theta;
        setDependentsDirty();
    }

    IndexType getNumEdges() const
    {
        return m_numEdges;
    }

    IndexType getNumVertices() const
    {
        return m_numEdges + 1;
    }

    virtual const char* name() const
    {
        return "DOFs";
    }

protected:
    virtual void compute() // Not implemented as this is an pure input node
    {
        std::cerr << "DegreesOfFreedom::compute() should never be called" << std::endl;
    }

private:
    IndexType m_numEdges;
};

/**
 * Unit: cm
 */
class Edges: public DependencyNode<Vec3Array>
{
public:
    Edges( DOFs& dofs ) :
            DependencyNode<Vec3Array>( 0, dofs.getNumEdges() ), 
            m_dofs( dofs )
    {
        m_dofs.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Edges";
    }

protected:
    virtual void compute();

    DOFs& m_dofs;
};

/**
 * Unit: cm
 */
class Lengths: public DependencyNode<std::vector<scalar> >
{
public:
    Lengths( Edges& edges ) :
            DependencyNode<std::vector<scalar> >( 0, edges.size() ), 
            m_edges( edges )
    {
        m_edges.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Lengths";
    }

protected:
    virtual void compute();

    Edges& m_edges;
};

/**
 * Unit: no dimension
 */
class Tangents: public DependencyNode<Vec3Array>
{
public:
    Tangents( Edges& edges, Lengths& lengths ) :
            DependencyNode<Vec3Array>( 0, edges.size() ), 
            m_edges( edges ), 
            m_lengths( lengths )
    {
        m_edges.addDependent( this );
        m_lengths.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Tangents";
    }

protected:
    virtual void compute();

    Edges& m_edges;
    Lengths& m_lengths;
};

/**
 * Unit: no dimension
 */
class CurvatureBinormals: public DependencyNode<Vec3Array>
{
public:
    CurvatureBinormals( Tangents& tangents ) :
            DependencyNode<Vec3Array>( 1, tangents.size() ), 
            m_tangents( tangents )
    {
        m_tangents.addDependent( this );
    }

    virtual const char* name() const
    {
        return "CurvatureBinormals";
    }

protected:
    virtual void compute();

    Tangents& m_tangents;
};

/**
 * Unit: no dimension
 */
class TrigThetas: public DependencyNode<std::pair<VecX, VecX> >
{
public:
    TrigThetas( DOFs& dofs ) :
            DependencyNode<std::pair<VecX, VecX> >( std::make_pair( VecX(), VecX() ) ), 
            m_dofs( dofs )
    {
        m_dofs.addDependent( this );
    }

    virtual const char* name() const
    {
        return "TrigThetas";
    }

    const VecX& getSines()
    {
        return get().first;
    }

    const VecX& getCosines()
    {
        return get().second;
    }

protected:
    virtual void compute();

    DOFs& m_dofs;

private:
    void vdSinCos( const int n, const double a[], double r1[], double r2[] );

};

#endif
