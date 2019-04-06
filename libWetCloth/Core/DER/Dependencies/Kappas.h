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

#ifndef KAPPAS_H_
#define KAPPAS_H_

//#define USE_APPROX_GRAD_KAPPA

#include "MaterialFrames.h"

/**
 * Unit: no dimension.
 */
class Kappas: public DependencyNode<Vec2Array>
{
public:
    Kappas( CurvatureBinormals& curvatureBinormals, MaterialFrames<1>& materialFrames1,
            MaterialFrames<2>& materialFrames2 ) :
            DependencyNode<Vec2Array>( 1, curvatureBinormals.size() ), m_curvatureBinormals(
                    curvatureBinormals ), m_materialFrames1( materialFrames1 ), m_materialFrames2(
                    materialFrames2 )
    {
        m_curvatureBinormals.addDependent( this );
        m_materialFrames1.addDependent( this );
        m_materialFrames2.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Kappas";
    }

protected:
    virtual void compute();

    CurvatureBinormals& m_curvatureBinormals;
    MaterialFrames<1>& m_materialFrames1;
    MaterialFrames<2>& m_materialFrames2;
};

typedef Eigen::Matrix<scalar, 11, 2> GradKType;
typedef std::vector<GradKType, Eigen::aligned_allocator<GradKType>> GradKArrayType;

/**
 * Unit: cm^-1 for position derivatives, no dimension for theta derivatives
 */
class GradKappas: public DependencyNode<GradKArrayType>
{
public:
    GradKappas( Lengths& lengths, Tangents& tangents, CurvatureBinormals& curvatureBinormals,
            MaterialFrames<1>& materialFrames1, MaterialFrames<2>& materialFrames2, Kappas& kappas ) :
            DependencyNode<GradKArrayType>( 1, curvatureBinormals.size() ), //
            m_lengths( lengths ), //
            m_tangents( tangents ), //
            m_curvatureBinormals( curvatureBinormals ), //
            m_materialFrames1( materialFrames1 ), //
            m_materialFrames2( materialFrames2 ), //
            m_kappas( kappas )
    {
        m_lengths.addDependent( this );
        m_tangents.addDependent( this );
        m_curvatureBinormals.addDependent( this );
        m_materialFrames1.addDependent( this );
        m_materialFrames2.addDependent( this );
        m_kappas.addDependent( this );
    }

    virtual const char* name() const
    {
        return "GradKappas";
    }

protected:
    virtual void compute();

    Lengths& m_lengths;
    Tangents& m_tangents;
    CurvatureBinormals& m_curvatureBinormals;
    MaterialFrames<1>& m_materialFrames1;
    MaterialFrames<2>& m_materialFrames2;
    Kappas& m_kappas;
};

typedef Mat11Pair HessKType;
typedef std::vector<HessKType> HessKArrayType;

inline std::ostream& operator<<( std::ostream& os, const HessKType& HessKappa )
{
    os << "Hess kappa1: " << HessKappa.first << '\n';
    os << "Hess kappa2: " << HessKappa.second;

    return os;
}

/**
 * Unit: cm^-2 for position derivatives, no dimension for theta derivatives
 */
class HessKappas: public DependencyNode<HessKArrayType>
{
public:
    HessKappas( Lengths& lengths, Tangents& tangents, CurvatureBinormals& curvatureBinormals,
            MaterialFrames<1>& materialFrames1, MaterialFrames<2>& materialFrames2, Kappas& kappas ) :
            DependencyNode<HessKArrayType>( 1, curvatureBinormals.size() ), //
            m_lengths( lengths ), //
            m_tangents( tangents ), //
            m_curvatureBinormals( curvatureBinormals ), //
            m_materialFrames1( materialFrames1 ), //
            m_materialFrames2( materialFrames2 ), //
            m_kappas( kappas )
    {
        m_lengths.addDependent( this );
        m_tangents.addDependent( this );
        m_curvatureBinormals.addDependent( this );
        m_materialFrames1.addDependent( this );
        m_materialFrames2.addDependent( this );
        m_kappas.addDependent( this );
    }

    virtual const char* name() const
    {
        return "HessKappas";
    }

protected:
    virtual void compute();

    Lengths& m_lengths;
    Tangents& m_tangents;
    CurvatureBinormals& m_curvatureBinormals;
    MaterialFrames<1>& m_materialFrames1;
    MaterialFrames<2>& m_materialFrames2;
    Kappas& m_kappas;
};

typedef std::pair<Mat2, Mat2> ThetaHessKType;
typedef std::vector<ThetaHessKType> ThetaHessKArrayType;

inline std::ostream& operator<<( std::ostream& os, const ThetaHessKType& HessKappa )
{
    os << "ThetaHess kappa1: " << HessKappa.first << '\n';
    os << "ThetaHess kappa2: " << HessKappa.second;

    return os;
}

/**
 * Unit: no dimension
 */
class ThetaHessKappas: public DependencyNode<ThetaHessKArrayType>
{
public:
    ThetaHessKappas( CurvatureBinormals& curvatureBinormals, MaterialFrames<1>& materialFrames1,
            MaterialFrames<2>& materialFrames2 ) :
            DependencyNode<ThetaHessKArrayType>( 1, curvatureBinormals.size() ), //
            m_curvatureBinormals( curvatureBinormals ), //
            m_materialFrames1( materialFrames1 ), //
            m_materialFrames2( materialFrames2 )
    {
        m_curvatureBinormals.addDependent( this );
        m_materialFrames1.addDependent( this );
        m_materialFrames2.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ThetaHessKappas";
    }

protected:
    virtual void compute();

    CurvatureBinormals& m_curvatureBinormals;
    MaterialFrames<1>& m_materialFrames1;
    MaterialFrames<2>& m_materialFrames2;
};

#endif
