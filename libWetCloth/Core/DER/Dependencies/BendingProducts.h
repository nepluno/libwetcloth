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

#ifndef BENDINGPRODUCTS_H_
#define BENDINGPRODUCTS_H_

#include "ElasticityParameters.h"
#include "Kappas.h"
#include "ElasticStrandUtils.h"

typedef std::vector<Mat2, Eigen::aligned_allocator<Mat2> > Mat2Array; ///< an array of 2d scalar matrices
typedef std::vector<Mat11, Eigen::aligned_allocator<Mat11> > Mat11Array; ///< an array of 11d scalar matrices

/**
 * \brief This class stores the products gradKappa^T B gradKappa that are used in both the viscous
 * and non-viscous bending forces.
 *
 * This product turned out to be one of the most costly operations in the collision-free simulation.
 *
 * The matrix B is not the same in the viscous and non-viscous cases, but differs only by a
 * proportionality factor that can be applied when computing the force Jacobian.
 *
 * Unit: cm^2
 */
class BendingProducts: public DependencyNode<Mat11Array>
{
public:
    BendingProducts( BendingMatrixBase& bendingMatrixBase, GradKappas& gradKappas ) :
            DependencyNode<Mat11Array>( 1, gradKappas.size() ), //
            m_bendingMatrixBase( bendingMatrixBase ), //
            m_gradKappas( gradKappas )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
        std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

        m_bendingMatrixBase.addDependent( this );
        m_gradKappas.addDependent( this );
    }

    virtual const char* name() const
    {
        return "BendingProducts";
    }

protected:
    virtual void compute()
    {
        m_value.resize( m_size );
        const MatX& bendingMatrix = m_bendingMatrixBase.get();
        const GradKArrayType& gradKappas = m_gradKappas.get();

        for( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
        {
            symBProduct<11>( m_value[vtx], bendingMatrix.block<2, 2>(vtx * 2, 0), gradKappas[vtx] );
        }

        setDependentsDirty();
    }

    BendingMatrixBase& m_bendingMatrixBase;
    GradKappas& m_gradKappas;
};

#endif
