//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
