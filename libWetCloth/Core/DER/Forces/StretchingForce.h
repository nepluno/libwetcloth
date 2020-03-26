//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef STRETCHINGFORCE_H
#define STRETCHINGFORCE_H

#include "ViscousOrNotViscous.h"
#include "../StrandForce.h"

template<typename ViscousT = NonViscous>
class StretchingForce
{
public:
    StretchingForce()
    {}

    virtual ~StretchingForce()
    {}

public:
    static const IndexType s_first = 0; // The first index on which this force can apply
    static const IndexType s_last = 1; // The last index (counting from the end)

    typedef Eigen::Matrix<scalar, 6, 1> LocalForceType;
    typedef Eigen::Matrix<scalar, 6, 6> LocalJacobianType;
    typedef VecX ForceVectorType;
    typedef scalar LocalMultiplierType;

    static std::string getName()
    {
        return ViscousT::getName() + "stretching";
    }

    static scalar localEnergy( const StrandForce& strand, const IndexType vtx );

    static void computeLocal( LocalMultiplierType& localL, const StrandForce& strand, const IndexType vtx, const scalar& dt );

    static void computeLocal( LocalForceType& localF, const StrandForce& strand, const IndexType vtx );

    static void computeLocal( LocalJacobianType& localJ, const StrandForce& strand,
                              const IndexType vtx );

    static void addInPosition( ForceVectorType& globalForce, const IndexType vtx,
                               const LocalForceType& localForce );

    static void addInPosition( VecX& globalMultiplier, const IndexType vtx, const LocalMultiplierType& localL );


    static void accumulateCurrentE( scalar& energy, StrandForce& strand );
    static void accumulateCurrentF( VecX& force, StrandForce& strand );
};

#endif
