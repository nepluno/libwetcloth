//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FORCEACCUMULATOR_H
#define FORCEACCUMULATOR_H

#include "../Definitions.h"
#include "../StrandForce.h"

template<typename ForceT> class ForceAccumulator
{
public:
    static void accumulate( scalar& energy, const StrandForce& strand )
    {
        for ( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            energy += ForceT::localEnergy( strand, vtx );
        }
    }

    static void accumulate( VecX& force, const StrandForce& strand )
    {
        typename ForceT::LocalForceType localF;
        for ( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localF, strand, vtx );
            ForceT::addInPosition( force, vtx, localF );
        }
    }

    static void accumulateMultipliers( VecX& multipliers, const StrandForce& strand, const scalar& dt )
    {
        typename ForceT::LocalMultiplierType localL;
        for ( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localL, strand, vtx, dt );
            ForceT::addInPosition( multipliers, vtx, localL );
        }
    }

    // Jacobian of the Force <==>  - Hessian of the Energy
    static void accumulate( TripletXs& hessianOfEnergy, TripletXs& angularhessianOfEnergy, const StrandForce& strand )
    {
        typename ForceT::LocalJacobianType localJ;
        for ( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localJ, strand, vtx );

            if ( localJ.rows() > 6 ) { // (Bending & Twisting)
                for ( IndexType r = 0; r < localJ.rows(); ++r )
                {
                    if (r % 4 == 3) {
                        for ( IndexType c = 0; c < localJ.cols(); ++c )
                        {
                            if ( c % 4 != 3 || isSmall( localJ(r, c) )  ) continue;
                            angularhessianOfEnergy.push_back( Triplets( (vtx - 1) * 4 + r, (vtx - 1) * 4 + c, localJ(r, c) ) );
                        }
                    } else {
                        for ( IndexType c = 0; c < localJ.cols(); ++c )
                        {
                            if ( c % 4 == 3 || isSmall( localJ(r, c) )  ) continue;
                            hessianOfEnergy.push_back( Triplets( (vtx - 1) * 4 + r, (vtx - 1) * 4 + c, localJ(r, c) ) );
                        }
                    }
                }
            }
            else { // Stretch
                int trCount = 0;
                for ( IndexType r = 0; r < localJ.rows(); ++r ) {
                    if ( r == 3 ) { // skip twist dof
                        ++trCount;
                    }
                    int tcCount = 0;
                    for ( IndexType c = 0; c < localJ.cols(); ++c ) {
                        if ( c == 3 ) { // skip twist dof
                            ++tcCount;
                        }
                        if ( isSmall( localJ(r, c) )  ) continue;
                        hessianOfEnergy.push_back( Triplets( vtx * 4 + r + trCount, vtx * 4 + c + tcCount, localJ(r, c) ) );
                    }
                }
            }
        }
    }
};

#endif
