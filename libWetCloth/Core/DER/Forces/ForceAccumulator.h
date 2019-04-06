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

#ifndef FORCEACCUMULATOR_H_
#define FORCEACCUMULATOR_H_

#include "../Definitions.h"
#include "../StrandForce.h"

template<typename ForceT> class ForceAccumulator
{
public:
    static void accumulate( scalar& energy, const StrandForce& strand )
    {
        for( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            energy += ForceT::localEnergy( strand, vtx );
        }
    }

    static void accumulate( VecX& force, const StrandForce& strand )
    {
        typename ForceT::LocalForceType localF;
        for( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localF, strand, vtx );
            ForceT::addInPosition( force, vtx, localF );
        }
    }

    static void accumulateMultipliers( VecX& multipliers, const StrandForce& strand, const scalar& dt )
    {
        typename ForceT::LocalMultiplierType localL;
        for( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localL, strand, vtx, dt );
            ForceT::addInPosition( multipliers, vtx, localL );
        }        
    }

    // Jacobian of the Force <==>  - Hessian of the Energy
    static void accumulate( TripletXs& hessianOfEnergy, TripletXs& angularhessianOfEnergy, const StrandForce& strand )
    {
        typename ForceT::LocalJacobianType localJ;
        for( IndexType vtx = ForceT::s_first; vtx < strand.getNumVertices() - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localJ, strand, vtx );

            if( localJ.rows() > 6 ){ // (Bending & Twisting)
                for( IndexType r = 0; r < localJ.rows(); ++r )
                {
                    if(r % 4 == 3) {
                        for( IndexType c = 0; c < localJ.cols(); ++c )
                        {
                            if( c % 4 != 3 || isSmall( localJ(r,c) )  ) continue;
                            angularhessianOfEnergy.push_back( Triplets( (vtx - 1) * 4 + r, (vtx - 1) * 4 + c, localJ(r,c) ) );
                        }
                    } else {
                        for( IndexType c = 0; c < localJ.cols(); ++c )
                        {
                            if( c % 4 == 3 || isSmall( localJ(r,c) )  ) continue;
                            hessianOfEnergy.push_back( Triplets( (vtx - 1) * 4 + r, (vtx - 1) * 4 + c, localJ(r,c) ) );
                        }
                    }
                }
            }
            else{ // Stretch
                int trCount = 0;
                for( IndexType r = 0; r < localJ.rows(); ++r ){
                    if( r == 3 ){ // skip twist dof
                        ++trCount;
                    }                    
                    int tcCount = 0;
                    for( IndexType c = 0; c < localJ.cols(); ++c ){
                        if( c == 3 ){ // skip twist dof
                            ++tcCount;
                        }                        
                        if( isSmall( localJ(r,c) )  ) continue;
                        hessianOfEnergy.push_back( Triplets( vtx * 4 + r + trCount, vtx * 4 + c + tcCount, localJ(r,c) ) );
                    }
                }
            }
        }
    }
};

#endif
