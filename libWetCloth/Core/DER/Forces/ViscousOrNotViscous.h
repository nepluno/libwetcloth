//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef VISCOUS_OR_NOT_VISCOUS_H
#define VISCOUS_OR_NOT_VISCOUS_H

#include "../StrandForce.h"
#include "../StrandParameters.h"
// These classes are taken as template arguments for the internal forces,
// indicating whether we want the non-viscous or the viscous version.
// The forces call their ViscousT's static methods returning the appropriate
// stiffness and "rest shape" (the actual rest-shape for non-viscous or the
// shape at the beginning of time step for viscous).

class NonViscous
{
protected:
    NonViscous()
    {}

    virtual ~NonViscous()
    {}

public:
    static std::string getName()
    {
        return "";
    }

    static scalar bendingCoefficient( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->bendingCoefficient( );
    }

    static Mat2 bendingMatrix( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->bendingMatrix( vtx );
    }

    static const Vec2 kappaBar( const StrandForce& strand, int vtx )
    {
        return strand.m_restKappas[vtx];
    }

    static scalar kt( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->getKt( vtx );
    }

    static scalar thetaBar( const StrandForce& strand, int vtx )
    {
        return strand.m_restTwists[vtx];
    }

    static scalar ks( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->getKs( vtx );
    }

    static scalar ellBar( const StrandForce& strand, int vtx )
    {
        return strand.m_restLengths[vtx];
    }

    static const Vec2 bendingMultiplier( const StrandForce& strand, int vtx )
    {
        return strand.m_bending_multipliers.segment<2>( vtx * 2 );
    }

    static scalar stretchingMultiplier( const StrandForce& strand, int vtx )
    {
        return strand.m_stretching_multipliers( vtx );
    }

    static scalar twistingMultiplier( const StrandForce& strand, int vtx )
    {
        return strand.m_twisting_multipliers( vtx );
    }

    class NonDissipativeForce {};
};

class Viscous
{
protected:
    Viscous()
    {}

    virtual ~Viscous()
    {}

public:
    static std::string getName()
    {
        return "viscous ";
    }

    static scalar bendingCoefficient( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->viscousBendingCoefficient( );
    }

    static Mat2 bendingMatrix( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->viscousBendingMatrix( vtx );
    }

    static const Vec2 kappaBar( const StrandForce& strand, int vtx )
    {
        return strand.m_startState->m_kappas[vtx];
    }

    static scalar kt( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->getViscousKt( vtx );
    }

    static scalar thetaBar( const StrandForce& strand, int vtx )
    {
        return strand.m_startState->m_twists[vtx];
    }

    static scalar ks( const StrandForce& strand, int vtx )
    {
        return strand.m_strandParams->getViscousKs( vtx );
    }

    static scalar ellBar( const StrandForce& strand, int vtx )
    {
        return strand.m_startState->m_lengths[vtx];
    }

    static const Vec2 bendingMultiplier( const StrandForce& strand, int vtx )
    {
        return strand.m_viscous_bending_multipliers.segment<2>( vtx * 2 );
    }

    static scalar stretchingMultiplier( const StrandForce& strand, int vtx )
    {
        return strand.m_viscous_stretching_multipliers( vtx );
    }

    static scalar twistingMultiplier( const StrandForce& strand, int vtx )
    {
        return strand.m_viscous_twisting_multipliers( vtx );
    }

    class DissipativeForce {};
};

#endif
