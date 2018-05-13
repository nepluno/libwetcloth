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

#ifndef __COHESION_FORCE_H__
#define __COHESION_FORCE_H__

#include <Eigen/Core>
#include "Force.h"
#include <iostream>
#include <memory>

class TwoDScene;

class CohesionForce : public Force
{
public:
	
    CohesionForce(const std::shared_ptr<TwoDScene>& scene_ptr);
	
	virtual ~CohesionForce();
	
	virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E );
	
	virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE );
	
	virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt );
	
	virtual void preCompute();
	
	virtual void updateStartState();
	
	virtual Force* createNewCopy();
	
	virtual int numHessX();
	
	virtual int flag() const;
	
private:
    std::shared_ptr<TwoDScene> m_scene;
    std::vector< std::vector<int> > m_hess_offsets;
};

#endif
