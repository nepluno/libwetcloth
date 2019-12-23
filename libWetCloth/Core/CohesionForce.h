//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

	virtual void updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt );
	
	virtual void preCompute();
	
	virtual void updateStartState();
	
	virtual Force* createNewCopy();
	
	virtual int numHessX();
	
	virtual int flag() const;

	virtual bool parallelized() const;
	
private:
    std::shared_ptr<TwoDScene> m_scene;
    std::vector< std::vector<int> > m_hess_offsets;
};

#endif
