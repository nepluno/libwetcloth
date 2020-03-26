//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ATTACH_FORCE_H
#define ATTACH_FORCE_H

#include <Eigen/Core>
#include "Force.h"
#include <iostream>
#include <memory>

class TwoDScene;

class AttachForce : public Force
{
public:

	AttachForce( const int pidx, const std::shared_ptr<TwoDScene>& scene, const scalar& k, const scalar& k_twist, const scalar& b, const scalar& m_b_twist );

	virtual ~AttachForce();

	virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E );

	virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE );

	virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt );

	virtual void addAngularHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt );

	virtual void updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt );

	virtual void preCompute();

	virtual void updateStartState();

	virtual Force* createNewCopy();

	virtual int numHessX();

	virtual int numAngularHessX();

	virtual int flag() const;

	int getParticleIndex() const;

	scalar getKs() const;
	scalar getKt() const;

private:
	int m_pidx;
	std::shared_ptr<TwoDScene> m_scene;
	scalar m_k;
	scalar m_k_twist;

	scalar m_b;
	scalar m_b_twist;

	Vector4s m_start_x;
};

#endif
