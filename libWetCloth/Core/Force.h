//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FORCE_H
#define FORCE_H

#include <Eigen/Core>

#include "MathDefs.h"

class TwoDScene;

class Force {
 public:
  virtual ~Force();

  virtual void addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                const VectorXs& m, const VectorXs& psi,
                                const scalar& lambda, scalar& E) = 0;

  virtual void addGradEToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, const VectorXs& psi,
                               const scalar& lambda, VectorXs& gradE) = 0;

  virtual void addLiquidGradEToNode(const TwoDScene& scene,
                                    std::vector<VectorXs>& node_rhs_x,
                                    std::vector<VectorXs>& node_rhs_y,
                                    std::vector<VectorXs>& node_rhs_z,
                                    const scalar& coeff);

  virtual void addHessXToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, const VectorXs& psi,
                               const scalar& lambda, TripletXs& hessE,
                               int hessE_index, const scalar& dt) = 0;

  virtual void addAngularHessXToTotal(const VectorXs& x, const VectorXs& v,
                                      const VectorXs& m, const VectorXs& psi,
                                      const scalar& lambda, TripletXs& hessE,
                                      int hessE_index, const scalar& dt);

  virtual int numHessX() = 0;

  virtual int numAngularHessX();

  virtual void preCompute() = 0;

  virtual void updateMultipliers(const VectorXs& x, const VectorXs& vplus,
                                 const VectorXs& m, const VectorXs& psi,
                                 const scalar& lambda, const scalar& dt) = 0;

  virtual void updateStartState() = 0;

  virtual Force* createNewCopy() = 0;

  virtual void postCompute(VectorXs& v, const scalar& dt);

  virtual int flag() const = 0;

  virtual bool parallelized() const;
};

#endif
